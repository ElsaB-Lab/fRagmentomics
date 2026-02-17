#!/usr/bin/bash

# Initialize variables
BAM_IN=""
REF_FASTA=""
SAMPLE_TRUE=""
SAMPLE_ANON=""
TARGET_READS=-1
CHR=""
START=""
END=""
START_FASTA=""
END_FASTA=""
BAM_OUT_SHIFTED=""
FASTA_OUT_SHIFTED=""
QNAME_MAP=""

while [[ "$#" -gt 0 ]]; do
  case $1 in
  --bam_in)
    BAM_IN="$2"
    shift 2
    ;;
  --ref_fasta)
    REF_FASTA="$2"
    shift 2
    ;;
  --sample_true)
    SAMPLE_TRUE="$2"
    shift 2
    ;;
  --sample_anon)
    SAMPLE_ANON="$2"
    shift 2
    ;;
  --target_reads)
    TARGET_READS="$2"
    shift 2
    ;;
  --chr)
    CHR="$2"
    shift 2
    ;;
  --start)
    START="$2"
    shift 2
    ;;
  --end)
    END="$2"
    shift 2
    ;;
  --start_fasta)
    START_FASTA="$2"
    shift 2
    ;;
  --end_fasta)
    END_FASTA="$2"
    shift 2
    ;;
  --bam_out_shifted)
    BAM_OUT_SHIFTED="$2"
    shift 2
    ;;
  --fasta_out_shifted)
    FASTA_OUT_SHIFTED="$2"
    shift 2
    ;;
  --qname_map)
    QNAME_MAP="$2"
    shift 2
    ;;
  *)
    echo "Unknown option: $1"
    exit 1
    ;;
  esac
done

TMP_DIR="bam/${SAMPLE_ANON}_${CHR}_${START}_${END}_tmp"
mkdir -p ${TMP_DIR}
trap 'rm -rf "$TMP_DIR"' EXIT

BAM_IN_SUB_TRUE="${TMP_DIR}/true.bam"
BAM_IN_SUB_TRUE_RAND="${TMP_DIR}/true.rand.bam"
BAM_IN_SUB_ANON_RG="${TMP_DIR}/anon.rg.bam"
BAM_IN_SUB_QNAME_UNIQ="${TMP_DIR}/qname.uniq.txt"
BAM_IN_SUB_QNAME_MAP="${TMP_DIR}/qname.map.txt"

# 1. Select reads covering the WIDER region (START to END)
echo -n "PROCESS: selecting reads covering ${CHR}:${START}-${END} ... "
samtools view -h ${BAM_IN} ${CHR}:${START}-${END} >${BAM_IN_SUB_TRUE}
echo "done!"

# 2. Subsample
READ_COUNT=$(samtools view -c ${BAM_IN_SUB_TRUE})
if [[ $TARGET_READS == -1 ]] || [[ $READ_COUNT -le $TARGET_READS ]]; then
  FRACTION=1.0
else
  FRACTION=$(awk -v total="$READ_COUNT" -v target="$TARGET_READS" 'BEGIN { printf "%.6f", target / total }')
fi
samtools view -h -s $FRACTION -b ${BAM_IN_SUB_TRUE} >${BAM_IN_SUB_TRUE_RAND}

# 3. Anonymize RG
module load singularity
PICARD_SINGULARITY="/mnt/beegfs01/scratch/y_pradat/resources/singularity/picard_3.2.0--hdfd78af_0.simg"
singularity exec --bind "$PWD" ${PICARD_SINGULARITY} picard AddOrReplaceReadGroups \
  I=${BAM_IN_SUB_TRUE_RAND} O=${BAM_IN_SUB_ANON_RG} \
  RGID=1 RGLB=lib1 RGPL=illumina RGSM=${SAMPLE_ANON} RGPU=unit1

# 4. Anonymize QNAMEs
samtools view "${BAM_IN_SUB_ANON_RG}" | cut -f1 | sort -u >${BAM_IN_SUB_QNAME_UNIQ}
TOTAL_QNAMES=$(wc -l <${BAM_IN_SUB_QNAME_UNIQ})
ID_WIDTH=$((${#TOTAL_QNAMES} + 1))
awk -v width="$ID_WIDTH" '{printf "%s\tREAD_%0*d\n", $0, width, NR}' ${BAM_IN_SUB_QNAME_UNIQ} >"${BAM_IN_SUB_QNAME_MAP}"
cp ${BAM_IN_SUB_QNAME_MAP} ${QNAME_MAP}

################################################################
# OUTPUT: SHIFTED COORDINATES
################################################################
echo "PROCESS: generating SHIFTED output (Loose Filter)..."
LN=$((END_FASTA - START_FASTA + 1))

samtools view -h "${BAM_IN_SUB_ANON_RG}" |
  awk -v OFS="\t" -v S_FASTA="$START_FASTA" -v E_FASTA="$END_FASTA" -v CHR="$CHR" -v LN="$LN" -v mapfile="${BAM_IN_SUB_QNAME_MAP}" '
  BEGIN {
      while ((getline < mapfile) > 0) qmap[$1] = $2
  }
  /^@PG/ { next }
  /^@/ {
      if ($1 == "@SQ" && $2 == "SN:" CHR) {
        print "@SQ\tSN:" CHR "\tLN:" LN;
      } else print;
      next;
  }
  {
    # 1. Filter: Check Read Position (Must be inside [START_FASTA, END_FASTA])
    if ($4 < S_FASTA || $4 > E_FASTA) next;

    # Apply Anonymized Name
    $1 = qmap[$1];

    # SHIFT COORDINATES (By START_FASTA)
    shift_val = S_FASTA - 1;

    # Shift Read Position
    $4 = $4 - shift_val;

    # Shift Mate Position ONLY if mate is on the same chromosome
    # $7 is RNEXT (Reference name of the mate/next read)
    # It can be "=" (same as current) or the actual name (e.g. "chr7")
    if ($7 == "=" || $7 == CHR) {
        $8 = $8 - shift_val;
    }

    print $0;
  }
' | samtools view -b -o "${BAM_OUT_SHIFTED}"

samtools index ${BAM_OUT_SHIFTED}

################################################################
# FASTA GENERATION
################################################################

#  Shifted Fasta (Using START_FASTA)
if [ -n "$FASTA_OUT_SHIFTED" ]; then
  samtools faidx ${REF_FASTA} "${CHR}:${START_FASTA}-${END_FASTA}" >${FASTA_OUT_SHIFTED}
  sed -i "s/>${CHR}:${START_FASTA}-${END_FASTA}/>${CHR}/" ${FASTA_OUT_SHIFTED}
  samtools faidx ${FASTA_OUT_SHIFTED}
fi

if [ -d "$TMP_DIR" ]; then rm -rf "$TMP_DIR"; fi
rmdir bam 2>/dev/null
