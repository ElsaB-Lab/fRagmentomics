#!/usr/bin/bash

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --bam_file)
            BAM_FILE="$2"
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
        --bam_out)
            BAM_OUT="$2"
            shift 2
            ;;
        --fasta_out)
            FASTA_OUT="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done


LN=$((END - START + 1))  # Compute new sequence length

TMP_DIR="bam/${SAMPLE_ANON}_${CHR}_${START}_${END}_tmp"
BAM_FILE_SUB_QNAME_UNIQ="${TMP_DIR}/${SAMPLE_ANON}_${CHR}_${START}_${END}.uniq.qname.txt"
BAM_FILE_SUB_QNAME_MAP="${TMP_DIR}/${SAMPLE_ANON}_${CHR}_${START}_${END}.map.qname.txt"
BAM_FILE_SUB_TRUE="${TMP_DIR}/${SAMPLE_TRUE}_${CHR}_${START}_${END}.true.bam"
BAM_FILE_SUB_TRUE_RAND="${TMP_DIR}/${SAMPLE_TRUE}_${CHR}_${START}_${END}.true.rand.bam"
BAM_FILE_SUB_ANON_RG="${TMP_DIR}/${SAMPLE_ANON}_${CHR}_${START}_${END}.anonymized.rg.bam"

# create temp dir
mkdir -p ${TMP_DIR}

# subset BAM file
module load samtools
echo -n "PROCESS: selecting reads covering ${CHR}:${START}-${END} ... "
samtools view -h ${BAM_FILE} ${CHR}:${START}-${END} > ${BAM_FILE_SUB_TRUE}
echo "done!"

# subset BAM file
READ_COUNT=$(samtools view -c ${BAM_FILE_SUB_TRUE})
TARGET_READS=10000

if [[ $READ_COUNT -le $TARGET_READS ]]; then
        FRACTION=1.0  # Keep all reads if there are fewer than 10,000
    else
	FRACTION=$(awk -v total="$READ_COUNT" -v target="$TARGET_READS" 'BEGIN { printf "%.6f", target / total }')
fi

echo -n "PROCESS: selecting $FRACTION of all reads at random ... "
samtools view -h -s $FRACTION -b ${BAM_FILE_SUB_TRUE} > ${BAM_FILE_SUB_TRUE_RAND}
echo "done!"

# Count the number of reads in the BAM file
READ_COUNT=$(samtools view -c ${BAM_FILE_SUB_TRUE_RAND})

# Get the file size of the BAM file in a human-readable format
FILE_SIZE=$(stat -c %s ${BAM_FILE_SUB_TRUE_RAND})
HUMAN_READABLE_SIZE=$(numfmt --to=iec-i --suffix=B ${FILE_SIZE})

# Print the results
echo "INFO: The test BAM file contains ${READ_COUNT} reads."
echo "INFO: The test BAM file size is ${HUMAN_READABLE_SIZE}."

# anonymize sample id
echo "PROCESS: Anonymizing read groups info using Picard AddOrReplaceReadGroups..."
module load singularity
PICARD_SINGULARITY="/mnt/beegfs02/scratch/y_pradat/resources/singularity/picard_3.2.0--hdfd78af_0.simg"
singularity exec --bind "$PWD" ${PICARD_SINGULARITY} picard AddOrReplaceReadGroups \
  I=${BAM_FILE_SUB_TRUE_RAND} \
  O=${BAM_FILE_SUB_ANON_RG} \
  RGID=1 \
  RGLB=lib1 \
  RGPL=illumina \
  RGSM=${SAMPLE_ANON} \
  RGPU=unit1

echo "INFO: Checking anonymization of read groups info..."
echo "INFO: RG infos before anonymization."
samtools view -H ${BAM_FILE_SUB_TRUE_RAND} | grep '^@RG'

echo "INFO: RG infos after anonymization."
samtools view -H ${BAM_FILE_SUB_ANON_RG} | grep '^@RG'

# Extract unique QNAME
samtools view "${BAM_FILE_SUB_ANON_RG}" | cut -f1 | sort -u > ${BAM_FILE_SUB_QNAME_UNIQ}

# Count unique QNAMEs and determine padding size
echo "PROCESS: Anonymizing read names..."
TOTAL_QNAMES=$(wc -l < ${BAM_FILE_SUB_QNAME_UNIQ})
ID_WIDTH=$(( ${#TOTAL_QNAMES} + 1 ))

# Generate anonymized QNAME mapping
awk -v width="$ID_WIDTH" '{printf "%s\tREAD_%0*d\n", $0, width, NR}' ${BAM_FILE_SUB_QNAME_UNIQ} > "${BAM_FILE_SUB_QNAME_MAP}"

# Replace QNAMEs in BAM and ajsut POS/PNEXT where necessary
samtools view -h "${BAM_FILE_SUB_ANON_RG}" | 
awk -v OFS="\t" -v START="$START" -v CHR="$CHR" -v LN="$LN" -v mapfile="${BAM_FILE_SUB_QNAME_MAP}" '
    BEGIN {
	while ((getline < mapfile) > 0) { qmap[$1] = $2 }
    }
    /^@PG/ { next }  # Skip @PG lines
    /^@/ { 
	if ($1 == "@SQ" && $2 == "SN:" CHR) { 
	    print "@SQ\tSN:" CHR "\tLN:" LN; 
	} else {
	    if ($1 != "@SQ") {
		print;
	    }
	}
	next;
    }
    { 
	$1 = qmap[$1]; 

	# Adjust column 4 (POS)
	new_pos = $4 - (START - 1);

        # Filter based on column 4: Keep only if new_pos is in range [1, LN]
        if (new_pos < 1 || new_pos > LN) next;
	
	# Update POS
	$4 = new_pos;

	# Adjust column 8 (PNEXT) only if column 7 (RNEXT) is "="
	if ($7 == "=" && $8 != 0 && $8 ~ /^[0-9]+$/) {
	    new_pnext = $8 - (START - 1);

	    # Ensure new_pnext is within [1, LN]
	    if (new_pnext < 1 || new_pnext > LN) next;

	   # Update PNEXT
	   $8 = new_pnext;
	}

	print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11;  # Keep only first 11 fields
    }
' | samtools view -b -o "${BAM_OUT}"

# Add index
samtools index ${BAM_OUT}
echo "Anonymized test BAM written to ${BAM_OUT}"

# create small fasta
samtools faidx ${REF_FASTA} "${CHR}:${START}-${END}" > ${FASTA_OUT}
sed -i "s/>${CHR}:${START}-${END}/>${CHR}/" ${FASTA_OUT}
samtools faidx ${FASTA_OUT}
echo "Test FASTA written to ${FASTA_OUT}"

# Clean
rm -rf ${TMP_DIR}
