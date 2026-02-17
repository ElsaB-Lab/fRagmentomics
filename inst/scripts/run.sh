#!/usr/bin/bash

format_number() {
  if [[ "$1" == "-1" ]]; then
    echo ""
    return 0
  fi
  if [[ -z "$1" || ! "$1" =~ ^[0-9]+$ ]]; then
    echo "Error: A valid integer argument is required." >&2
    return 1
  fi
  local num=$1
  if ((num >= 1000000000)); then
    echo "_$((num / 1000000000))B"
  elif ((num >= 1000000)); then
    echo "_$((num / 1000000))M"
  elif ((num >= 1000)); then
    echo "_$((num / 1000))k"
  else echo "_$num"; fi
}

REF_FASTA="/mnt/beegfs02/database/bioinfo/Index_DB/BWA/0.7.17/UCSC/hg19.fa"
REF_FASTA_BASENAME=$(basename "${REF_FASTA}" | sed -E 's/\.(fa|fasta)$//')

OUT_DIR="../extdata"
mkdir -p $OUT_DIR/fasta $OUT_DIR/bam $OUT_DIR/mutations

DATA_DIR="/path/to/data/dir/with/bam"
SAMPLE_TRUE="sample"
SAMPLE_ANON="cfdna-egfr-del"

# Coordinates
CHR="chr7"
POS=55242464
REF="AGGAATTAAGAGAAGC"
ALT="A"

# Region of Interest (ROI)
START=$(($POS - 600))
END=$(($POS + 600))

# FASTA / BAM Boundaries (ROI + Padding)
PADDING=10000
START_FASTA=$((START - PADDING))
END_FASTA=$((END + PADDING))

# 1. Shifted (Mini) FASTA
FASTA_OUT_SHIFTED="$OUT_DIR/fasta/${REF_FASTA_BASENAME}_${CHR}_${START_FASTA}_${END_FASTA}.fa"

BAM_IN="${DATA_DIR}/${SAMPLE_TRUE}.bam"

# The target read counts used for the benchmarking of performances
# READ_COUNTS=(100 500 1000 2000 5000 10000 20000 50000 100000 -1)

# The target read counts of the test data provided
READ_COUNTS=(10000)

for TARGET_READS in "${READ_COUNTS[@]}"; do
  TARGET_READS_STR=$(format_number "$TARGET_READS")
  echo "Processing ${SAMPLE_ANON} with target reads: ${TARGET_READS} (${TARGET_READS_STR})"

  MUTATIONS_OUT_SHIFTED="$OUT_DIR/mutations/${SAMPLE_ANON}_${CHR}_${START}_${END}${TARGET_READS_STR}.mutations.tsv"
  BAM_TMP_SHIFTED="${OUT_DIR}/bam/${SAMPLE_ANON}_${CHR}_${START}_${END}${TARGET_READS_STR}.tmp.bam"
  BAM_OUT_SHIFTED="${OUT_DIR}/bam/${SAMPLE_ANON}_${CHR}_${START}_${END}${TARGET_READS_STR}.bam"
  QNAME_MAP="${OUT_DIR}/bam/${SAMPLE_ANON}_${CHR}_${START}_${END}${TARGET_READS_STR}_qname_map.txt"

  # 1. Subsample and extract region
  bash subsample_and_anonymize_reads.sh \
    --bam_in "$BAM_IN" \
    --ref_fasta "$REF_FASTA" \
    --sample_true "$SAMPLE_TRUE" \
    --sample_anon "$SAMPLE_ANON" \
    --chr "$CHR" \
    --start "$START" \
    --end "$END" \
    --start_fasta "$START_FASTA" \
    --end_fasta "$END_FASTA" \
    --bam_out_shifted "$BAM_TMP_SHIFTED" \
    --fasta_out_shifted "$FASTA_OUT_SHIFTED" \
    --qname_map "$QNAME_MAP" \
    --target_reads "$TARGET_READS"

  # 2. Anonymize SNPs (Shifted)
  bash anonymize_snps.sh \
    --bam_in "$BAM_TMP_SHIFTED" \
    --chr "$CHR" \
    --start "$START" \
    --start_fasta "$START_FASTA" \
    --end "$END" \
    --bam_out "$BAM_OUT_SHIFTED" \
    --preview 5 \
    --shift_coord

  # 3. Mutations Files
  # Shifted: Calculate new position relative to START_FASTA
  SHIFT_OFFSET=$((START_FASTA - 1))
  echo -e "CHROM\tPOS\tREF\tALT" >"${MUTATIONS_OUT_SHIFTED}"
  echo -e "${CHR}\t$((POS - SHIFT_OFFSET))\t${REF}\t${ALT}" >>"${MUTATIONS_OUT_SHIFTED}"

  # Cleanup
  if [ -f "$BAM_TMP_SHIFTED" ]; then rm "$BAM_TMP_SHIFTED"; fi

  echo "---------------------------------------------------"
done
