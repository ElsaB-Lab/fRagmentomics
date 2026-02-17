#!/usr/bin/bash

# ==============================================================================
# Script Name: anonymize_snps.sh
# Description: Automates BAM anonymization by reverting SNPs to reference bases.
#              1. Downloads/Subsets gnomAD VCF for the target region.
#              2. Generates a list of SNPs.
#              3. Runs python script to sanitize BAM (and optionally shifts coords).
#              4. Indexes the output BAM.
#
# Usage:       bash anonymize_snps.sh --bam_in input.bam --bam_out output.bam ...
# ==============================================================================

# 1. Initialize variables with default values
SHIFT_COORD=false
BAM_IN=""
BAM_OUT=""
CHR=""
START=""
START_FASTA=""
END=""
PYTHON_SCRIPT="anonymize_snps.py" # Default assumes script is in current dir
CONDA_ENV="pysam"                 # Name of the conda environment to activate
PREVIEW_LIMIT=5                   # Number of reads to preview in logs

# 2. Parse arguments
while [[ "$#" -gt 0 ]]; do
  case $1 in
  --bam_in)
    BAM_IN="$2"
    shift 2
    ;;
  --bam_out)
    BAM_OUT="$2"
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
  --start_fasta)
    START_FASTA="$2"
    shift 2
    ;;
  --end)
    END="$2"
    shift 2
    ;;
  --python_script)
    PYTHON_SCRIPT="$2"
    shift 2
    ;;
  --conda_env)
    CONDA_ENV="$2"
    shift 2
    ;;
  --preview)
    PREVIEW_LIMIT="$2"
    shift 2
    ;;
  --shift_coord)
    SHIFT_COORD=true
    shift
    ;;
  *)
    echo "Unknown option: $1"
    echo "Usage: $0 --bam_in <file> --bam_out <file> --chr <chr> --start <int> --start_fasta <int> --end <int> [--shift_coord]"
    exit 1
    ;;
  esac
done

# 3. Validation
if [[ -z "$BAM_IN" || -z "$BAM_OUT" || -z "$CHR" || -z "$START" || -z "$END" ]]; then
  echo "Error: Missing required arguments."
  echo "Usage: $0 --bam_in <file> --bam_out <file> --chr <chr> --start <int> --start_fasta <int> --end <int> [--shift_coord]"
  exit 1
fi

echo "=========================================="
echo "Starting BAM Anonymization Pipeline"
echo "Region:      $CHR:$START-$END"
echo "Shift Coord: $SHIFT_COORD"
echo "Input:       $BAM_IN"
echo "Output:      $BAM_OUT"
echo "=========================================="

# 4. Activate Conda Environment
if [ -f "$(conda info --base)/etc/profile.d/conda.sh" ]; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$CONDA_ENV" || {
    echo "Error: Failed to activate conda env '$CONDA_ENV'"
    exit 1
  }
else
  echo "Warning: Conda not found at standard path. Assuming environment is already active."
fi

# 5. Prepare gnomAD VCF (Download & Subset)
QUERY_CHR=${CHR#chr} # Handle "chr" prefix mismatch

# Define filenames
RESOURCES="resources"
GNOMAD_VCF_NAME="gnomad.exomes.r2.1.1.sites.${QUERY_CHR}.vcf.bgz"
GNOMAD_VCF_URI="gs://gcp-public-data--gnomad/release/2.1.1/vcf/exomes/${GNOMAD_VCF_NAME}"
GNOMAD_VCF_FILE="${RESOURCES}/gnomad.exomes.r2.1.1.sites.${QUERY_CHR}.vcf.bgz"
GNOMAD_VCF_FILE_INDEX="${GNOMAD_VCF_FILE}.tbi"
REGION_VCF="${RESOURCES}/gnomad_exomes_${QUERY_CHR}_${START}-${END}.vcf.gz"
SITES_FILE="${RESOURCES}/gnomad_exomes_${QUERY_CHR}_${START}-${END}.snps.txt"

mkdir -p ${RESOURCES}

# Check and Download VCF
if [[ ! -f "$GNOMAD_VCF_FILE" ]]; then
  echo "[1/4] Downloading gnomAD VCF (Chromosome ${QUERY_CHR})..."
  gsutil cp "${GNOMAD_VCF_URI}" ${RESOURCES} || {
    echo "Error: gsutil failed to download VCF"
    exit 1
  }
else
  echo "[1/4] gnomAD VCF found locally."
fi

# Check and Download Index (.tbi)
if [[ ! -f "${GNOMAD_VCF_FILE_INDEX}" ]]; then
  echo "[1/4] Downloading gnomAD VCF Index..."
  gsutil cp "${GNOMAD_VCF_URI}.tbi" ${RESOURCES} || {
    echo "Error: gsutil failed to download Index"
    exit 1
  }
fi

# 6. Generate Sites File (bcftools)
echo "[2/4] Generating SNP list for region..."

if [[ ! -f "$REGION_VCF" ]]; then
  bcftools view \
    -r ${QUERY_CHR}:${START}-${END} \
    "${GNOMAD_VCF_FILE}" \
    -O z -o "$REGION_VCF"
fi

bcftools query \
  -i 'TYPE="snp"' \
  -f '%CHROM\t%POS\t%REF\t%ALT\n' \
  -r ${QUERY_CHR}:${START}-${END} \
  "${GNOMAD_VCF_FILE}" \
  >"${SITES_FILE}"

NUM_SITES=$(wc -l <"${SITES_FILE}")
echo "      Found ${NUM_SITES} SNPs in region."

# 7. Run Python Anonymization Script
echo "[3/4] Running Python Anonymization (Shift: $SHIFT_COORD)..."

if [[ ! -f "$PYTHON_SCRIPT" ]]; then
  echo "Error: Python script '$PYTHON_SCRIPT' not found."
  exit 1
fi

# Construct command args
PY_ARGS="--input_bam ${BAM_IN} --output_bam ${BAM_OUT} --sites_file ${SITES_FILE} --preview ${PREVIEW_LIMIT}"

# If shift coord is requested, pass the extra arguments to Python
if [ "$SHIFT_COORD" == true ]; then
  PY_ARGS="$PY_ARGS --shift_coord --chrom ${CHR} --start_fasta ${START_FASTA}"
fi

python "${PYTHON_SCRIPT}" $PY_ARGS

if [[ $? -ne 0 ]]; then
  echo "Error: Python script failed."
  exit 1
fi

# 8. Index the Output BAM
echo "[4/4] Indexing output BAM..."
samtools index "${BAM_OUT}"

echo "=========================================="
echo "Anonymization Complete."
echo "Result: ${BAM_OUT}"
echo "Index:  ${BAM_OUT}.bai"
echo "=========================================="
