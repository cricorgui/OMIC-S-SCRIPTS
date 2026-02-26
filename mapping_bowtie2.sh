#!/usr/bin/env bash
set -o errexit
set -o nounset
set -o pipefail

SCRIPT_NAME="$(basename "$0")"

uso() {
    cat <<EOF

Usage:
  $SCRIPT_NAME <directory> <fastq_fw> <fastq_rv> <genome_file> <chrom_sizes> <threads>

Description:
  Script for mapping paired samples with Bowtie2 filtering bam files by:
    - Properly pairs
    - Deleting secondary and suplementary alignments 
    - Reads with quality higher than 30
    - Deleting PCR duplicates

Also, index and stats from BAM files are created, as bedgraph and bigwig.

Extra directories will be created.

Mandatory arguments:
  directory      Directory where script is going to be executed
  fastq_fw       FASTQ forward (sample_Fw.fastq.gz)
  fastq_rv       FASTQ reverse (sample_Rv.fastq.gz)
  genome_file    Genome file built by Bowtie2 (usually more than one file are built, add only the name without variable extension)
  chrom_sizes    Chromosome sizes' file
  threads        Number of avaliable threads to run the tools

Options:
  -h, --help

Example:
  $SCRIPT_NAME results ./raw/A_Fw.fastq.gz ./raw/A_Rv.fastq.gz genome/hg38 genome/mm39.sizes 6


_______________________________________________________________________________________________________________

Script developed by Cristóbal Coronel Guisado (ORCID: https://orcid.org/0000-0003-4157-8231) in February 2026.

EOF
}

# ==================== OPTION -h / --help ====================

if [[ "${1-}" == "-h" || "${1-}" == "--help" ]]; then
    uso
    exit 0
fi

# ==================== GENERAL VALIDATION ====================

# Function for printing help in case of error
fatal() {
    echo "Error: $1" >&2
    echo >&2
    uso >&2           # help en stderr
    exit 1
}

# Verifying minimum number of arguments
if [[ $# -lt 6 ]]; then
    fatal "Missing mandatory arguments. Check script help:"
fi

DIRECTORY="$1"
FW_FILE="$2"
RV_FILE="$3"
GENOME_FILE="$4"
CHROM_SIZES="$5"
THREADS="$6"

# ==================== SPECIFIC VALIDATIONS ====================

[[ -d "$DIRECTORY" ]] || fatal "Working directory '$DIRECTORY' doesn't exit."
[[ -f "$FW_FILE"    ]] || fatal "Missing FASTQ forward file: $FW_FILE"
[[ -f "$RV_FILE"    ]] || fatal "Missing FASTQ reverse file: $RV_FILE"
[[ -f "$GENOME_FILE" ]] || fatal "Genome file not found: $GENOME_FILE"
[[ -f "$CHROM_SIZES" ]] || fatal "Chromosome sizes' file not found: $CHROM_SIZES"
[[ -f "$THREADS" ]] || fatal "Proper number of threads required: $THREADS"

# ==================== CHANGING DIRECTORY ====================

cd "$DIRECTORY" || fatal "Access to directory '$DIRECTORY' is not possible"
echo "Working directory: $(pwd)"

# ==================== SHOWING PARAMETERS ====================

echo "Forward FASTQ:      $FW_FILE"
echo "Reverse FASTQ:      $RV_FILE"
echo "Genome file:        $GENOME_FILE"
echo "Chromosome sizes:   $CHROM_SIZES"
echo "Number of threads:   $THREADS"

# ====================  CREATING DIRECTORIES  ====================

mkdir -p ./raw_data/sam
mkdir -p ./raw_data/bam/bam_stats/
mkdir -p ./data/bedgraph/separated/
mkdir -p ./data/bigwig/separated/

# ====================  MAPPING  ====================

bowtie2 -p ${THREADS} -x ${GENOME_FILE} -1 ${FW_FILE} -2 ${RV_FILE} --fr -k 1 -S ./raw_data/sam/$(basename ${FW_FILE} _Fw.fastq.gz).sam
echo "Mapping finished, starting conversion"

# ====================  FILTER SAM FILE AND CONVERSION TO BAM  ====================

# Filtering by properly pairs (-f 2), deleting secondary (-F 256) and suplementary (-F 2048) alignments and choosing only by quality reads (MAPQ>30)
samtools view -@ ${THREADS} -h -q 30 -f 2 -F 256 -F 2048 -b ./raw_data/sam/$(basename ${FW_FILE} _Fw.fastq.gz).sam | samtools sort -@ ${THREADS} -o ./raw_data/bam/$(basename ${FW_FILE} _Fw.fastq.gz)_filterQuality.bam
echo "SAM filtering, sorting and conversion to BAM conversion finished succesfully"

# Deleting PCR duplicates
sambamba markdup -r -t ${THREADS} -p ./raw_data/bam/$(basename ${FW_FILE} _Fw.fastq.gz)_filterQuality.bam ./raw_data/bam/$(basename ${FW_FILE} _Fw.fastq.gz)_filterQuality_dedup.bam
echo "Deleting PCR duplicates finished succesfully"

# Indexing BAM files
samtools index -@ ${THREADS} ./raw_data/bam/$(basename ${FW_FILE} _Fw.fastq.gz)_filterQuality_dedup.bam
echo "Indexing BAM finished succesfully"

# Stats from BAM files
samtools flagstat -@ ${THREADS} ./raw_data/bam/$(basename ${FW_FILE} _Fw.fastq.gz)_filterQuality_dedup.bam > ./raw_data/bam/bam_stats/$(basename ${FW_FILE} _Fw.fastq.gz)_filterQuality_dedup_bam_stats.txt
echo "BAM stats obtained successfully"

# Deleting unsorted bam
rm rm -rf ./raw_data/sam/$(basename ${FW_FILE} _Fw.fastq.gz).sam 
rm ./raw_data/bam/$(basename ${FW_FILE} _Fw.fastq.gz)_filterQuality.bam

# ====================  BEDGRAPH AND BIGWIG  ====================

# Bin size = 1bp and CPM normalization
bamCoverage -p ${THREADS} -b ./raw_data/bam/$(basename ${FW_FILE} _Fw.fastq.gz)_filterQuality_dedup.bam -of bedgraph --binSize 1 --normalizeUsing CPM -o ./data/bedgraph/separated/$(basename ${FW_FILE} _Fw.fastq.gz)_filterQuality_CPM_unsorted.bedgraph
bedtools sort -i ./data/bedgraph/separated/$(basename ${FW_FILE} _Fw.fastq.gz)_filterQuality_CPM_unsorted.bedgraph > ./data/bedgraph/separated/$(basename ${FW_FILE} _Fw.fastq.gz)_filterQuality_CPM.bedgraph
rm ./data/bedgraph/separated/$(basename ${FW_FILE} _Fw.fastq.gz)_filterQuality_CPM_unsorted.bedgraph
echo "Bedgraph file obtained successfully"

# Bigwig file
bamCoverage -p ${THREADS} -b ./raw_data/bam/$(basename ${FW_FILE} _Fw.fastq.gz)_filterQuality_dedup.bam -of bigwig --binSize 1 --normalizeUsing CPM -o ./data/bigwig/separated/$(basename ${FW_FILE} _Fw.fastq.gz)_filterQuality_CPM.bigwig
echo "Bigwig file obteined finished succesfully"

echo "Analysis of $(basename ${FW_FILE} _Fw.fastq.gz) sample finished"