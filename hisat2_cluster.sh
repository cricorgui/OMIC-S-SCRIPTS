#!/usr/bin/env bash
#SBATCH --cpus-per-task=20
#SBATCH --mem=40G
#SBATCH --time=8:00:00
#SBATCH --job-name=hisat2_map
#SBATCH --output=/mnt/beegfs/cristobal/working_files/hisat2_map_%j.out
#SBATCH --error=/mnt/beegfs/cristobal/working_files/hisat2_map_%j.err

set -o errexit
set -o nounset
set -o pipefail

SCRIPT_NAME="$(basename "$0")"

uso() {
    cat <<EOF

Usage:
  sbatch $SCRIPT_NAME <directory> <fastq_fw> <fastq_rv> <genome_file> <chrom_sizes> 

Description:
  SLURM script for mapping paired samples with HISAT2 using Dockers.
      - CPU per task: 20GB
      - Memory per task: 40GB

  NOTE: Remember to create all directories before running the code.

Mandatory arguments:
  directory      Directory where script is going to be executed
  fastq_fw       FASTQ forward (sample_Fw.fastq.gz)
  fastq_rv       FASTQ reverse (sample_Rv.fastq.gz)
  genome_file    Genome file built by Bowtie2 (usually more than one file are built, add only the name without variable extension)
  chrom_sizes    Chromosome sizes' file

Options:
  -h, --help

Example:
  sbatch $SCRIPT_NAME results ./raw/A_Fw.fastq.gz ./raw/A_Rv.fastq.gz genome/hg38 genome/hg38.sizes 


_______________________________________________________________________________________________________________

Script developed by Cristóbal Coronel Guisado (ORCID: https://orcid.org/0000-0003-4157-8231) in February 2026.


EOF
}

# ======== OPTION -h / --help ==========
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
if [[ $# -lt 5 ]]; then
    fatal "Missing mandatory arguments. Check script help:"
fi

DIRECTORY="$1"
FW_FILE="$2"
RV_FILE="$3"
GENOME_FILE="$4"
CHROM_SIZES="$5"

# ==================== SPECIFIC VALIDATIONS ====================

[[ -d "$DIRECTORY" ]] || fatal "Working directory '$DIRECTORY' doesn't exit."
[[ -f "$FW_FILE"    ]] || fatal "Missing FASTQ forward file: $FW_FILE"
[[ -f "$RV_FILE"    ]] || fatal "Missing FASTQ reverse file: $RV_FILE"
[[ -f "$GENOME_FILE" ]] || fatal "Genome file not found: $GENOME_FILE"
[[ -f "$CHROM_SIZES" ]] || fatal "Chromosome sizes' file not found: $CHROM_SIZES"

# ==================== CHANGING DIRECTORY ====================

cd "$DIRECTORY" || fatal "Access to directory '$DIRECTORY' is not possible"
echo "Working directory: $(pwd)"

# ==================== SHOWING PARAMETERS ====================

echo "Forward FASTQ:      $FW_FILE"
echo "Reverse FASTQ:      $RV_FILE"
echo "Genome file:        $GENOME_FILE"
echo "Chromosome sizes:   $CHROM_SIZES"

# ====================  MAPPING  ====================
docker run --rm -v /mnt/beegfs/cristobal/genomes/hg38/hg38_index_HISAT2:/genome -v ${DIRECTORY}/raw_data:/raw_data -v ${DIRECTORY}/raw_data/sam:/results nanozoo/hisat2:latest hisat2 -p 20 -x /genome/hg38_index_HISAT2 -1 /raw_data/fastq/concatenated/$(basename "$FW_FILE") -2 /raw_data/fastq/concatenated/$(basename "$RV_FILE") --rna-strandness RF -S /results/"$(basename ${FW_FILE} _Fw.fastq.gz).sam"

echo "Mapping finished, starting conversion"

# ====================  CONVERSION TO BAM  ====================
# Convert SAM a BAM	
docker run --rm --user root -v ${DIRECTORY}:/workspace biocontainers/samtools:v1.9-4-deb_cv1 samtools view -@ 20 -bS -O BAM "/workspace/raw_data/sam/$(basename ${FW_FILE} _Fw.fastq.gz).sam" -o "/workspace/raw_data/bam/$(basename ${FW_FILE} _Fw.fastq.gz)_unsorted.bam"
echo "Sam to bam conversion finished succesfully"

# Sort BAM files
docker run --rm --user root -v ${DIRECTORY}:/workspace biocontainers/samtools:v1.9-4-deb_cv1 samtools sort -@ 20 "/workspace/raw_data/bam/$(basename ${FW_FILE} _Fw.fastq.gz)_unsorted.bam" -o "/workspace/raw_data/bam/$(basename ${FW_FILE} _Fw.fastq.gz).bam"
echo "Bam sorting finished succesfully"

# Indexing BAM files
docker run --rm --user root -v ${DIRECTORY}:/workspace biocontainers/samtools:v1.9-4-deb_cv1 samtools index -@ 20 "/workspace/raw_data/bam/$(basename ${FW_FILE} _Fw.fastq.gz).bam"
echo "Indexing bam finished succesfully"

# Stats from BAM files
docker run --rm --user root -v ${DIRECTORY}:/workspace biocontainers/samtools:v1.9-4-deb_cv1 samtools flagstat -@ 20 "/workspace/raw_data/bam/$(basename ${FW_FILE} _Fw.fastq.gz).bam" | tee "${DIRECTORY}/raw_data/bam/bam_stats/$(basename ${FW_FILE} _Fw.fastq.gz)_bam_stats.txt"
echo "Bam stats obtained successfully"

# Multi bam stats
docker run --rm --user root -v ${DIRECTORY}:/workspace multiqc/multiqc "/workspace/raw_data/bam/bam_stats/"
echo "Multi bam stats performed successfully"

# Deleting unsorted bam
rm -f ${DIRECTORY}/raw_data/bam/"$(basename ${FW_FILE} _Fw.fastq.gz)_unsorted.bam"


# ====================  BEDGRAPH AND BIGWIG  ====================
# Bin size = 1bp and CPM normalization
docker run --rm --user root -v ${DIRECTORY}:/workspace ctomlins/bamcoverage bamCoverage -p 20 -b "/workspace/raw_data/bam/$(basename ${FW_FILE} _Fw.fastq.gz).bam" -of bedgraph --binSize 1 --normalizeUsing CPM --filterRNAstrand forward -o "/workspace/data/bedgraph/$(basename ${FW_FILE} .fastq.gz)_CPM_unsorted.bedgraph"
docker run --rm --user root -v ${DIRECTORY}:/workspace ctomlins/bamcoverage bamCoverage -p 20 -b "/workspace/raw_data/bam/$(basename ${RV_FILE} _Fw.fastq.gz).bam" -of bedgraph --binSize 1 --normalizeUsing CPM --filterRNAstrand reverse -o "/workspace/data/bedgraph/$(basename ${RV_FILE} .fastq.gz)_CPM_unsorted.bedgraph"

# Sort bedgraph
docker run --rm --user root -v ${DIRECTORY}:/workspace biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1 bedtools sort -i "/workspace/data/bedgraph/$(basename ${FW_FILE} .fastq.gz)_CPM_unsorted.bedgraph" | tee "${DIRECTORY}/data/bedgraph/$(basename ${FW_FILE} .fastq.gz)_CPM.bedgraph"
docker run --rm --user root -v ${DIRECTORY}:/workspace biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1 bedtools sort -i "/workspace/data/bedgraph/$(basename ${RV_FILE} .fastq.gz)_CPM_unsorted.bedgraph" | tee "${DIRECTORY}/data/bedgraph/$(basename ${RV_FILE} .fastq.gz)_CPM.bedgraph"
echo "Bedgraph file obtained successfully"

rm  -f ${DIRECTORY}/data/bedgraph/"$(basename ${FW_FILE} .fastq.gz)_CPM_unsorted.bedgraph"
rm  -f ${DIRECTORY}/data/bedgraph/"$(basename ${RV_FILE} .fastq.gz)_CPM_unsorted.bedgraph"

# Bigwig file
docker run --rm --user root -v ${DIRECTORY}:/workspace ctomlins/bamcoverage bamCoverage -p 20 -b "/workspace/raw_data/bam/$(basename ${FW_FILE} .fastq.gz).bam" -of bigwig --binSize 1 --normalizeUsing CPM --filterRNAstrand forward -o "/workspace/data/bigwig/$(basename ${FW_FILE} .fastq.gz)_CPM.bigwig"
docker run --rm --user root -v ${DIRECTORY}:/workspace ctomlins/bamcoverage bamCoverage -p 20 -b "/workspace/raw_data/bam/$(basename ${RV_FILE} .fastq.gz).bam" -of bigwig --binSize 1 --normalizeUsing CPM --filterRNAstrand reverse -o "/workspace/data/bigwig/$(basename ${RV_FILE} .fastq.gz)_CPM.bigwig"

echo "Bigwig file obteined finished succesfully"

echo "Analysis of $(basename ${FW_FILE} _Fw.fastq.gz) sample finished"










