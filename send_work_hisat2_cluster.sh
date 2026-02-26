#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=8:00:00
#SBATCH --job-name=send_work
#SBATCH --output=/mnt/beegfs/cristobal/working_files/send_work.out
#SBATCH --error=/mnt/beegfs/cristobal/working_files/send_work.err

# ====================  COMMAND LINE  ====================
# sbatch ~/scripts/send_work_hisat2_cluster.sh <directory> <samples_names_file> <genome_file> <chrom_sizes_file>

# ====================  PARAMETERS  ====================

# Verify if there is any argument
if [ -z "$1" ]; then
  echo "Use: $0 <directory>"
  exit 1
fi

# Try to change working directoy
if cd "$1" 2>/dev/null; then
  # Print current working directory
  echo "Working directory: $(pwd)"
else
  # Deal with wrong path to working directory
  echo "Error: Impossible to change directory to '$1'"
  exit 1
fi

# Samples' names file.
SAMPLES_NAMES=$2
echo "Samples' names file is: '${SAMPLES_NAMES}'"

# Path to genome file.
GENOME_FILE=$3
echo "Genome file: ${GENOME_FILE}"

# Path to genome file.
CHROM_SIZES=$4
echo "Chromosome sizes file: ${CHROM_SIZES}"

# ====================  NAMES' ARRAY  ====================

# Names file must be a tabulated file with two columns with the original name in the first column and the final one in the second

# Create an array with samples' names from provided file 
NAMES=() 

# Read file and add to arrays
while IFS=$'\t' read -r col1 col2; do
	NAMES+=("$col2")
done < ${SAMPLES_NAMES}

# Check samples' names 
echo "Samples' names are:"

for i in ${NAMES[@]}; do 
	echo ${i}; 
done

# ====================  SEND TO WORK  ====================
for i in ${NAMES[@]}; do 
	sbatch ~/scripts/hisat2_cluster.sh $1 ./raw_data/fastq/concatenated/${i}_Fw.fastq.gz ./raw_data/fastq/concatenated/${i}_Rv.fastq.gz ${GENOME_FILE} ${CHROM_SIZES}; 
done

# Ending
echo "All work sent"






