#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --job-name="bam-sort-index"
#SBATCH --output=bam-sort-index.out
#SBATCH --mail-user=emilyise@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc


###############################################################################
############ SET UP ############
# set up temporary environment
echo "Establishing temporary environment..."
python3 -m venv .tempenv
source .tempenv/bin/activate
echo "Temporary environment activated."

# load packages
echo "Checking for required modules." 
## define paths
MYLIB_DIR="./mylib"
BISMARK_DIR="${MYLIB_DIR}/Bismark"
BISMARK_URL="https://github.com/FelixKrueger/Bismark.git"

## check if the Bismark directory exists
if [ -d "$BISMARK_DIR" ]; then
    echo "Bismark directory already exists. Proceeding to the next step."
else
    echo "Bismark directory not found. Downloading it now..."
    
    ### create the mylib directory if it doesn't exist
    mkdir -p "$MYLIB_DIR"
    cd "$MYLIB_DIR"

    ### download and extract the Bismark directory
    git clone "$BISMARK_URL"
    export PATH=$PWD/Bismark:$PATH

    echo "Bismark directory downloaded and extracted successfully."
fi

## load bowtie and samtools
echo "Loading required modules."
module load gcc/11.2.0
module load bowtie2/2.4.4
module load samtools/1.16.1
echo "Modules loaded."

# set project directory
PRJ_DIR="/projects/rpci/joyceohm/Emily/RRBS-pipeline/Filtered_Human/"
cd $PRJ_DIR 

###############################################################################
############ PULL BAM FILES ############
# pull bam files
bam=$(find . -name '*.bam')

# Print the found file paths (for debugging)
## echo "$bam"

############ SORT ############
# sort
printf "%s\n" "$bam" | \
  parallel -j 6 'outname=$(echo {} | sed "s/.bam/sorted.bam/"); \
  samtools sort -o "$outname" {}'

############ INDEX ############
# index
## pull sorted bam files
sbam=$(find . -name '*sorted.bam')

## run indexing 
printf "%s\n" "$sbam" | \
  parallel -j 16 'samtools index {}'
