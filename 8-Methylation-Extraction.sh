#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --mem=50G
#SBATCH --job-name="methylation-_extraction"
#SBATCH --output=methylation_extraction.out
#SBATCH --mail-user=emilyise@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

###############################################################################
############ SET UP TEMPORARY ENVIRONMENT ############
echo "Establishing temporary environment..."

# set up temporary environment
python3 -m venv .tempenv3
source .tempenv3/bin/activate

echo "Temporary environment activated successfully."

# load bowtie and samtools
echo "Loading required modules."

module load gcc/11.2.0
module load bowtie2/2.4.4
module load samtools/1.16.1
echo "Modules loaded."

###############################################################################
############ LOAD PACKAGES ############
echo "Checking for required modules..."

# Define paths
PROJECT_DIR="/projects/rpci/joyceohm/Emily/RRBS-pipeline"
cd $PROJECT_DIR
MYLIB_DIR="./mylib"
BISMARK_DIR="${MYLIB_DIR}/Bismark"
BISMARK_URL="https://github.com/FelixKrueger/Bismark.git"
PYTHON2_URL="https://www.python.org/ftp/python/2.7.18/Python-2.7.18.tgz"
PYTHON2_DIR="${MYLIB_DIR}/python2"

# Create the mylib directory if it doesn't exist
mkdir -p "$MYLIB_DIR"


############ BISMARK SETUP ############
# Check if the Bismark directory exists
cd $PROJECT_DIR
cd $MYLIB_DIR

if [ -d "$BISMARK_DIR" ]; then
    echo "Bismark directory already exists. Proceeding to the next step."
else
    echo "Bismark directory not found. Downloading it now..."

    # Clone the Bismark repository
    git clone "$BISMARK_URL" "$BISMARK_DIR"

    # Add Bismark to PATH within the temporary environment
    export PATH=$PWD/Bismark:$PATH

    echo "Bismark directory downloaded and extracted successfully."
fi

export PATH=~/mylib/Bismark:$PATH
chmod +x ~/mylib/Bismark/*

###############################################################################
############ ENVIRONMENT READY ############
echo "Setup complete. Bismark is ready to use in the temporary environment."
echo "To exit, run 'deactivate'."

PROJECT_DIR="/projects/rpci/joyceohm/Emily/RRBS-pipeline"
cd $PROJECT_DIR

###############################################################################
# Methylation extraction with NuGen Ovation documentation
#https://felixkrueger.github.io/Bismark/bismark/library_types/
# run with --ingore 3
## this will remove the first 3 bp from the 5' end of read 1 - which will
## remove the mspl cut site from read 1 in our data and avoid artificial 
## methylation calls
# run with --ignore_r2 2 
## this will remove the the first 2 bp from the 5' end of read 1 of pe seq 
## results- the first couple bases of BS-seq have bias toard non-mehylation
## see documentation of mehylation extraction options for more details
# -p is for paired end data 

cd Filtered_Human
files=$(find -name '*mouse.bam')

mkdir ../H2_Methylation_Extraction

# fix sort by name 
printf "%s\n" "$files"| \
  parallel -j 1 \
  samtools sort -n {} -o {}_nsort.bam

files=$(find -name '*_nsort.bam')

printf "%s\n" "$files"| \
  parallel -j 1 \
  bismark_methylation_extractor --gzip --ignore 3 --ignore_r2 2 --p \
  -o ../H2_Methylation_Extraction {}

