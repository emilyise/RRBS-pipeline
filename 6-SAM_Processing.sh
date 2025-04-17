#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=120G
#SBATCH --job-name="RRBS_SAM_Processing"
#SBATCH --output=RRBS_SAM_Processing.out
#SBATCH --mail-user=emilyise@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute

###############################################################################
set -ex # 'e' for exit on error, 'x' for command tracing

###############################################################################
############ SET UP TEMPORARY ENVIRONMENT ############
echo "Establishing temporary environment..."

# set up temporary environment
python3 -m venv .tempenv2
source .tempenv2/bin/activate

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

############ PYTHON 2 SETUP ############
# Check for Python 2
if command -v python2.7 &>/dev/null; then
    echo "Python 2.7 is already installed. Proceeding to the next step."
else
    echo "Python 2.7 not found. Installing it locally..."

    # Navigate to the Python 2 directory
    mkdir -p "$PYTHON2_DIR"
    cd "$PYTHON2_DIR"

    # Download Python 2.7 source code
    echo "Downloading Python 2.7..."
    curl -O "$PYTHON2_URL"

    # Extract and install Python 2.7 locally
    echo "Extracting Python 2.7..."
    tar -xzf Python-2.7.18.tgz
    cd Python-2.7.18

    echo "Configuring Python 2.7 installation..."
    ./configure --prefix="$(pwd)/localpython2"

    echo "Building and installing Python 2.7 (this may take a while)..."
    make && make install

    # Add local Python 2.7 to PATH within the temporary environment
    export PATH="$(pwd)/localpython2/bin:$PATH"
    echo "Python 2.7 installed successfully."
    cd $PROJECT_DIR
    cd "$MYLIB_DIR"
fi

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

###############################################################################
############ ENVIRONMENT READY ############
echo "Setup complete. Python 2.7 and Bismark are ready to use in the temporary environment."
echo "To exit, run 'deactivate'."

cd $PROJECT_DIR

# pull NuDup script and sam stripper script
curl -o Scripts/nudup.py \
https://raw.githubusercontent.com/tecangenomics/nudup/master/nudup.py

curl -o Scripts/strip_bs.sh \
https://raw.githubusercontent.com/nugentechnologies/NuMetRRBS/master/\
strip_bismark_sam.sh

###############################################################################
############ CONVERT TO SAM ############
# Define the input and output directories
INPUT_DIR="/vscratch/grp-joyceohm/Emily/3-Filtered_Human"  # Directory containing the BAM files
OUTPUT_DIR="/vscratch/grp-joyceohm/Emily/4-Filtered_Human_SAM"  # Directory to store the SAM files
#
# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"
# 
echo "Converting all BAM files in $INPUT_DIR to SAM files in $OUTPUT_DIR..."
# 
# Loop through each BAM file in the input directory
# for BAM_FILE in "$INPUT_DIR"/*.bam; do
#     # Check if there are any BAM files in the directory
#     if [ ! -e "$BAM_FILE" ]; then
#        echo "No BAM files found in $INPUT_DIR."
#       exit 1
#     fi
# 
#     # Extract the base name of the BAM file
#     BASENAME=$(basename "$BAM_FILE" .bam)
#     
#     # Define the output SAM file path
#     SAM_FILE="$OUTPUT_DIR/${BASENAME}.sam"
#  
#     # Convert BAM to SAM using samtools
#     echo "Converting $BAM_FILE to $SAM_FILE..."
#     samtools view -h -o "$SAM_FILE" "$BAM_FILE" || { 
#     echo "Error converting $BAM_FILE"; 
#     exit 1; 
#     }
# 
# 
#     # Check if the conversion was successful
#     if [ $? -eq 0 ]; then
#         echo "Successfully converted $BAM_FILE to $SAM_FILE."
#     else
#         echo "Failed to convert $BAM_FILE. Skipping..."
#     fi
# done
# 
# echo "Conversion process completed. All SAM files are in $OUTPUT_DIR."

###############################################################################
############ SETUP AND PREPARATION ############

PROJECT_DIR="/projects/rpci/joyceohm/Emily/RRBS-pipeline"
cd $PROJECT_DIR
# 
# Make the SAM stripper script executable
echo "Making SAM stripper script executable..."
chmod +x Scripts/strip_bs.sh
# 
# Define working and output directories
WORKING_DIR="/vscratch/grp-joyceohm/Emily/4-Filtered_Human_SAM"
TEMP_DIR="./temp_files"
 
echo "Changing to working directory: $WORKING_DIR"
cd "$WORKING_DIR"

# Create necessary directories
mkdir -p "$TEMP_DIR"
    
###############################################################################
############ SORT SAM FILES ############
find . -name '*.sam' | parallel -j 16 '\
    SAM_FILE={};\
    BASENAME=$(basename "$SAM_FILE" .sam);\
    SORTED_SAM="./temp_files/${BASENAME}_sorted.sam";\
    
    echo "Sorting $SAM_FILE to $SORTED_SAM...";\
    samtools sort -O SAM -o "$SORTED_SAM" "$SAM_FILE"\
    || { echo "Error sorting $SAM_FILE"; exit 1; };\
    
    echo "Sorting completed for $SAM_FILE."\
'


  
 
 
