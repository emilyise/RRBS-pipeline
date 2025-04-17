#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --cpus-per-task=1
#SBATCH --mem=300G
#SBATCH --job-name="RRBS_Deduplication"
#SBATCH --output=RRBS_Deduplication.out
#SBATCH --mail-user=emilyise@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute

###############################################################################
############ SET UP SCRIPT ############
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
############ DEDUPLICATION SETUP ############

# Ensure all necessary tools and directories are in place
echo "Setting up deduplication environment..."

PROJECT_DIR="/projects/rpci/joyceohm/Emily/RRBS-pipeline"
SCRATCH="/vscratch/grp-joyceohm/Emily"
DEDUP_OUTPUT_DIR="$SCRATCH/5-Deduplicated"
TEMP_FILES="$DEDUP_OUTPUT_DIR/temp_files"
INDEX_FOLDER="$SCRATCH/RQ025150-Ohm"
SCRIPTS="$PROJECT_DIR/Scripts"

# Ensure deduplication script is executable
chmod +x "$SCRIPTS/nudup.py"

# Create deduplication output directory
mkdir -p "$DEDUP_OUTPUT_DIR"

# Set Directory
cd $SCRATCH/4-Filtered_Human_SAM/temp_files

###############################################################################
############ DEDUPLICATION FUNCTION ############

dedup_sam() {
    local file="$1"
    local index_folder="$2"
    local scripts="$3"
    local output_dir="$4"
    local temp_files="$5"

    # Extract sample name from the SAM file
    sample_name=$(echo "$file" | cut -d "_" -f4)
    echo "Processing file: $file"
    echo "Sample name: $sample_name"

    # Locate the index file
    local index=$(find "$INDEX_FOLDER" -type f -name "*$sample_name*")
    if [ -z "$index" ]; then
        echo "Error: Index file not found for $sample_name"
        return 1
    fi
    echo "Using index file: $index"

    # Create a temporary directory for the sample
    local temp_dir="$temp_files/$sample_name"
    mkdir -p "$temp_dir"

    # Run deduplication
    echo "Running deduplication on $file..."
    python2.7 "$PROJECT_DIR/Scripts/nudup.py" -2 -f "$index" -o "$output_dir/${sample_name}_deduplicated" \
    -T "$temp_dir" --rmdup-only "$file"

    if [ $? -eq 0 ]; then
        echo "Successfully deduplicated $file."
    else
        echo "Deduplication failed for $file!"
        return 1
    fi

    # Cleanup temporary files for this sample after successful processing
    echo "Cleaning up temporary files for $sample_name..."
    rm -rf "$temp_dir"
    if [ $? -eq 0 ]; then
        echo "Temporary files for $sample_name cleaned up successfully."
    else
        echo "Failed to clean up temporary files for $sample_name."
    fi

    echo "Finished processing: $sample_name"
}

# Export variables and function for parallel execution
export -f dedup_sam
export PROJECT_DIR INDEX_FOLDER SCRIPTS DEDUP_OUTPUT_DIR TEMP_FILES

###############################################################################
############ PROCESS SAM FILES FOR DEDUPLICATION ############

# Change to the directory containing sorted and cleaned SAM files
cd "$SCRATCH/4-Filtered_Human_SAM/temp_files"
echo "Current working directory: $PWD"

# Process all unique SAM files in parallel
echo "Starting deduplication of SAM files..."
find . -name '*_cleaned.sam' | parallel -j 5 dedup_sam {} "$INDEX_FOLDER" "$SCRIPTS" "$DEDUP_OUTPUT_DIR" "$TEMP_FILES"

echo "Deduplication process completed. Deduplicated files are in $DEDUP_OUTPUT_DIR."

###############################################################################
############ SUMMARY ############

echo "All SAM files deduplicated."

echo "Final cleanup of temporary files..."
# rm -rf "$TEMP_DIR"
# echo "Temporary files cleaned up."

