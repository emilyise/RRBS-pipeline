#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --job-name="remove-mouse-reads"
#SBATCH --output=remove-mouse-reads.out
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

# load required modules
echo "Loading required modules."
module load gcc/11.2.0
module load samtools/1.16.1
echo "Modules loaded."

# set project directories
PRJ_DIR="/projects/rpci/joyceohm/Emily/RRBS-pipeline/"
cd $PRJ_DIR 

WORKING_DIR="/vscratch/grp-joyceohm/Emily"
HUMAN_DIR="$WORKING_DIR/2-Aligned_Human"
MOUSE_DIR="$WORKING_DIR/2-Aligned_Mouse"
OUTPUT_DIR="$WORKING_DIR/3-Filtered_Human"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

###############################################################################
############ REMOVE MOUSE READS FOR EACH SAMPLE ############
echo "Starting filtering process for multiple samples..."

# Loop through human BAM files
for HUMAN_BAM in "$HUMAN_DIR"/*_sorted.pe.bam; do
    # Extract the base name to find corresponding mouse BAM
    BASENAME=$(basename "$HUMAN_BAM" | sed 's/_sorted.pe.bam//')

    # Define output BAM and its index file
    FILTERED_HUMAN_BAM="$OUTPUT_DIR/${BASENAME}_filtered_no_mouse.bam"
    FILTERED_HUMAN_BAI="${FILTERED_HUMAN_BAM}.bai"

    # Skip processing if the filtered BAM index file already exists
    if [ -f "$FILTERED_HUMAN_BAI" ]; then
        echo "Filtered BAM for sample $BASENAME already completed. Skipping..."
        continue
    fi

    # Define corresponding mouse BAM file
    MOUSE_BAM="$MOUSE_DIR/${BASENAME}_sorted.pe.bam"

    # Check if mouse BAM file exists
    if [ -f "$MOUSE_BAM" ]; then
        echo "Processing sample: $BASENAME"

        # Preprocess human and mouse BAM files to keep only aligned reads
        HUMAN_ALIGNED_BAM="${HUMAN_BAM%.bam}_aligned.bam"
        MOUSE_ALIGNED_BAM="${MOUSE_BAM%.bam}_aligned.bam"

        echo "Filtering aligned reads for human BAM: $HUMAN_BAM"
        samtools view -b -F 4 "$HUMAN_BAM" > "$HUMAN_ALIGNED_BAM"

        echo "Filtering aligned reads for mouse BAM: $MOUSE_BAM"
        samtools view -b -F 4 "$MOUSE_BAM" > "$MOUSE_ALIGNED_BAM"

        # Index the aligned BAM files
        samtools index "$HUMAN_ALIGNED_BAM"
        samtools index "$MOUSE_ALIGNED_BAM"

        # Remove mouse-aligned reads from human BAM
        echo "Filtering mouse reads from $HUMAN_ALIGNED_BAM..."
        samtools view -h "$HUMAN_ALIGNED_BAM" | \
            grep -v -f <(samtools view "$MOUSE_ALIGNED_BAM" | cut -f1) | \
            samtools view -Sb - > "$FILTERED_HUMAN_BAM"

        # Index the filtered BAM file
        echo "Indexing filtered human BAM: $FILTERED_HUMAN_BAM"
        samtools index "$FILTERED_HUMAN_BAM"
    else
        echo "Mouse BAM file not found for sample: $BASENAME. Skipping..."
    fi
done

echo "All samples processed. Filtered files are in $OUTPUT_DIR."

###############################################################################
############ CLEAN UP ############
# Deactivate and clean up the temporary environment
echo "Cleaning up temporary environment..."
deactivate
rm -rf .tempenv
echo "Temporary environment removed."



















