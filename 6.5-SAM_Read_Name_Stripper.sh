#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=10
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --job-name="SAM_Stripping"
#SBATCH --output=SAM_Stripping.out
#SBATCH --mail-user=emilyise@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute

###############################################################################
############ SET UP SCRIPT ############
# Exit on error
set -euo pipefail

###############################################################################
############ SET UP TEMPORARY ENVIRONMENT ############
echo "Loading required modules..."

module load gcc/11.2.0
module load samtools

PROJECT_DIR="/vscratch/grp-joyceohm/Emily"
DATA_DIR="${PROJECT_DIR}/4-Filtered_Human_SAM/temp_files"
###############################################################################
############ STRIP SAM READ NAMES WITH PARALLEL ############
echo "Stripping read names in parallel..."

cd $DATA_DIR

# Function to clean SAM file
clean_sam_file() {
    samfile=$1
    echo "Processing $samfile..."

    # Extract shorter sample name — up to "S##_L###"
    shortname=$(echo "$samfile" | grep -oP 'RS-[^_]+_[^_]+_RS-[^_]+_S[0-9]+_L[0-9]+')
    cleaned_sam="${shortname}_cleaned.sam"

    # Clean SAM read names
    samtools view -h "$samfile" \
    | sed -E 's/^(@?[A-Z0-9]+:[0-9]+:[A-Z0-9]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+)_[12]:N:0:[A-Z]+/\1/' \
    > "$cleaned_sam"

    echo "Cleaned SAM written to $cleaned_sam"
}

export -f clean_sam_file  # Export the function to use it with parallel

# Run the cleaning process in parallel
ls *_sorted.sam | parallel clean_sam_file
done

# Clean out sorted.sam files from directory
rm *_sorted.sam