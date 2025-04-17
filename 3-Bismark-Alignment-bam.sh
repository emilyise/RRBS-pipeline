#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5G
#SBATCH --job-name="bam-alignment-mouse"
#SBATCH --output=bam-alignment.out
#SBATCH --mail-user=emilyise@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

###############################################################################
############ SET UP ############
# set up temporary environment
python3 -m venv .tempenv
source .tempenv/bin/activate

# set library for packages
mkdir mylib
cd mylib/

# load packages
## load bismark
git clone https://github.com/FelixKrueger/Bismark.git
export PATH=$PWD/Bismark:$PATH
### check
### bismark --version
## load bowtie and samtools
module load gcc/11.2.0
module load bowtie2/2.4.4
module load samtools/1.16.1

# set up conda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install bismark
conda install samtools

# set project directory
PRJ_DIR="/projects/rpci/joyceohm/Emily/RRBS-pipeline"
cd $PRJ_DIR 

# run fastqc
# mkdir Fastqc
# module load fastqc
# find Trimmed -name '*fq_trimmed.fq.gz' | \
# parallel -j 6 fastqc -o Fastqc {}

# for alignment to human genome
## mkdir 2_Aligned_Human
## F_Aligned="/projects/rpci/joyceohm/Emily/RRBS-pipeline/2-Aligned_Human/"
## BISMARK_PATH="/projects/rpci/joyceohm/reference_genome/human_GRCh38/"

# for alignment to mouse genome
mkdir 2-Aligned_Mouse
F_Aligned="/projects/rpci/joyceohm/Emily/RRBS-pipeline/2-Aligned_Mouse/"
BISMARK_PATH="/projects/rpci/joyceohm/reference_genome/mouse_GRCm39/"

# export real paths for later
export BISMARK_PATH=$(realpath $BISMARK_PATH)
export F_Aligned=$(realpath $F_Aligned)

###############################################################################
############ GETTING TRIMMED SAMPLES ############
# pull trimmed files 
tfiles=$(find 1-Trimmed -name '*R1_001_val_1.fq_trimmed.fq.gz')

# Print the found file paths (for debugging)
## echo "Found files:" \
## echo "$tfiles"

# Pull trimmed sample numbers
## Initialize an empty string to store info
s_trimmed=""

## Loop through each found file path
for file_path in $tfiles; do
    ### Use grep and sed to extract the S** portion
    s_string=$(echo "$file_path" | grep -o '_S[0-9]\+' | sed 's/_//')
    ### Append the extracted S** to the s_portions string
    s_trimmed="${s_trimmed}${s_string} "
done

# Remove trailing spaces/newlines
s_trimmed=$(echo "$s_trimmed" | sed 's/[[:space:]]\+$//')

############ GETTING ALIGNED SAMPLES ############
# Pull aligned files
afiles=$(find $F_Aligned -name '*.bam')

# Pull aligned sample numbers
## Initialize an empty string to store info 
s_aligned=""

## Loop through each found file path
for file_path in $afiles; do
    ### Use grep and sed to extract the S** portion
    s_string=$(echo "$file_path" | grep -o '_S[0-9]\+' | sed 's/_//')
    ### Append the extracted S** to the s_portions string
    s_aligned="${s_aligned}${s_string} "
done

# Remove trailing spaces/newlines
s_aligned=$(echo "$s_aligned" | sed 's/[[:space:]]\+$//')

############ GET NOT YET ALIGNED SAMPES ############
s_not_aligned=""
for s in $s_trimmed; do
    if ! echo "$s_aligned" | grep -q "$s"; then
        s_not_aligned="${s_not_aligned}${s} "
    fi
done

# Remove trailing spaces/newlines
s_not_aligned=$(echo "$s_not_aligned" | sed 's/[[:space:]]\+$//')

############ SUBSET TO NEEDED ALIGNMENTS ############
# find the files 
matching_files=()
index=0

for s in $s_not_aligned; do
    # Find files in '1-Trimmed' that match each S**
    files=$(find 1-Trimmed -name "*_${s}_*")
    files=$(find $files -name '*R1_001_val_1.fq_trimmed.fq.gz')

    # Loop through each found file and assign to the array using an index
    while read -r file; do
        matching_files[index]="$file"  # Assign to current index
        index=$((index + 1))            # Increment the index
    done <<< "$files"
done

############ RUN ALIGNMENTS ############
# run alignment 
## this is where we need the real paths we made earlier
printf "%s\n" "${matching_files[@]}" | \
  parallel -j 12 'read2=$(echo {} | sed "s/R1_001_val_1.fq_trimmed.fq.gz/R3_001_val_2.fq_trimmed.fq.gz/"); \
  bismark --score_min L,0,-0.4 --genome_folder $BISMARK_PATH --output_dir $F_Aligned -1 {} -2 $read2'

###############################################################################