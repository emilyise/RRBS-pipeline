#!/bin/sh
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --job-name="RRBS_Deduplication"
#SBATCH --output=RRBS_Deduplication.out
#SBATCH --mail-user=emilyise@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

# set up the environment
eval "$(/util/common/python/py38/anaconda-2021.05/bin/conda shell.bash hook)"
conda create --prefix /tmp/eminenv2 python=2.7
conda activate /tmp/eminenv2

## load Bioconda
## documentation: https://bioconda.github.io/
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# set directory
export PROJECT_PATH="/projects/rpci/joyceohm/Emily/RRBS-pipeline"
cd $PROJECT_PATH

# pull NuDup script and sam stripper script
curl -o Scripts/nudup.py \
https://raw.githubusercontent.com/tecangenomics/nudup/master/nudup.py

curl -o Scripts/strip_bs.sh \
https://raw.githubusercontent.com/nugentechnologies/NuMetRRBS/master/\
strip_bismark_sam.sh

# run sam stripper
#chmod +x Scripts/strip_bs.sh
#find B1-Aligned-Done/ -name '*pe.sam' | \
#  parallel -j 9  \
#  'Scripts/strip_bs.sh' {}

# install necessary modules
conda install "samtools>=1.10" GNU-coreutils bowtie2

cd ./B1-Aligned-Done/

# sort sam files for deduplication
find -name '*sam_stripped.sam' | \
    parallel -j 3 '\
    myfile={};\
    base_name=$(basename {} .sam_stripped.sam);\
    samtools sort -O SAM -o "${base_name}_sorted.sam" "$myfile"'
    
# pull only unique alignments for deduplication
find -name '*_sorted.sam' | \
    grep "S31" | \
    parallel -j 1 '\
    myfile={};\
    base_name=$(basename {} .sam);\
    samtools view -h "$myfile" | grep "NH:i:1"'
    > "${base_name}_unique.sam";\
    echo done $base_name'

#run deduplication
chmod +x Scripts/nudup.py
find -name '*sorted.sam' | \
  grep "S31" | \
  parallel -j 1 '\
  myfile={};\
  base_name=$(basename {} .sam);\
  index_name=$(echo $base_name | sed "s/R[1-3]_001_val_[1-2].fq_trimmed_bismark_bt2_pe_sorted//");\
  index_name="${index_name}R2_001.fastq.gz";\
  index_name=$(find -name $index_name);\
  python Scripts/nudup.py -2 -f "$index_name" -o "$base_name" "$myfile"'
  
 
 
