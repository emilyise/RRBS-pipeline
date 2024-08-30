#!/bin/sh
#SBATCH --time=72:00:00
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=20G
#SBATCH --job-name="sam-alignment"
#SBATCH --output=sam-alignment.out
#SBATCH --mail-user=emilyise@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

# set up the environment
eval "$(/util/common/python/py38/anaconda-2021.05/bin/conda shell.bash hook)"
# set project directory

export PRJ_DIR="/projects/rpci/joyceohm/Emily/RRBS-pipeline"
cd $PRJ_DIR 

# run fastqc
# mkdir Fastqc
# module load fastqc
# find Trimmed -name '*fq_trimmed.fq.gz' | \
#  parallel -j 6 fastqc -o Fastqc {}

# run fastqscreen
## downloading and configuring fastqscreen: 
## wget https://github.com/StevenWingett/FastQ-Screen/archive/refs/tags/v0.15.3.tar.gz
## then v0.15.3.tar.gz
## cd FastQ-Screen-0.15.3/
## cp fastq_screen.conf.example fastq_screen.conf
## I'll come back to this later because I simply could not get it it find my 
## executable files for bowtie2 and bismark in my temp env

conda activate bismark

mkdir Aligned

export BISMARK_PATH="/projects/rpci/joyceohm/Shared_Resources/ICRG/"

find Trimmed -name '*R1_001_val_1.fq_trimmed.fq.gz' | \
  parallel -j 3 'read2=$(echo {} | sed "s/R1_001_val_1.fq_trimmed.fq.gz/R3_001_val_2.fq_trimmed.fq.gz/"); \
  bismark --sam --score_min L,0,-0.4 -o Aligned $BISMARK_PATH -1 {} -2 $read2'
