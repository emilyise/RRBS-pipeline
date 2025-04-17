#!/bin/sh
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10
#SBATCH --job-name="illumina-adapter-trimming"
#SBATCH --output=illumina-adapter-trimming.out
#SBATCH --mail-user=emilyise@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

# set up the environment
eval "$(/util/common/python/py38/anaconda-2021.05/bin/conda shell.bash hook)"
conda create --prefix /tmp/test-env python=3
conda activate /tmp/test-env

## load Bioconda
## documentation: https://bioconda.github.io/
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

## install trim galore to temporary environment 
conda install trim-galore

#set directory
PROJECT_DIR=$1
cd $PROJECT_DIR

mkdir 1-Trimmed

find -name "*R1_001.fastq.gz" | cut -d "_" -f1 | 
parallel -j 1 trim_galore \
--paired -o ./1-Trimmed {}*R1_001.fastq.gz {}*R3_001.fastq.gz

