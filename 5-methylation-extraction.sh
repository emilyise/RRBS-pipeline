#!/bin/sh
#SBATCH --time=24:00:00
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --job-name="methylation-_extraction"
#SBATCH --output=methylation_extraction.out
#SBATCH --mail-user=emilyise@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

# set up the environment
eval "$(/util/common/python/py38/anaconda-2021.05/bin/conda shell.bash hook)"
conda create --prefix /tmp/eminenv2 python=2.7
conda activate /tmp/eminenv2

conda install bismark samtools

# set directory
export PROJECT_PATH="/projects/rpci/joyceohm/Emily/RRBS-pipeline"
cd $PROJECT_PATH

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

cd ./B1-Aligned-Done

mkdir ./B1-Methylation_Extraction

find -name '*_pe.sam' | \
  parallel -j 6 \
  bismark_methylation_extractor --gzip --bedGraph --ignore 3 --ignore_r2 2 --p \
  -o ../B1-Methylation_Extraction {}
  
  




