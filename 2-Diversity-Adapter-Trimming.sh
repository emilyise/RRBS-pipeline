###############################################################################
# trim diversity adapters
## pull the script from NuGen's git repo 
## https://github.com/nugentechnologies/NuMetRRBS
curl -o Scripts/trimRRBSdiversity.py \
https://raw.githubusercontent.com/nugentechnologies/NuMetRRBS/master/\
trimRRBSdiversityAdaptCustomers.py

# activate temporary environment
conda activate /tmp/test-env

# change to py 2.7 to run nugen scripts
conda install python=2.7

## run NuGen diversity trimming script

python Scripts/trimRRBSdiversity.py -1 \
'Trimmed/*R1_001_val_1.fq.gz' -2 'Trimmed/*R3_001_val_2.fq.gz' -o 1-Trimmed