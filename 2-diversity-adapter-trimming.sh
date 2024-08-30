###############################################################################
# trim diversity adapters
## pull the script from NuGen's git repo 
## https://github.com/nugentechnologies/NuMetRRBS
cd ..
curl -o Scripts/trimRRBSdiversity.py \
https://raw.githubusercontent.com/nugentechnologies/NuMetRRBS/master/\
trimRRBSdiversityAdaptCustomers.py

# change to py 2.7 to run nugen scripts
conda install python=2.7

## run NuGen diversity trimming script
##for read1 in Trimmed/*R1_001_val_1.fq.gz; do \
##read2=$(echo $read1| sed 's/R1_001_val_1.fq.gz/R3_001_val_2.fq.gz/'); \
##python Scripts/trimRRBSdiversity.py -1 $read1 -2 $read2 -o Trimmed & done

python Scripts/trimRRBSdiversity.py -1 \
'Trimmed/*R1_001_val_1.fq.gz' -2 'Trimmed/*R3_001_val_2.fq.gz' -o Trimmed