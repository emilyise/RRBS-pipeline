# RRBS-pipeline

## Introduction
This directory includes notes and code used to pre-process reduced representation bisulfite sequencing (RRBS) data. Analysis was performed during May 2024 during Emily Isenhart's PhD in the Ohm Lab at RPCC. 

## Background 
This is a ipeline built for NuGEN Ovation RRBS Methyl-Seq with heavy influence from [NuGEN Technologies](https://github.com/nugentechnologies/NuMetRRBS). I wrote this pipeline in an R markdown format for ease of execution. Shell commands can be run directly to to terminal with cmd + optn + enter. In this instance, as I frequently use shared computational resources, I activate a temporary directory to install neccessary packages for the pipeline. This step can be skipped if the user has write permissions to the default directory. 

* module load with [Bioconda](https://bioconda.github.io/)
* Adapter trimming with [TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)([github](https://github.com/FelixKrueger/TrimGalore))
* Diversity adapter trimming with [NuGEN](https://github.com/nugentechnologies/NuMetRRBS) custom script
* QC with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)([github](https://github.com/s-andrews/FastQC))
* Alignment to genome with [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)
* Deduplication with [NuGen](https://github.com/nugentechnologies/NuMetRRBS) custom script

## Pipeline overview: 

### I. Adapter Trimming
First, adapters need to be trimmed from the 3' end of read. We use Trim Galore with Illumina adapter sequences: AGATCGGAAGAGC (read 1) & AAATCAAAAAAAC (read 2). If no sequence is supplied, TrimGalore will try to auto-detec the adapter, or other assay specific sequences can be supplied (ex. Nextera: CTGTCTCTTATA (single)). For NuGEN Ovation, this should be run without the --rrbs option. NuGEN utilizes diversity adapters (describled in II.) which will be errantly trimmed if run with --rrbs. 

### II. Diversity Adapter Trimming 
Sequence complexity is an important consideration on Illimina sequences. If sequences lack complexity, cluster identification and color matrix calculations can be negatively impacted. RRBS libraries all begin with CGG or TCC. To ensure that all four bases are present in the first few cycles of sequencing, NuGEN Ovation includes the addition of 0-3 bp diversity adapters in between the sequencing primer and the fragment. Thus, the diversity adapters must be trimmed after illumina adapter sequences are trimmed. NuGEN has published a script for this purpose, which we apply in this pipeline. 

### III. FastQC
After adapter trimming, we perfrom FastQC to perform quality control checks on raw sequence data. We expect per base sequence content to be flaged in a FastQC report for RRBS. Position 2 and 3 (as previously discussed) should contain 100% G sequences and the overall content should be sketwed with higher T content and lower C content than other non-targeted assays. Per sequence GC content should be higher than the theoretical distriubition, as RRBS is biased towed CpG islands. Likely to also identify warnings for overrepresented sequencences and sequence duplication level. 

### IV. Alignment
Aligned with bismark. Note that genome preparation for bismark is NOT performed with this pipeline- files were previously prepared using GrCh38. Run with --sam for use with NuGEN specific deduplicaiton. (if not run with --sam, will have to strip bam files for use with nudup). 

### V. Deduplication
Performed with NuGEN custom script. Can't run genomic repeat based deduplication on RRBS (targeted). But, NuGEN introduces unique molecular tags for this purpose, therfore deduplication can be performed with the custom script. 


