# RRBS-pipeline

## Introduction
This directory includes notes and code used to pre-process reduced representation bisulfite sequencing (RRBS) data. Analysis was performed during May 2024 during Emily Isenhart's PhD in the Ohm Lab at RPCC. 

## Background 
Pipeline built with influence from: [NuGEN Technologies](https://github.com/nugentechnologies/NuMetRRBS) and [Peter Fiorka](https://github.com/pfiorica/PDX_RRBS_Processing). I wrote this pipeline entirely in an R markdown format for ease of execution. Shell commands can be run directly to to terminal with cmd + opt + enter. In this instance, as I frequently use shared computational resources, I activate a temporary directory to install neccessary packages for the pipeline. This step can be skipped if the user has write permissions to the default directory. 

* Illumina adapter and diversity adapter trimming with [TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)([github](https://github.com/FelixKrueger/TrimGalore))

## Pipeline overview: 

### I. Illumina Adapter Trimming 


### II. Diversity Adapter Trimming 
