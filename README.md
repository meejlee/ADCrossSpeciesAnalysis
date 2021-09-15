# ADCrossSpeciesAnalysis
Cross-species analysis of AD data via Translatable Component Regression and mixed linear modeling

# Overview
Publicly available datasets analyzed: GSE48350, GSE1297 (Human postmortem hippocampal transcriptomics), GSE64398 (Mouse hippocampal transcriptomics)

# Getting started
## Import Data
Start with .R files in ../scripts to import samples and normalize data
Scripts for importing are named "ImportGSE48350.R", "ImportGSE1297.R", and "ImportGSE64398.R" and are not dependent on one another.
Imported data saves to ../data/external appended by the study name

NOTE: Manual downloads are required for GSE64398. Additional details are commented into the ImportGSE64398.R script.

## Run TransComp-R
Move to ../scripts/MATLAB to implement TransComp-R
### Importing data into MATLAB
Start with Run_InitMATLAB.m to load data from R and conduct pre-processing (human-to-mouse homolog matching, median collapsing). This script calls all .m files that start with Import_ while GeneFiltering_ and HomologMatching .m files are called as necessary

NOTE: Run_InitMATLAB.m has flags at the top to specify the human and mouse study of interest. Case study 1 uses the human flag 'h48350', and case study 2 uses the human flag 'h1297'.

.mat files are saved in ../data/interim as generated

### Running Models
After loading data into MATLAB, run Run_Models.m

This script calls all .m files that start with Model_ and executes TransComp-R modeling, age vs disease linear modeling, and null model testing

Modeling results are saved as .mat files in ../modeling outputs for subsequent plotting

### Run Plotting
This script calls all .m files that start with Plot_ and regenerates figures in the manuscript

NOTE: An additional .R scripts is provided in ../scripts to regenerate Figures 2B,D (HumanPhenotypeAnova.R)

# Additional Information
Human phenotype information is directly from each respective public study via GEO and saved in ../data as .csv files to match the formatting anticipated in the scripts provided. Human-to-mouse homolog conversion was conducted using the same method as Kumar, et al. Cell, 2017 (https://doi.org/10.1016/j.celrep.2018.10.047) using the same .txt file provided in the original paper's supplement (included).
