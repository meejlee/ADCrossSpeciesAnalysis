%% Import_MouseData.m
% Overview: Import mouse phenotype and gene expression data from GEO/R to
% MATLAB

%% Mouse
cd '../../data/external/GSE64398/data'
GSE64398 = importdata('GSE64398_raw_norm.txt');
GSE64398_rma = GSE64398.data;
GSE64398_genename = GSE64398.textdata(2:end,5);

%%%% Phenotype
cd '../../..'
GSE64398_phno = readtable('GSMList_GSE64398mouse.csv');

cd '../scripts/MATLAB'