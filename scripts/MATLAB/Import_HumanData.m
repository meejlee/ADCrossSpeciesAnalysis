%% Run_ImportHumanData.m
% Overview: Import normalized gene expression data from R and phenotype information from GEO into MATLAB

%% GSE1297 (human)
%%%% Gene Expression:
cd '../../data/external/GSE1297/data'
GSE1297 = importdata('GSE1297_rmaprobes.txt'); %probe level expression
GSE1297_genename = GSE1297.textdata(2:end,4); %probe gene names
GSE1297_ENSEMBL = GSE1297.textdata(2:end,7); %probe ENSEMBL IDs
GSE1297_rma = GSE1297.data; %pull expression value only into numeric matrix

%% GSE48350 (human)
%%%% Gene Expression:
cd '../../GSE48350/data'
GSE48350 = importdata('GSE48350_rmaprobes.txt'); %probe level expression
GSE48350_genename = GSE48350.textdata(2:end,4);
GSE48350_ENSEMBL = GSE48350.textdata(2:end,7);
GSE48350_rma = GSE48350.data;

GSE48350_Dis = importdata('GSE48350degs_coef1_pvalonly_Dis.txt');
    % file: list of DEGs (adjusted p-value < 0.20 between control and AD
    % patients (binary label)

%% Phenotype:
cd '../../..'
GSE1297_phno = readtable('GSMList_GSE1297human.csv');
GSE48350_phno = readtable('GSMList_GSE48350human.csv');

%%
cd '../scripts/MATLAB'