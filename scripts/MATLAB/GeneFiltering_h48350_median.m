%% Run_GeneFiltering_h48350_Median.m
% Overview: Match human and mouse probes, Collapse probe level measurements into a single value per gene (median)

%% Filtering human DEG list to ID unique list of genes
Hu48350_idx = GSE48350_Dis.textdata(2:end,1); %gene index in full matrix
Hu48350_name = GSE48350_Dis.textdata(2:end,2); %gene name

% sort to unique
[sorted_HuGenes_unique,idx,~] = unique(Hu48350_name,'stable');  

%% Trim if # of genes is greater than 2 orders of magnitude against # of samples
if length(sorted_HuGenes_unique) > size(GSE48350_rma,2)*100
    sorted_HuGenes_unique = sorted_HuGenes_unique(1:size(GSE48350_rma,2)*100);
end
    
%% Matching mouse and human gene names
% Init variable to hold indices for each gene
idxMm_Hu48350Unique = cell(length(sorted_HuGenes_unique),1);

% Variable type reformatting
mouse_homologs_updated = cell(size(mouse_homologs));
for i = 1:length(mouse_homologs) 
    mouse_homologs_updated{i} = char(mouse_homologs{i});
end

% Find index for each human gene from initial DEG list
for iter_hu = 1:length(sorted_HuGenes_unique)
    loop_idx = find(ismember(mouse_homologs_updated,sorted_HuGenes_unique{iter_hu}));
    idxMm_Hu48350Unique{iter_hu} = loop_idx;
end

%% Median collapse probes that match for each species
%% Human Samples
Hu48350_rma_median = zeros(length(sorted_HuGenes_unique),size(GSE48350_rma,2));
for iter_gene = 1:length(sorted_HuGenes_unique)
   loop_idx = find(ismember(GSE48350_genename,sorted_HuGenes_unique{iter_gene}));
   for px = 1:size(GSE48350_rma,2)
       temp_val = median(GSE48350_rma(loop_idx,px),'omitnan');
       Hu48350_rma_median(iter_gene,px) = temp_val;
   end
end

%% Mouse Samples
%%%% Init variable
Mm64398_rma_median = zeros(length(sorted_HuGenes_unique),size(GSE64398_rma,2));
%%%% Median collapse
for i = 1:length(sorted_HuGenes_unique)
   loop_idx = idxMm_Hu48350Unique{i}; 
   for px = 1:size(GSE64398_rma,2)
       Mm64398_rma_median(i,px) = median(GSE64398_rma(loop_idx,px),'omitnan');
   end
end
%%%% Remove any samples with NaN values in the median (from unmapped
%%%% genes)
findNaN = sum(isnan(Mm64398_rma_median),2);
idxNaN = find(findNaN == size(Mm64398_rma_median,2));
Mm64398_rma_median(idxNaN,:) = [];

%% Remove genes that didn't match from human list
% Remove from gene list
sorted_HuGenes_unique_noNan = sorted_HuGenes_unique;
sorted_HuGenes_unique_noNan(idxNaN,:) = [];
% Remove from dataset
Hu48350_rma_median(idxNaN,:) = [];

%% Initialize variable
X_hu = Hu48350_rma_median;