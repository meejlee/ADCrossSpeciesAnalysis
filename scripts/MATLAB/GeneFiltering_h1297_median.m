%% Run_GeneFiltering_h1297_median.m
% Overview: Match human and mouse probes, Collapse probe level measurements into a single value per gene (median)

%%
% Identify unique gene names
sorted_HuGenes_unique = unique(GSE1297_genename,'stable');

% Variable type reformatting
mouse_homologs_updated = cell(size(mouse_homologs));
for iter_gene = 1:length(mouse_homologs) 
    mouse_homologs_updated{iter_gene} = char(mouse_homologs{iter_gene});
end

% Init variable to hold indices for each gene
idxMm_Hu1297Unique = cell(length(sorted_HuGenes_unique),1);

% Find index for each human gene from initial DEG list
for iter_hu = 1:length(sorted_HuGenes_unique)
    loop_idx = find(ismember(mouse_homologs_updated,sorted_HuGenes_unique{iter_hu}));
    idxMm_Hu1297Unique{iter_hu} = loop_idx; %store multiple values if multiple hits
end

%% Median collapse probes that match
%% Human Samples 
Hu1297_rma_median = zeros(length(sorted_HuGenes_unique),size(GSE1297_rma,2)); %# Genes x # samples
for iter_gene = 1:length(sorted_HuGenes_unique)
   loop_idx = find(ismember(GSE1297_genename,sorted_HuGenes_unique{iter_gene}));
   for px = 1:size(GSE1297_rma,2)
       temp_val = median(GSE1297_rma(loop_idx,px),'omitnan');
       Hu1297_rma_median(iter_gene,px) = temp_val;
   end
end

%% Mouse Samples
%%%% Init variable
Mm64398_rma_median = zeros(length(sorted_HuGenes_unique),size(GSE64398_rma,2)); % # DEGs x # Mouse Samples
%%%% Median collapse
for iter_gene = 1:length(sorted_HuGenes_unique)
   loop_idx = idxMm_Hu1297Unique{iter_gene}; 
   for px = 1:size(GSE64398_rma,2)
       Mm64398_rma_median(iter_gene,px) = median(GSE64398_rma(loop_idx,px),'omitnan');
   end
end
%%%% Remove any samples with NaN values in the median (from unmapped
%%%% genes)
findNaN = sum(isnan(Mm64398_rma_median),2);
idxNaN = find(findNaN == size(Mm64398_rma_median,2));
Mm64398_rma_median(idxNaN,:) = [];

%% Remove genes that didn't match from the human list
% Remove from gene list
sorted_HuGenes_unique_noNan = sorted_HuGenes_unique;
sorted_HuGenes_unique_noNan(idxNaN,:) = [];

% Remove from dataset
Hu1297_rma_median(idxNaN,:) = [];

%% Initialize variable
X_hu = Hu1297_rma_median;
