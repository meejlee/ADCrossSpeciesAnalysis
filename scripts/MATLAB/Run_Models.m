%% Run_Models.m
% Overview: Conduct all modeling
% 1. Generate human PCA space
% 2. Run linear model age vs disease comparison
% 3. Conduct LASSO + PCR

%% Run initial PCA for human data
% Reshape matrix if necessary
if (size(X_hu,1) > size(X_hu,2)) == 1
    X_hu = transpose(X_hu);
end
[coeffHu, scoresHu, ~, ~, varExplainedHu, pca_centersHu] = pca(X_hu,'Algorithm','svd');

%% Calculate # PCs to keep
%%%% ID # of components to keep
% EITHER: 95% cumulative variance explained
% OR: each retained PC explains at least 1% of variance
% Use more restrictive criteria

explained_variance_to_keep = 95/100;
nComponents95 = find(cumsum(varExplainedHu)/sum(varExplainedHu) >= explained_variance_to_keep, 1);
for iter_nC = 1:length(varExplainedHu)
   if varExplainedHu(iter_nC) < 1
       break
   end
end
if iter_nC < nComponents95
    nComponents95 = iter_nC - 1;
end

%% Format mouse phenotype data
idxMmDis = vertcat(idxMouseGroups{2},idxMouseGroups{3}); 
filenameMmStrain = 'm64398TASTPM';
    
% merge control and disease indexing
idxMmFinal = vertcat(idxMouseGroups{1},idxMmDis);

%%%% Format Y (mouse outcomes)
YMmage = phnoMm(idxMmFinal,2)/72; %age in months (normalized to max age; age originally in weeks)
YMmdisease = phnoMm(idxMmFinal,1);

%% Project mouse data into human PCA
%%%% Z-score within the samples used and then project
data_mm = zscore(transpose(log2(Mm64398_rma_median(:,idxMmFinal))));
X_mm = data_mm * coeffHu(:,1:nComponents95);

%% Run Linear modeling (age vs disease)
run 'Model_Linear_AgeVSDisease.m';

%% Run Principal Component Regression (with LASSO)
run 'Model_RepeatLASSOLM.m';
% To run once:
% run 'Model_LASSOandLM.m';

run 'Model_LMfixedInputs.m';

%% Testing true model against null models
run 'Model_NullModelPCs.m'
run 'Model_NullPhenotype.m'