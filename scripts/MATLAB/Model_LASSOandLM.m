%% Model_LASSOandLM.m
% Build linear model using mouse samples projected on human PCs (X) and mouse
% disease background (Y)
% LASSO (alpha = 1.0) for feature selection prior to linear model buildling

%% LASSO
if (exist('idxTrim','var') == 1); clear idxTrim; end % remove if exists from previous runs; can cause script to throw errors if variable is already initialized
Y_mm = YMmdisease;
noSamples = length(Y_mm);
nfold = 5;
nThresh = 0.4;
alphaFinal = 1;

% Init variables
yPredVTrueLM = nan(2,noSamples);
yPredVTrueLM(1,:) = Y_mm;
PCFqCumulative = nan(nComponents95,1);
PCFqKfold = nan(nComponents95,nfold);
PCLassoVals = nan(nComponents95,2);

CVO = cvpartition(Y_mm,'kfold',nfold);
WeightsLASSO_store = cell(nfold,1);

for iterCV = 1:nfold
    % Split test + train
    idx_te = find(CVO.test(iterCV)==1);
    idx_tr = find(CVO.training(iterCV)==1);

    Y_tr = Y_mm(idx_tr);
    X_LASSOtr = X_mm(idx_tr,:);

    nStab=10; nCut=nThresh*nStab; Sel_matrix=zeros(nStab,nComponents95);
    WeightsLASSO_loop = zeros(nComponents95,nStab);
    for k=1:nStab
        [B, FitInfo] = lasso(X_LASSOtr, Y_tr,'CV',5,'Alpha',alphaFinal); %Fold-specific variable selection using LASSO using ONLY fold-specific training data. Lambda tuned using an internal CV on the fold-specific training data ONLY.
        chosen_b=max(FitInfo.Index1SE,1);
        WeightsLASSO_loop(:,k) = B(:,chosen_b);
        B_fin=B(:,chosen_b); var_selec_pre=find(abs(B_fin)>0); Sel_matrix(k,var_selec_pre)=1;
    end

    % Store LASSO weights for each of the 10 repeats
    WeightsLASSO_store{iterCV} = WeightsLASSO_loop;
end

%%%% Storing info from multiple rounds of LASSO CV
% Pull weights from all LASSO runs (nStab times per nfold CV folds)
% Reformatting to pull from cell array into a single matrix
WeightsLASSO = nan(nComponents95,nfold*nStab);
counter = 1;
for iterCV = 1:nfold
    matLoop = WeightsLASSO_store{iterCV};
    idxloop = (counter):1:(counter+nStab-1);
    WeightsLASSO(:,idxloop) = matLoop;
    counter = counter + nStab;
end

% Frequency of PC selection across nStab * nfold # of LASSO runs
for i = 1:nComponents95
    PCFqCumulative(i,1) = nnz(WeightsLASSO(i,:));
end

% Frequency of PC selection by nfold
for i = 1:nComponents95
    for j = 1:nfold
        PCFqKfold(i,j) = nnz(WeightsLASSO(i,((j-1)*nStab+1):(j)*nStab));
    end
end

% Mean and StDev values for LASSO weights across all runs
for i = 1:nComponents95
    PCLassoVals(i,1) = mean(WeightsLASSO(i,:));
    PCLassoVals(i,2) = std(WeightsLASSO(i,:))/sqrt(length(WeightsLASSO));
end

%% Linear Model
%% Final Model with all samples
% Pull PCs that were selected 40% + of the time for the final model
counter = 1;
for m = 1:nComponents95
   if PCFqCumulative(m) >= (nfold*nStab*nThresh)
       idxTrim(counter) = m;
       counter = counter + 1;
   end
end

% Build linear model
try
    Model_LM = fitlm(X_mm(:,idxTrim),Y_mm);
    for LOO = 1:length(idxMmFinal)
        % Split test + train
        idx_te = LOO;
        idx_tr = 1:1:length(idxMmFinal); 
        idx_tr(LOO) = [];
        Y_tr = Y_mm(idx_tr);
        X_LMtr = X_mm(idx_tr,idxTrim);
        % Build model and predict
        Model_LMtr = fitlm(X_LMtr,Y_tr);
        X_LMte = X_mm(idx_te,idxTrim);
        Y_predict = predict(Model_LMtr,X_LMte);
        yPredVTrueLM(2,idx_te) = Y_predict;
    end
catch
% Catch statement to let the code finish running but check back and do not proceed if models have manually assigned PC
    idxTrim = -1;
    Model_LM = fitlm(X_mm(:,1),Y_mm);
    for LOO = 1:length(idxMmFinal)
        % Split test + train
        idx_te = LOO;
        idx_tr = 1:1:length(idxMmFinal); 
        idx_tr(LOO) = [];
        Y_tr = Y_mm(idx_tr);
        X_LMtr = X_mm(idx_tr,1);
        % Build model and predict
        Model_LMtr = fitlm(X_LMtr,Y_tr);
        X_LMte = X_mm(idx_te,1);
        Y_predict = predict(Model_LMtr,X_LMte);
        yPredVTrueLM(2,idx_te) = Y_predict;
    end
end

RMSE_qsquared = sqrt(sum((yPredVTrueLM(1,:) - yPredVTrueLM(2,:)).^2)/length(idxMmFinal));

clear WeightsLASSO_store    
