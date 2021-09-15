%% Model_NullModelTesting.m
% Test true model against null models

%% TEST 2: Scrambled phenotype
%% Initialize permuted phenotypes
max_null = 100;
nThresh = 0.4;
null_perms = zeros(max_null,length(Y_mm));
for i = 1:1.5*max_null %1.5 buffer in case permutations overlap/aren't unique
    null_perms(i,:) = Y_mm(randperm(length(Y_mm)));
end
null_perms = unique(null_perms,'rows');

rmse_NullPheno = zeros(max_null,2); %intialize to score

for iter_null = 1:max_null
    iter_phno = null_perms(iter_null,:); %pull scrambled Y vector

    %% LASSO
    alphaFinal = 1.0;
    nfold = 5;    

    yPredVTrue_null = nan(2,noSamples);
    PCFqCumulative_null = nan(nComponents95,1);

    CVO_null = cvpartition(iter_phno,'kfold',nfold);
    WeightsLASSO_store_null = cell(nfold,1);

    for iterCV_null = 1:nfold
        % Split test + train
        idx_te_null = find(CVO_null.test(iterCV_null)==1);
        idx_tr_null = find(CVO_null.training(iterCV_null)==1);

        phno_tr_null = iter_phno(idx_tr_null);
        X_LASSOtr_null = X_mm(idx_tr_null,:);

        nStab=10; ncut=nThresh*nStab; Sel_matrix_null=zeros(nStab,nComponents95);
        WeightsLASSO_loop_null = zeros(nComponents95,nStab);
        for k=1:nStab
            [B_null, FitInfo_null] = lasso(X_LASSOtr_null, phno_tr_null,'CV',5,'Alpha',alphaFinal); %Fold-specific variable selection using LASSO using ONLY fold-specific training data. Lambda tuned using an internal CV on the fold-specific training data ONLY.
            chosen_b_null=max(FitInfo_null.Index1SE,1);
            WeightsLASSO_loop_null(:,k) = B_null(:,chosen_b_null);
            B_fin_null=B_null(:,chosen_b_null); var_selec_pre=find(abs(B_fin_null)>0); Sel_matrix_null(k,var_selec_pre)=1;
        end

        % Store LASSO weights for each of the 10 repeats
        WeightsLASSO_store_null{iterCV_null} = WeightsLASSO_loop_null;
    end

    %%%% Storing info from multiple rounds of LASSO CV
    % Pull weights from all LASSO runs (nStab times per nfold CV folds)
    % Reformatting to pull from cell array into a single matrix
    WeightsLASSO_null = nan(nComponents95,nfold*nStab);
    counter = 1;
    for iterCV_null = 1:nfold
        matLoop = WeightsLASSO_store_null{iterCV_null};
        idxloop = (counter):1:(counter+nStab-1);
        WeightsLASSO_null(:,idxloop) = matLoop;
        counter = counter + nStab;
    end

    % Frequency of PC selection across nStab * nfold # of LASSO runs
    for i = 1:nComponents95
        PCFqCumulative_null(i,1) = nnz(WeightsLASSO_null(i,:));
    end

    %% Linear Model
    %% Final Model with all samples
    % Pull PCs that were selected 40% + of the time for the final model
    counter = 1;
    for m = 1:nComponents95
       if PCFqCumulative_null(m) >= (nfold*nStab*nThresh)
           idxTrim_null(counter) = m;
           counter = counter + 1;
       end
    end

    % Build linear model
    try
        yPredict_null = nan(length(idxFinal),1);
        for LOO = 1:length(idxFinal)
            % Split test + train
            idx_te_null = LOO;
            idx_tr_null = 1:1:length(idxFinal); 
            idx_tr_null(LOO) = [];
            Y_tr_null = Y_mm(idx_tr_null);
            X_LMtr_null = X_mm(idx_tr_null,idxTrim_null);
            % Build model and predict
            Model_LMtr_null = fitlm(X_LMtr_null,Y_tr_null);
            X_LMte_null = X_mm(idx_te_null,idxTrim_null);
            yPredict_null(LOO) = predict(Model_LMtr_null,X_LMte_null);
        end
        rmse_NullPheno(iter_null) = sqrt(sum((iter_phno - yPredict_null).^2)/length(idxFinal));
    catch
        idxTrim_null = -1;
        rmse_NullPheno(iter_null,2) = -1;
    end
    clear idxTrim_null
    disp(iter_null);
end

rmse_NullPhenoSorted = sort(rmse_NullPheno,'ascend');

%% Save results as .mat file
cd '../../modeling outputs'
empiricalPVal = length(find(rmse_NullPCsSorted < RMSE_qsquared));
filename = strcat('NullModelTesting_',nameHu,filenameMmStrain,'.mat');
save(filename,'rmse_NullPCsSorted','rmse_NullPhenoSorted','empiricalPVal');
cd '../scripts/MATLAB'
