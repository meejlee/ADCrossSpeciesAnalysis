%% Model_NullModelPCs.m
% Generate 1000 (or max possible #) random/null models with size-matched PCs to calculate an empirical p-value for the true linear rmodel

%% TEST 1: Random PCs
%% Init variables
scoresm = X_mm;
max_iter = nComponents95;
flag_nopcs = length(idxTrim);
noSamples = length(idxMmFinal);

%% Run null testing
if flag_nopcs == 1
    % Init
    rmse_NullPCs = zeros(max_iter,1);
    % Null models
    for iter = 1:max_iter
        yPredictNull = nan(length(idxMmFinal),1);
        for LOO = 1:length(idxMmFinal)
            % Split test + train
            idx_te = LOO;
            idx_tr = 1:1:length(idxMmFinal); 
            idx_tr(LOO) = [];
            Y_tr = Y_mm(idx_tr);
            % Build model and predict
            mdlNull = fitlm(scoresm(idx_tr,iter),Y_mm(idx_tr));
            yPredictNull(LOO) = predict(mdlNull,scoresm(idx_te,iter));
        end
        rmse_NullPCs(iter) = sqrt(sum((Y_mm - yPredictNull).^2)/length(idxMmFinal));
    end
elseif flag_nopcs == 2
    counter = 1;
    for iter = 1:max_iter
        for iter2 = 1:max_iter
            if (isequal(iter,iter2) == 1); break; end
            yPredictNull = nan(length(idxMmFinal),1);
            for LOO = 1:length(idxMmFinal)
                % Split test + train
                idx_te = LOO;
                idx_tr = 1:1:length(idxMmFinal); 
                idx_tr(LOO) = [];
                Y_tr = Y_mm(idx_tr);
                % Build model and predict
                mdlNull = fitlm(scoresm(idx_tr,[iter,iter2]),Y_mm(idx_tr));
                yPredictNull(LOO) = predict(mdlNull,scoresm(idx_te,[iter,iter2]));
            end
            rmse_NullPCs(counter) = sqrt(sum((Y_mm - yPredictNull).^2)/length(idxMmFinal));
            counter = counter + 1;
        end
    end
elseif flag_nopcs == 3
    counter = 1;
    for iter = 1:max_iter
        for iter2 = 1:max_iter
            for iter3 = 1:max_iter
                if (isequal(iter,iter3) == 1); break; end
                if (isequal(iter,iter2) == 1); break; end
                if (isequal(iter2,iter3) == 1); break; end
                yPredictNull = nan(length(idxMmFinal),1);
                for LOO = 1:length(idxMmFinal)
                    % Split test + train
                    idx_te = LOO;
                    idx_tr = 1:1:length(idxMmFinal); 
                    idx_tr(LOO) = [];
                    Y_tr = Y_mm(idx_tr);
                    % Build model and predict
                    mdlNull = fitlm(scoresm(idx_tr,[iter,iter2,iter3]),Y_mm(idx_tr));
                    yPredictNull(LOO) = predict(mdlNull,scoresm(idx_te,[iter,iter2,iter3]));
                end                    
                rmse_NullPCs(counter) = sqrt(sum((Y_mm - yPredictNull).^2)/length(idxMmFinal));
                counter = counter + 1;
            end
        end
    end
else %i.e., the number of PCs selected is greater than three and if the original 
    % # of PCs is greater than 10, then looping through all possible combinations 
    % should be greater than 1000 random combos
    max_iter = 1000;
    pcs_scram = zeros(max_iter,flag_nopcs);
    rmse_NullPCs = zeros(max_iter,1);
    for iter = 1:max_iter
        fold_flag = 1;
        while fold_flag > 0
            pcs_fold = randperm(nComponents95);
            pcs_fold = pcs_fold(1:flag_nopcs);
            if sum(ismember(pcs_scram,pcs_fold,'rows')) == 0
                pcs_scram(iter,:) = pcs_fold;
                fold_flag = -1;
            end 
        end

        yPredVTrueLM = nan(2,noSamples);
        yPredVTrueLM(1,:) = Y_mm;
        for LOO = 1:length(idxMmFinal)
            % Split test + train
            idx_te = LOO;
            idx_tr = 1:1:length(idxMmFinal); 
            idx_tr(LOO) = [];
            Y_tr = Y_mm(idx_tr);
            % Build model and predict
            Model_LMtr = fitlm(scoresm(idx_tr,pcs_fold),Y_tr);
            yPredVTrueLM(2,idx_te) = predict(Model_LMtr,scoresm(idx_te,pcs_fold));
        end
        rmse_NullPCs(iter) = sqrt(sum((yPredVTrueLM(1,:) - yPredVTrueLM(2,:)).^2)/length(idxMmFinal));
    end
end

[rmse_NullPCsSorted,~] = sort(rmse_NullPCs,'ascend');
rmse_NullPCsSorted = rmse_NullPCsSorted';
