%% Model_LASSOandLM.m
% Build linear model using consensus PCs from 10 rounds of TransComp-R

%% Fixed inputs to re-generate exact paper figures
% LASSO outcome, esp GSE1297, can slightly shift
if isequal(nameHu,'h48350') == 1
   idxTrim = [2;5;9;11];
elseif isequal(nameHu,'h1297') == 1
    idxTrim = [2;4;5;6;7;8;9;11;14;15;16];
end
%% Linear Model
% Init variables
yPredVTrueLM = nan(2,noSamples);
yPredVTrueLM(1,:) = Y_mm;

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
        
RMSE_qsquared = sqrt(sum((yPredVTrueLM(1,:) - yPredVTrueLM(2,:)).^2)/length(idxMmFinal));
% disp(idxTrim);

cd '../../modeling outputs'
filename = strcat('TransCompR_',nameHu,filenameMmStrain,'.mat');
save(filename,'idxTrim','Model_LM','RMSE_qsquared');
cd '../scripts/MATLAB'