%% Model_Linear_AgeVSDisease.m
% Overview: Generate linear models for AGE vs (AGE + DISEASE) as factors on a
% pc-by-pc basis

%% 
modelPFRval = nan(nComponents95,9);

for i = 1:nComponents95
    % reformat data for fitlm function
    tblNull = table(YMmage,X_mm(:,i),'VariableNames',{'Age','PC'});
    tblAlt = table(YMmage,YMmdisease,X_mm(:,i),'VariableNames',{'Age','Disease','PC'});
    % Null and alternate models
    mdlNull = fitlm(tblNull);
    mdlAlt = fitlm(tblAlt);
    
    statsNull = anova(mdlNull,'summary');
    statsAlt = anova(mdlAlt,'summary');

    % Store/format outputs for visualization
    modelPFRval(i,1) = statsNull{2,4}; % Fstat
    modelPFRval(i,2) = statsNull{2,5}; %pVal
    modelPFRval(i,3) = mdlNull.Rsquared.Ordinary;

    modelPFRval(i,4) = statsAlt{2,4};
    modelPFRval(i,5) = statsAlt{2,5};
    modelPFRval(i,6) = mdlAlt.Rsquared.Ordinary;

    modelPFRval(i,7) = ((statsNull{3,1} - statsAlt{3,1})/(statsAlt{3,1}/(length(X_mm)-size(tblAlt,2)-1)));
    modelPFRval(i,8) = 1-fcdf(modelPFRval(i,7),1,length(X_mm));        
end
modelPFRval(:,9) = round(mafdr(modelPFRval(:,8),'BHFDR','true'),4);

cd '../../modeling outputs'
filename = strcat('LMagevsdis_',nameHu,filenameMmStrain,'.mat');
save(filename,'modelPFRval');
cd '../scripts/MATLAB'