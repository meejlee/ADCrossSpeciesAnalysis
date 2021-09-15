%% Model_RepeatLASSOLM.m
% Overview: Run LASSO + Linear Modeling TransComp-R workflow multiple times and
% identify PC combination that is selected most frequently

%% Init + Run LASSO
if isequal(nameHu,'h48350') == 1
    noRepeat = 10;
elseif isequal(nameHu,'h1297') == 1
    noRepeat = 50;
end

idxTrim_repeats = cell(noRepeat,1);
lassostates_repeats = nan(noRepeat,2);
for iter = 1:noRepeat
    run 'Model_LASSOandLM.m'
    if iter == 1
        PCFqKfold_repeats = PCFqKfold;
    elseif iter > 1
        PCFqKfold_repeats = [PCFqKfold_repeats,PCFqKfold];
    end    
    idxTrim_repeats{iter} = idxTrim;
    lassostates_repeats(iter,1) = length(idxTrim);
    lassostates_repeats(iter,2) = Model_LM.Rsquared.Adjusted;
end
PCFqKfold_summary = sum(PCFqKfold_repeats,2);
idxTrim = find(PCFqKfold_summary > nStab*noRepeat*nfold*nThresh);

%% Plot frequency of selection
figLASSOselecfq = figure;
heatmap(PCFqKfold_repeats);
xlabel(strcat('CV folds, ',num2str(noRepeat),' repeats'));
ylabel('PC number')

%% Save figure
cd '../../figures'
filename = strcat('LASSOselecfq_',nameHu,filenameMmStrain,'.fig');
saveas(figLASSOselecfq,filename);
close(figLASSOselecfq);

%% Save outcomes
cd  '../modeling outputs'
filename = strcat('LASSORepeats_',nameHu,filenameMmStrain,'.mat'); 
save(filename,'PCFqKfold_repeats','lassostates_repeats','idxTrim');

cd '../scripts/MATLAB'
