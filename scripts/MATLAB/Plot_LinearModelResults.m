%% Plot_LinearModelResults.m

%%
summaryTable = Model_LM.Coefficients;

rowNames = cell(length(idxTrim),1);
for iter = 1:length(idxTrim)
    rowNames{iter} = strcat('pc',num2str(idxTrim(iter)));
end
rowNames = ['Intercept';rowNames];

summaryTable.Properties.RowNames = rowNames;

%%
Pvalue = [0;modelPFRval(idxTrim,8)];
FDR = [0;modelPFRval(idxTrim,9)];
appendTable = table(Pvalue,FDR);

summaryTable = [summaryTable, appendTable];

clear Pvalue, clear FDR, clear appendTable

%% 
figLMCoefficients = figure;
hold on
plotX = 1:1:length(idxTrim);
bar(plotX,summaryTable.Estimate(2:end));
er = errorbar(plotX,summaryTable.Estimate(2:end),summaryTable.SE(2:end),summaryTable.SE(2:end));
er.Color = [0,0,0];
er.LineStyle = 'none';
rowNames = rowNames';
xticks(plotX)
xticklabels(rowNames(2:end))
ylabel('Linear Model Coefficient')
figLMPValues = figure;
plotX = 1:1:length(idxTrim);
bar(plotX,-log10(summaryTable.Pvalue(2:end)));
xticks(plotX)
xticklabels(rowNames(2:end))
ylabel('-log10(FDR)');
