%% Run_Plotting.m
% Plot results from TransComp-R and age vs disease linear modeling
% Re-generate figures in manuscript

%% Calculate Mouse Variance
varExplained_Mm = var(X_mm);
varExplained_Mm = 100*varExplained_Mm'/sum(varExplained_Mm);
varExplained_combined = horzcat(transpose(1:1:nComponents95),varExplainedHu(1:nComponents95),varExplained_Mm);
corrMm = corr(X_mm);

%% Plot
figHuMmVariance = figure;
bar(varExplained_combined(:,2:3));
xticks(varExplained_combined(:,1));
legend('Human','Mouse');
xlabel('Human PC Number');
ylabel('Percent Variance Explained');
ax = gca; ax.FontSize = 20;

%% Generate Plots
run 'Plot_TransCompRPredictions.m'
run 'Plot_HumanPCScores.m'
run 'Plot_MousePCScores.m'
run 'Plot_LinearModelResults.m'

%% Save figures
cd '../../figures'

filename = strcat('TransCompRPredictions_',nameHu,filenameMmStrain,'.fig');
saveas(figTransCompRPredictions,filename);
filename = strcat('HumanPCScores_',nameHu,filenameMmStrain,'.fig');
saveas(figHumanPCScores,filename);
filename = strcat('MousePCScores_',nameHu,filenameMmStrain,'.fig');
saveas(figMousePCScores,filename);
filename = strcat('TransCompRModelCoefficients_',nameHu,filenameMmStrain,'.fig');
saveas(figLMCoefficients,filename);
filename = strcat('AgevsDiseaseLMPValues_',nameHu,filenameMmStrain,'.fig');
saveas(figLMPValues,filename);
filename = strcat('HumanAndMouseVarianceExpl_',nameHu,filenameMmStrain,'.fig');
saveas(figHuMmVariance,filename);
close all
cd '../scripts/MATLAB'

%% Generate xls of PC scores for Anova
run 'Export_PCScoresForAnova.m'