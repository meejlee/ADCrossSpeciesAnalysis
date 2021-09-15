%% Generate linear model prediction
Xpred_mm = predict(Model_LM,X_mm(:,idxTrim));

%% Format phenotype input
phnoMmPlot = phnoMm(idxMmFinal,1);

%% Plot based on human cohort
orange = [0.8500, 0.3250, 0.0980];
yellow = [0.9290, 0.6940, 0.1250];
blue = [0, 0, 1];

figTransCompRPredictions = figure;
hold on
for iter = 0:2
    if (iter == 0); p = blue;
    elseif (iter == 1); p = yellow;
    elseif (iter == 2); p = orange;
    end
    idxPlot = find(phnoMmPlot == iter);
    scatter(phnoMmPlot(idxPlot),Xpred_mm(idxPlot),100,p,'filled');
end