%% Plot_MousePCScores.m

%% PCs?
% Note: Hard coded based on results in manuscript to explicitly re-generate
% supplementary figures
if isequal(nameHu,'h48350') == 1
    PCs_all = [2,5;2,9;2,11;5,9;5,11;9,11];
    sbpltR = 2; sbpltC = 3;
elseif isequal(nameHu,'h1297') == 1
    PCs_all = [2,4;2,8;2,9;2,16;4,8;4,9;4,16;8,9;8,16;9,16];
    sbpltR = 2; sbpltC = 5;
end
flagLegend = 'disease';

%% Plotting
figMousePCScores = figure;

for iter = 1:length(PCs_all)
    subplot(sbpltR,sbpltC,iter)
    PCs = PCs_all(iter,:);
    %% Format phenotype input
    if isequal(flagLegend,'disease') == 1
        phnoMmPlot = phnoMm(idxMmFinal,1);
    end

    %% Plot based on human cohort
    str1 = strcat('PC',num2str(PCs(1)));
    str2 = strcat('PC',num2str(PCs(2)));
    PCAplot_cat(X_mm(:,PCs(1)),X_mm(:,PCs(2)),phnoMmPlot,3,{'WT','HET','HO',...
        strcat(str1,'-',num2str(round(varExplained_combined(PCs(1),3))),'%'),...
        strcat(str2,'-',num2str(round(varExplained_combined(PCs(2),3))),'%')});
    title(strcat(str1,'/',str2),'FontSize',12);
end

%% Functions
function output = PCAplot_cat(PC1,PC2,labels_num,noCats,labels_char)    
    red = [1, 0, 0];
    orange = [0.8500, 0.3250, 0.0980];
    yellow = [0.9290, 0.6940, 0.1250];
    blue = [0, 0, 1];
    marksize = 5;
    fontsize = 10;
    
    hold on
    h1 = plot(1,1,'o','MarkerEdgeColor','k','MarkerFaceColor',blue,'MarkerSize',marksize);
    h2 = plot(1,1,'o','MarkerEdgeColor','k','MarkerFaceColor',yellow,'MarkerSize',marksize);
    h3 = plot(1,1,'o','MarkerEdgeColor','k','MarkerFaceColor',orange,'MarkerSize',marksize);
    if noCats == 3
        leg = legend([labels_char(1),labels_char(2),labels_char(3)],'location','best');
    elseif noCats == 4
        h4 = plot(1,1,'o','MarkerEdgeColor','k','MarkerFaceColor',red,'MarkerSize',marksize);
        leg = legend([labels_char(1),labels_char(2),labels_char(3),labels_char(4)],'location','best');
    end
    leg.FontSize = 13; 
    leg.AutoUpdate = 'off';
    delete(h1); delete(h2); delete(h3); 
    if (noCats == 4); delete(h4); end
    for i = 1:length(PC1)
        if labels_num(i) == 3
            p = red;
        elseif labels_num(i) == 2
            p = orange;
        elseif labels_num(i) == 1
            p = yellow;
        elseif labels_num(i) == 0
            p = blue;
        end
        plot(PC1(i),PC2(i),'o','MarkerEdgeColor',p,'MarkerFaceColor',p,'MarkerSize',marksize);
    end
    xlabel(labels_char(length(labels_char)-1),'FontSize',fontsize);
    ylabel(labels_char(end),'FontSize',fontsize);
    ax = gca;
    ax.FontSize = fontsize;
    output = 0;
end
