%% Plot_HumanPCScores.m
% Plot human patient sample scores in human PC space

%%
xPC = 1;
yPC = 2;
flagLegend = 'disease';

%% Format phenotype input
phnoHuPlot = nan(size(scoresHu,1),1);
if (isequal(flagLegend,'disease') == 1) && isequal(nameHu,'h48350') == 1
    for iter = 1:size(scoresHu,1)
       if (isequal(GSE48350_phno.Disease{iter},'AD') == 1); phnoHuPlot(iter) = 1;
       elseif (isequal(GSE48350_phno.Disease{iter},'Control') == 1); phnoHuPlot(iter) = 0; end 
    end
elseif (isequal(flagLegend,'disease') == 1) && isequal(nameHu,'h1297') == 1
    for iter = 1:size(scoresHu,1)
       if (isequal(GSE1297_phno.Group{iter},'Control') == 1); phnoHuPlot(iter) = 0;
       elseif (isequal(GSE1297_phno.Group{iter},'Incipient') == 1); phnoHuPlot(iter) = 1; 
       elseif (isequal(GSE1297_phno.Group{iter},'Moderate') == 1); phnoHuPlot(iter) = 2; 
       elseif (isequal(GSE1297_phno.Group{iter},'Severe') == 1); phnoHuPlot(iter) = 3; 
       end 
    end
end

%% Plot based on human cohort
figHumanPCScores = figure;
str1 = strcat('PC',num2str(xPC));
str2 = strcat('PC',num2str(yPC));
if isequal(nameHu,'h48350') == 1
    PCAplot(scoresHu(:,xPC),scoresHu(:,yPC),phnoHuPlot,{'ctrl','AD',...
        strcat(str1,'-',num2str(round(varExplainedHu(xPC),1)),'%'),...
        strcat(str2,'-',num2str(round(varExplainedHu(yPC),1)),'%')});
elseif isequal(nameHu,'h1297') == 1
    PCAplot_cat(scoresHu(:,xPC),scoresHu(:,yPC),phnoHuPlot,4,{'Ctrl','Incip','Mod','Sev',...
        strcat(str1,'-',num2str(round(varExplainedHu(xPC),1)),'%'),...
        strcat(str2,'-',num2str(round(varExplainedHu(yPC),1)),'%')});

end
title(strcat(nameHu,': ',str1,'/',str2),'FontSize',16);

%% Functions
function output = PCAplot(PC1,PC2,labels_num,labels_char)
    data0 = zeros(length(labels_num)-nnz(labels_num),2);
    data1 = zeros(nnz(labels_num),2);
    counter0 = 1;
    counter1 = 1;
    for i = 1:length(PC1)
        if labels_num(i) == 0
            data0(counter0,:) = [PC1(i),PC2(i)];
            counter0 = counter0 + 1;
        elseif labels_num(i) == 1
            data1(counter1,:) = [PC1(i),PC2(i)];  
            counter1 = counter1 + 1;
        end
    end
    hold on
    constant_size = 15;
    plot(data0(1,1),data0(1,2),'bo','MarkerFaceColor','b','MarkerSize',constant_size);
    plot(data1(1,1),data1(1,2),'ro','MarkerFaceColor','r','MarkerSize',constant_size);
    leg = legend([labels_char(1),labels_char(2)],'location','best');
    leg.FontSize = 13;
    leg.AutoUpdate = 'off';
    plot(data0(2:end,1),data0(2:end,2),'bo','MarkerFaceColor','b','MarkerSize',constant_size);
    plot(data1(2:end,1),data1(2:end,2),'ro','MarkerFaceColor','r','MarkerSize',constant_size);
    xlabel(labels_char(3),'FontSize',16);
    ylabel(labels_char(4),'FontSize',16);
    ax = gca;
    ax.FontSize = 16;
    output = 0;
end

function output = PCAplot_cat(PC1,PC2,labels_num,noCats,labels_char)    
    red = [1, 0, 0];
    orange = [0.8500, 0.3250, 0.0980];
    yellow = [0.9290, 0.6940, 0.1250];
    blue = [0, 0, 1];
    marksize = 10;
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
    xlabel(labels_char(length(labels_char)-1),'FontSize',16);
    ylabel(labels_char(end),'FontSize',16);
    ax = gca;
    ax.FontSize = 16;
    output = 0;
end
