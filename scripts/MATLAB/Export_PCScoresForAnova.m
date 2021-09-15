%% Export_PCScoresForAnova.m

filename = strcat('PCScoresForAnova_',nameHu,nameMm,'.xls');
printTable = table(scoresHu(:,1:nComponents95));
writetable(printTable,filename,'WriteVariableNames',0,'Sheet',1);

printTable = table(scoresm);
writetable(printTable,filename,'WriteVariableNames',0,'Sheet',2);

if isequal(nameHu,'h48350') == 1
    printTable = table(GSE1297_phno.Group,GSE1297_phno.Sex,GSE1297_phno.Age,'VariableNames',{'disease','sex','age'});
elseif isequal(nameHu,'h1297') == 1
    printTable = table(GSE48350_phno.Disease,GSE48350_phno.Sex,GSE48350_phno.Age,'VariableNames',{'disease','sex','age'});
end
writetable(printTable,filename,'Sheet',3);

printTable = table(YMmdisease,YMmage,'VariableNames',{'disease','age'});
writetable(printTable,filename,'Sheet',4);
movefile(filename,'../../modeling outputs');
cd '../../scripts/MATLAB'
