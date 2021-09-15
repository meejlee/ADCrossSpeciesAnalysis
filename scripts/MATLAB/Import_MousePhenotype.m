%% Import_MousePhenotype.m
% Generate phenotype matrix for GSE64398 samples

%% Phenotype formatting
%%%% Format phenotypes into numeric matrix
phnoMm = zeros(size(GSE64398_phno,1),2);
% 2 rows:
    % 1: Strain
    % 2: Age
for i = 1:length(phnoMm)
    loop_str = GSE64398_phno{i,8};
    if isequal(loop_str,"HET_TASTPM") == 1
        phnoMm(i,1) = 1;
    elseif isequal(loop_str,"HO_TASTPM") == 1
        phnoMm(i,1) = 2;
    elseif isequal(loop_str,"TAU") == 1
        phnoMm(i,1) = 3;    
    elseif isequal(loop_str,"TPM") == 1
        phnoMm(i,1) = 4;    
    elseif isequal(loop_str,"TAS10") == 1
        phnoMm(i,1) = 5;    
    elseif isequal(loop_str,"WILD") == 1
        phnoMm(i,1) = 0;    
    end

    phnoMm(i,2) = GSE64398_phno{i,7};
end

%%%% How many samples/group?
idxMouseGroups = cell(6,1);
counter = 1;
for i = 0:5
    idxMouseGroups{counter} = find(phnoMm(:,1) == i);
    counter = counter + 1;
end