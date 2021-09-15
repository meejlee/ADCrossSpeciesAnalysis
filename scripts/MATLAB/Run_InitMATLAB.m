%% Run_InitMATLAB.m
% Overview: Imports data into MATLAB to start TransComp-R workflow

% Note: Select flags at the start of the code to specify proper human
% cohort

%% Generate additional directories
% (Github doesn't keep track of empty folders)
try
   cd ../../data/interim
catch
   mkdir ../../data/ interim
end
   
try
   cd ../../data/external
catch
   mkdir ../../data/ external
end
   
try
   cd ../../figures
catch
   mkdir ../../ figures
end
   
try
   cd ../../modeling outputs
catch
   mkdir ../../ 'modeling outputs'
end
   
%% Set study flags
%% HUMAN STUDY ID?
% nameHu = 'h48350'; % Cohort 1
% nameHu = 'h1297'; % Cohort 2

%% MOUSE STUDY ID?
nameMm = 'm64398';

%% Import data
run 'Import_HumanData.m'
run 'Import_MouseData.m'
run 'Import_MousePhenotype.m'

%% Mouse-to-Human homolog matching
try
    cd '../data/interim'
    load 'mouse_homologs.mat'
catch
    run 'HomologMatching.m'
end

%% Normalize and 1:1 match genes based on human study flag
if isequal(nameHu,'h48350') == 1
    run 'GeneFiltering_h48350_median.m';
elseif isequal(nameHu,'h1297') == 1
    run 'GeneFiltering_h1297_median.m' 
end
