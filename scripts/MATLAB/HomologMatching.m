%% HomologMatching.m
% Overview: Identify mouse and human gene names for homologous matches

% Notes: Output is a saved map of homologs -- should only need to run once and
% then will just load saved .mat subsequently

%% Homolog matching + saving .mat file
%%%% HOM_MouseHumanSequence.rpt, script, function from Kumar et
%%%% al Cell 2018)
cd '../../data'
mouse_homologs = find_homologs(GSE64398_genename,'mouse, laboratory');

counter = 1;
for i = 1:length(mouse_homologs)
    [geneL,geneW] = size(mouse_homologs{i});
   if geneL == 0
   else
       check_mouse(counter) = i;
       counter = counter + 1;
   end
end

cd '../data/interim'
save('mouse_homologs.mat','mouse_homologs'); 
cd '../../scripts/MATLAB'

%% Function: Find Homologs
function [ return_symbols ] = find_homologs(search_symbols, search_species )
% Find human-mouse homologs
%[ return_ID ] = find_homologs(search_symbols, search_species )
%   INPUTS
%       search_symbols - gene symbols to be converted. 
%       search_species - species of input symbols. Either 'human' or
%       'mouse, laboratory'
%   OUTPUTS: 
%       return_symbols: symbols of identified homologs
    homolog_table = readtable('HOM_MouseHumanSequence.rpt.txt',...
                              'ReadVariableNames',true,...
                              'ReadRowNames',false,...
                              'Delimiter','\t');

    search_organism_idx = strcmp(homolog_table{:,'CommonOrganismName'},search_species);
    find_input = @(input_symbol) ismember(homolog_table{:,'Symbol'},input_symbol) & search_organism_idx;
    get_homolog_ID = @(input_idx) ismember(homolog_table{:,'HomoloGeneID'}, homolog_table{input_idx,'HomoloGeneID'});
    get_homolog_idx = @(homolog_idx)  homolog_idx & ~search_organism_idx ;
    get_homologs = @(return_idx) homolog_table{return_idx,'Symbol'} ;

    return_symbols = cellfun( @(x) get_homologs(get_homolog_idx(get_homolog_ID(find_input(x)))), search_symbols, 'UniformOutput', false );
    
end