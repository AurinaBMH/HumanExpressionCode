% make a list of genes to test for enrichment. 
% we're interested if genes selected for the analysis using RNAseq are
% enriched in something
clear all; 

load('MicroarrayDataProbesUpdated.mat')
load('MicroarrayDataProbesUpdatedRNAseq.mat')%- %last version kept only genes that were above cackground (lessnoise)
% and were correlated with RNAseq >0.3
[A, ind1] = intersect(DataTableProbe.ProbeName{1}, probeInformation.ProbeName);
list = DataTableProbe.ProbeName{1}; 
% label genes that were selected with a small number. Keep others as
% background. 
for i=1:length(ind1)
list{ind1(i),2} = 0.001; 
end

% save this file as txt and replace [] with 1. Make sure to copy this to
% atom and see if formating is right. 
% Then run enrichment using ermineJ
