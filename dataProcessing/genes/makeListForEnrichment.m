% make a list of genes to test for enrichment. 
% we're interested if genes selected for the analysis using RNAseq are
% enriched in something
clear all; 
cd ('data/genes/processedData');
load('MicroarrayDataProbesUpdated.mat')
load('MicroarrayDataProbesUpdatedRNAseq.mat')%- %last version kept only genes that were above cackground (lessnoise)
% and were correlated with RNAseq >0.3

%[genes, indGenes] = unique(DataTableProbe.EntrezID{1}); 
%listID = DataTableProbe.EntrezID{1}(indGenes); 
%listName = DataTableProbe.ProbeName{1}(indGenes); 

[A, ind1, ind2] = intersect(genes, probeInformationALL.EntrezID);
listName = probeInformationALL.ProbeName(ind2); 
% label genes that were selected with a small number. Keep others as
% background. 
for i=1:length(ind1)
listName{i,2} = max(avgCorr{i}); 
end


listName(any(cellfun(@(x) any(isnan(x)),listName),2),:) = [];
% save this file as txt and replace [] with 1. Make sure to copy this to
% atom and see if formating is right. 
% Then run enrichment using ermineJ
