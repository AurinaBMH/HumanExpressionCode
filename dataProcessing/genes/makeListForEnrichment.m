% make a list of genes to test for enrichment. 
% we're interested if genes selected for the analysis using RNAseq are
% enriched in something

% import data
clear all; 
cd ('data/genes/processedData')
load('MicroarrayDataWITHcustProbesUpdatedXXX.mat')

% remove CUST probes because they can't be used in the enrichment in
% ermineJ
probeName = DataTableProbe.ProbeName{1};
entrezID = DataTableProbe.EntrezID{1}; 
% 
% cust = strfind(probeName, 'CUST');
% remInd = find(~cellfun(@isempty,cust));
% 
% probeName(remInd) = []; 
% entrezID(remInd) = []; 
% noiseall(remInd,:) = []; 

% keep one probe per gene - doesn't matter the rule - all should represent
% the same gene; 

%probeList = probeName(uind); 
%noiseList = noiseall(uind,:); 

% do QC
signalLevel = sum(noiseall,2)./size(noiseall,2);
indKeepProbes = find(signalLevel>=0.5);
indRemProbes = find(signalLevel<0.5);

geneList = zeros(size(noiseall,1),2); 
geneList(:,1) = entrezID; 
geneList(indKeepProbes,2) = 1; 
geneList(indRemProbes,2) = 0.0001; 

% find unique genes and keep one value per gene
[ugenes, uind] = unique(entrezID); 

geneList2 = geneList(uind,:); 
T = table(geneList2(:,1), geneList2(:,2)); 
T.Properties.VariableNames = {'GeneEntrezID' 'Score'}; 

writetable(T, 'QCprobes_enrichment.txt', 'Delimiter', '\t'); 


% probeList(indKeepProbes,2) = {1}; 
% probeList(indRemProbes,2) = {0.0001}; 
% 
% % keep one probe per gene
% 
% 
% 
% T = table(probeList(:,1), probeList(:,2)); 
% T.Properties.VariableNames = {'Probe' 'Score'}; 
% 
% % save table
% writetable(T, 'QCprobes_enrichment.txt', 'Delimiter', '\t'); 
% 
% 
% 
% 
% 
% 
% 
% 
% 
%     fprintf(1,'%d CUST probes removed\n', length(remInd))
%     DataTableProbe.ProbeName{1}(remInd) = {NaN};
%     DataTableProbe.ProbeID{1}(remInd) = NaN;
%     
% probeName = DataTableProbe.ProbeName{1}; 
% probeName(remInd) = [];
% 
% entrezID = DataTableProbe.EntrezID{1}; 
% entrezID(remInd) = []; 
% 





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
