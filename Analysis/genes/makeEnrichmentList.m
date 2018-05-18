% make a list of genes with scores for enrichment
load('MicroarrayDataWITHcustProbesUpdatedXXXRNAseqnoQC.mat')

for i=1:size(avgCorr,1)
    c(i) = max(avgCorr{i}); 
end
c=c'; 
A = [genes, c];
A(any(isnan(A), 2),:)=[];

load('MicroarrayDataWITHcustProbesUpdatedXXX.mat')

% get a score for each probe
score = mean(noiseall,2); 
entrezID = DataTableProbe.EntrezID{1}; 
ugene = unique(entrezID, 'stable'); 
scorelist = zeros(length(ugene),2); 
% for each gene select max score value
for i=1:length(ugene)
    ind = find(entrezID==ugene(i)); 
    scorelist(i,1) = ugene(i); 
    scorelist(i,2) = max(score(ind)); 
end


scorelist(any(isnan(scorelist), 2),:)=[];






