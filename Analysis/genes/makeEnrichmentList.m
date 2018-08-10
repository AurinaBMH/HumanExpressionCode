% make a list of genes with scores for enrichment
load('MicroarrayDataWITHcustProbesUpdatedXXXRNAseqnoQC.mat')

for i=1:size(avgCorr,1)
    c(i) = max(avgCorr{i}); 
end
c=c'; 
RNAseqCorr = [genes, c];
RNAseqCorr(any(isnan(RNAseqCorr), 2),:)=[]; % this is used for GSR

% make a list for ORA for low correlated genes
corVals = RNAseqCorr(:,2); 
indRNA = corVals <0.2; 
RNAseqORA = zeros(length(indRNA),1); 
RNAseqORA(indRNA==1) = 0.0001; 
RNAseqORA(indRNA==0) = 1; 
ORAlistgenes = [RNAseqCorr(:,1), RNAseqORA]; % this is used for ORA


load('MicroarrayDataWITHcustProbesUpdatedXXX.mat')

% get a score for each probe
score = mean(noiseall,2); 
entrezID = DataTableProbe.EntrezID{1}; 
ugene = unique(entrezID, 'stable'); 
backgroundGSR = zeros(length(ugene),2); 
% for each gene select max score value
for i=1:length(ugene)
    ind = find(entrezID==ugene(i)); 
    backgroundGSR(i,1) = ugene(i); 
    backgroundGSR(i,2) = max(score(ind)); 
end
backgroundGSR(any(isnan(backgroundGSR), 2),:)=[];

%% make a list of gene for enrichment
% testing if genes that don't overlap between RNA-seq and microarray data
% are overrepresented in any functional groups. 
% load data selected based on RNAseq
load('MicroarrayDataWITHcustProbesUpdatedXXXRNAseqnoQC.mat')
% in av Corr file NaN values correspond to gene that do not overlap between
% microarray and RNAseq files; make a list for enrichment where selected
% genes are the ones excluded
for i=1:size(avgCorr,1)
    c(i) = max(avgCorr{i}); 
end
c=c'; 

RNAseqCorr = [genes, c];
ind = isnan(c); 
RNAseqCorr(ind,2) = 0.0001; 
RNAseqCorr(~ind,2) = 1;










