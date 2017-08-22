%% Author: Aurina

%Last modiffied: 2017-07-31
%Last modiffied: 2017-08-22
%close all;
%clear all;
%------------------------------------------------------------------------------
% Choose options
%------------------------------------------------------------------------------
useCUSTprobes = false; % choose if you want to use data with CUST probes
probeSelection = 'PC';% (Variance', LessNoise', 'Mean', 'PC')
parcellation = 'aparcaseg';%, 'cust100', 'cust250'};
distanceThreshold = 2; % first run 30, then with the final threshold 2
normMethod = 'zscore';
normaliseWhat = 'Lcortex'; %(LcortexSubcortex, wholeBrain, LRcortex)
% choose Lcortex if want to normalise samples assigned to left cortex separately;
% choose LcortexSubcortex if want to normalise LEFT cortex + left subcortex together
% choose wholeBrain if you want to normalise the whole brain.
% choose LRcortex if you want to normalise left cortex + right cortex.
%------------------------------------------------------------------------------
% Assign variables
%------------------------------------------------------------------------------
if strcmp(parcellation, 'aparcaseg')
    NumNodes = 82;
elseif strcmp(parcellation, 'cust100')
    NumNodes = 220;
elseif strcmp(parcellation, 'cust250')
    NumNodes = 530;
end

if useCUSTprobes
    startFileName = 'MicroarrayDataWITHcust';
else
    startFileName = 'MicroarrayData';
end

cd ('data/genes/processedData');
load(sprintf('%s%s%dDistThresh%d_CoordsAssigned.mat', startFileName, probeSelection, NumNodes, distanceThreshold));
%------------------------------------------------------------------------------
% Generate 1000 distributions for DS for each gene
%------------------------------------------------------------------------------
doRandom = true; 
numGenes = size(DataExpression{1,1},2)-2; 
DSall = zeros(1000,numGenes); 
for i=1:1000
       
DSall(i,:) = calculateDS(DataExpression,DataCoordinatesMNI, parcellation, normaliseWhat, normMethod, doRandom); 

end
%------------------------------------------------------------------------------
% Generate empirical DS scores 
%------------------------------------------------------------------------------
DSempirical = calculateDS(DataExpression,DataCoordinatesMNI, parcellation, normaliseWhat, normMethod); 
Pvals = zeros(numGenes,1); 
Tvals = zeros(numGenes,1); 
%------------------------------------------------------------------------------
% Evaluate DS value for each gene by performing ttest
%------------------------------------------------------------------------------
for gene = 1:numGenes
    
[h,p,ci,stats] = ttest2(DSempirical(gene),DSall(:,gene));
Pvals(gene) = p; 
Tvals(gene) = stats.tstat; 
    
end
