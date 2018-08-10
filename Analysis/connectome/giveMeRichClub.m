% Rich club in the whole brain

function [PhiNormMean, G] = giveMeRichClub(matrices, COG, groupMatrixType,densThreshold,numIter, numRepeats, whatTypeNetwork ,whatNullModel,scaling)

if nargin < 3
    groupMatrixType = 'consistency';
    densThreshold = 0.3; 
    fprintf('Making consistency based group matrix BY DEFAULT\n')
end

if nargin < 5
    numIter = 50;
    numRepeats = 100;
    whatTypeNetwork = 'bu';
    whatNullModel = 'randmio_und';
    scaling = 0.4;
    fprintf('Calculating RC using %d iterations with %d repeats on %s network with %s null model with %d distance scaling BY DEFAULT\n', numIter, numRepeats, whatTypeNetwork, whatNullModel, scaling)
end

numNodes = size(matrices{1},1);
numSubjects = length(matrices);

A = zeros(numNodes, numNodes,numSubjects);
dist = zeros(numNodes, numNodes,numSubjects);
for m=1:numSubjects
    A(:,:,m) = matrices{m};
    coords = COG{m};
    dist(:,:,m) = pdist2(coords,coords);
end

mdist = mean(dist,3);
% get the weights for only the existing edges
nanA = A;
nanA(nanA==0) = NaN;
avWeight = nanmean(nanA,3);
avWeight(isnan(avWeight)) = 0;

% make a group matrix
if strcmp(groupMatrixType, 'lengthCV')
    hemiid = zeros(numNodes,1);
    hemiid(1:numNodes/2) = 1;
    hemiid(numNodes/2+1: numNodes) = 2;
    
    Gr = fcn_group_average(A,mdist,hemiid);
    G = Gr.*avWeight;
elseif strcmp(groupMatrixType, 'variance')
    Gr = giveMeGroupAdj_variance(matrices, densThreshold);
    Gr = logical(Gr); 
    G = Gr.*avWeight;

elseif strcmp(groupMatrixType, 'consistency')
    G = giveMeGroupAdj_consistency(matrices, densThreshold);
end

[dMiddleNorm, dMiddle, PhiNormMean, PhiTrue, PhiRand] = PlotRichClub(G,mdist,whatTypeNetwork,whatNullModel,numIter,numRepeats,scaling);

end