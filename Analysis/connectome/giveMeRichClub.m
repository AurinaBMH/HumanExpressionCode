% Rich club in the whole brain

function [G,  PhiNormMean, dMiddle] = giveMeRichClub(matrices, COG, groupMatrixType, densThreshold, giveRC, numIter, numRepeats, whatTypeNetwork ,whatNullModel,scaling)

if nargin < 3
    groupMatrixType = 'consistency';
    densThreshold = 0.6; 
    fprintf('Making consistency based group matrix BY DEFAULT\n')
end

if nargin < 6
    numIter = 50;
    numRepeats = 100;
    whatTypeNetwork = 'bu';
    whatNullModel = 'randmio_und';
    scaling = 0.4;
    if giveRC
    fprintf('Calculating RC using %d iterations with %d repeats on %s network with %s null model with %d distance scaling BY DEFAULT\n', numIter, numRepeats, whatTypeNetwork, whatNullModel, scaling)
    end
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
    if numNodes==180
        hemiid = ones(numNodes,1);
    else
        
    hemiid = zeros(numNodes,1);
    hemiid(1:numNodes/2) = 1;
    hemiid(numNodes/2+1: numNodes) = 2;
    end
    
    Gr = fcn_group_average(A,mdist,hemiid);
    G = Gr.*avWeight;
elseif strcmp(groupMatrixType, 'variance')
    Gr = giveMeGroupAdj_variance(matrices, densThreshold);
    Gr = logical(Gr); 
    G = Gr.*avWeight;
    
elseif strcmp(groupMatrixType, 'varianceAA')
    Gr = giveMeGroupAdj_variance_AA(matrices, densThreshold);
    Gr = logical(Gr); 
    G = Gr.*avWeight;

elseif strcmp(groupMatrixType, 'consistency')
    G = giveMeGroupAdj_consistency(matrices, densThreshold);
end

if giveRC

[dMiddleNorm, dMiddle, PhiNormMean, PhiTrue, PhiRand] = PlotRichClub(G,mdist,whatTypeNetwork,whatNullModel,numIter,numRepeats,scaling);
end

end