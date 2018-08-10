
load('custom500ANDfslatlas20_acpc_FACT_SIFT2_length_structnets.mat')
DIST = zeros(size(ADJS,2), size(ADJS{1},1), size(ADJS{1},2)); 
for s=1:size(ADJS,2)
    DIST(s,:,:) = ADJS{s}; 
end
 
DISTmean = squeeze(mean(DIST,1)); 

load('custom500ANDfslatlas20_acpc_FACT_SIFT2_standard_structnets.mat')
[Adj, ~] = giveMeGroupAdj(ADJS);

WhatTypeNetwork = 'bu'; 
whatNullModel = 'randmio_und'; 
numIter = 50; 
numRepeats = 100; 

PlotRichClub(Adj,DISTmean, WhatTypeNetwork,whatNullModel,numIter,numRepeats); 
% plot C elegans
load('CelegansData.mat')

WhatTypeNetwork = 'bd'; 
whatNullModel = 'randmio_dir'; 

Adj = C.Adj_B{1,3}; 
DISTmean = C.Eucl_Dist_2D_noLR; 

PlotRichClub(Adj,DISTmean, WhatTypeNetwork,whatNullModel,numIter,numRepeats);
