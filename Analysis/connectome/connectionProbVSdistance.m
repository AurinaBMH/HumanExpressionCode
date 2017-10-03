% calculate connection probability VS distance
type = 'HCP';
doLC = false;

if strcmp(type, 'HCP')
    load('HCPMMP1ANDfslatlas20_acpc_connectome_data.mat')
elseif strcmp(type, 'GenCog')
    load('HCPMMP1ANDfslatlas20_GenCOG_connectome_data.mat');
end
load('HCPMMP1ANDfslatlas20_MNILinear_COGflippedX.mat')

[Adj, consist] = giveMeGroupAdj(standard);
Adj = double(logical(Adj));
dist = pdist2(coordinates,coordinates);

if doLC
    Adj = Adj(1:180,1:180);
    dist = dist(1:180,1:180);
end

BF_PlotQuantiles(dist(:), Adj(:), 50, false, true);
xlabel('Euclidean distance between regions'); ylabel('Connection probability');
if doLC
    title('Left cortex')
else
    title('Whole brain')
end