
parcellation = 'HCPMMP1ANDfslatlas20';
tract = 'FACT';
weight = 'standard';
brainPart = 'Lcortex'; % Lcortex; 

load(sprintf('%s_acpc_%s_SIFT2_%s_structnets.mat', parcellation, tract, weight));
% remove subject 299
ADJS{299} = []; ADJS = ADJS(~cellfun('isempty',ADJS));
COG{299} = []; COG = COG(~cellfun('isempty',COG));
SUBS(299) = [];

% stack matrices to 3D matrix
if strcmp(parcellation, 'aparc+aseg')
    LC = 1:34;
    RC = 42:75;
elseif strcmp(parcellation, 'HCPMMP1ANDfslatlas20')
    LC = 1:180;
    RC = 191:370;
elseif strcmp(parcellation, 'custom200ANDfslatlas20')
    LC = 1:100;
    RC = 111:210;
end

cort = [LC,RC];

if strcmp(brainPart, 'LRcortex')
 keepNodes = cort; 
elseif strcmp(brainPart, 'Lcortex')
 keepNodes = LC; 
elseif strcmp(brainPart, 'wholeBrain')
 keepNodes = 1:size(ADJS{1},1);
end

numNodes = length(keepNodes); 
matrices = cell(1,size(ADJS,2)); 
coordinates = cell(1,size(ADJS,2)); 
for m=1:length(ADJS)
    if ~strcmp(brainPart, 'wholeBrain')
        matrices{m} = ADJS{m}(keepNodes, keepNodes);
        coordinates{m} = COG{m}(keepNodes,:); 
    else
        matrices{m} = ADJS{m}; 
    end
end

for th=[0.1 0.15 0.2 0.25 0.3]

[PhiNormMean, G] = giveMeRichClub(matrices, coordinates, 'variance', th);
figure; imagesc(log(G)); 
end


[PhiNormMean, G] = giveMeRichClub(matrices, coordinates, 'lengthCV');
figure; imagesc(log(G));
[PhiNormMean, G] = giveMeRichClub(matrices, coordinates, 'consistency');
figure; imagesc(log(G));
