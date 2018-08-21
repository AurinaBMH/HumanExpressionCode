
function [coexpData, A, matrices, coordinates, avWeight] = giveConnCoexp(parc,tract,weight,brainPart)

if nargin<3
    weight = 'standard';
    brainPart = 'Lcortex'; 
end


if strcmp(parc, 'aparcaseg')
    parcellation = 'aparc+aseg'; 
    LC = 1:34;
    RC = 42:75;
    nROIs = 82;
elseif strcmp(parc, 'HCP')
    parcellation = 'HCPMMP1ANDfslatlas20'; 
    LC = 1:180;
    RC = 191:370;
    nROIs = 360;
elseif strcmp(parc, 'cust100')
    parcellation = 'custom200ANDfslatlas20'; 
    LC = 1:100;
    RC = 111:210;
    nROIs = 220;
end

load(sprintf('%s_acpc_%s_SIFT2_%s_structnets.mat', parcellation, tract, weight));

DS = 100;

% remove subject 299
ADJS{299} = []; ADJS = ADJS(~cellfun('isempty',ADJS));
COG{299} = []; COG = COG(~cellfun('isempty',COG));
SUBS(299) = [];

% stack matrices to 3D matrix
if strcmp(parcellation, 'aparc+aseg')
    LC = 1:34;
    RC = 42:75;
    nROIs = 82;
elseif strcmp(parcellation, 'HCPMMP1ANDfslatlas20')
    LC = 1:180;
    RC = 191:370;
    nROIs = 360;
elseif strcmp(parcellation, 'custom200ANDfslatlas20')
    LC = 1:100;
    RC = 111:210;
    nROIs = 220;
end

coexpData = load(sprintf('%dDS%dscaledRobustSigmoidNGRNAseqQC1Lcortex_ROI_NOdistCorrEuclidean.mat', DS, nROIs)); 

cort = [LC,RC];
if strcmp(brainPart, 'LRcortex')
    keepNodes = cort;
elseif strcmp(brainPart, 'Lcortex')
    keepNodes = LC;
elseif strcmp(brainPart, 'wholeBrain')
    keepNodes = 1:size(ADJS{1},1);
end

numNodes = length(keepNodes);
numSubjects = length(ADJS);
matrices = cell(1,size(ADJS,2));
coordinates = cell(1,size(ADJS,2));
A = zeros(numNodes, numNodes,numSubjects);
for m=1:length(ADJS)
    if ~strcmp(brainPart, 'wholeBrain')
        matrices{m} = ADJS{m}(keepNodes, keepNodes);
        if strcmp(weight, 'standard')
        matrices{m}(matrices{m}<10) = 0;  
        end
        coordinates{m} = COG{m}(keepNodes,:);
        A(:,:,m) = matrices{m};
    else
        matrices{m} = ADJS{m};
        A(:,:,m) = matrices{m};
        coordinates{m} = COG{m}(keepNodes,:);
    end
end

nanA = A;
nanA(nanA==0) = NaN;
avWeight = nanmean(nanA,3);
avWeight(isnan(avWeight)) = 0;
end