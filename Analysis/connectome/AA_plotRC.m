
parcellation = 'HCPMMP1ANDfslatlas20';
tract = 'FACT';
weight = 'density'; % density, FA
brainPart = 'Lcortex'; % Lcortex; LRcortex; wholeBrain; 

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
% 
for th=[0.1 0.15 0.2 0.25 0.3]

[~, G] = giveMeRichClub(matrices, coordinates, 'variance', th);
figure; imagesc(log(G)); 
end

[~, G] = giveMeRichClub(matrices, coordinates, 'lengthCV');
figure; imagesc(log(G));
[PhiNormMean, G] = giveMeRichClub(matrices, coordinates, 'consistency');
figure; imagesc(log(G));

load('10DS360scaledRobustSigmoidNSGRNAseqQC1Lcortex_ROI_NOdistCorrEuclidean.mat')

[~, G] = giveMeRichClub(matrices, coordinates, 'variance', 0.25);
Gr = logical(G);
RichClubHuman(G,averageCoexpression);

% see if there is the trlationship between CGE and distance for hubs
deg = degrees_und(Gr); 
isHub = deg>60; 

A = zeros(numNodes, numNodes); 
A(isHub, isHub) = 1; 
A = A.*Gr;
figure; imagesc(A); 
B = averageCoexpression.*A; 
figure; imagesc(B); 
di = pdist2(COG{1}, COG{1}); 
di = di(1:180,1:180); 
C = di.*A; figure; imagesc(C); 
C(C==0) = NaN; B(B==0) = NaN; 
distVect = C(:); coexpVect = B(:); 
distVect(isnan(distVect)) = []; 
coexpVect(isnan(coexpVect)) = [];
figure; scatter(distVect, coexpVect); 


% see CGE-distance relationship for all connected regions
CGEcon = averageCoexpression.*Gr; 
DISTcon = di.*Gr; 
CGEcon(CGEcon==0) = NaN; DISTcon(DISTcon==0) = NaN; 
distExpVect(:,1) = DISTcon(:); 
distExpVect(:,2) = CGEcon(:);

distExpVect(any(isnan(distExpVect), 2),:)=[];
figure; scatter(distExpVect(:,1),distExpVect(:,2)); 




