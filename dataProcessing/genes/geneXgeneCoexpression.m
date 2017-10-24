% to create gene x gene coexpression matrix:
%% 1. run a first part of S5NormalisationDS script. 
%% Author: Aurina

%Last modiffied: 2017-07-31
%Last modiffied: 2017-08-01
%close all;
%clear all;
%------------------------------------------------------------------------------
% Choose options
%------------------------------------------------------------------------------
useCUSTprobes = true; % choose if you want to use data with CUST probes
probeSelection = 'Variance';% (Variance', LessNoise', 'Mean', 'PC')
parcellation = 'HCP';%, 'cust100', 'cust250'};
distanceThreshold = 30; % first run 30, then with the final threshold 2
percentDS = 100;
multipleProbes = false; % it this is true, only genes that have multiple probes will be selected. 
correctDistance = false; 
distanceCorrection = 'Euclidean';
coexpressionFor = 'all';
Fit = {'removeMean'};
normMethod = 'scaledRobustSigmoid'; %'scaledRobustSigmoid', 'zscore';
normaliseWhat = 'LcortexSubcortex'; %(LcortexSubcortex, wholeBrain, LRcortex)
% choose Lcortex if want to normalise samples assigned to left cortex separately;
% choose LcortexSubcortex if want to normalise LEFT cortex + left subcortex together
% choose wholeBrain if you want to normalise the whole brain.
% choose LRcortex if you want to normalise left cortex + right cortex.

%------------------------------------------------------------------------------
% Define number of subjects and parcellation details based on choises
%------------------------------------------------------------------------------
if strcmp(parcellation, 'aparcaseg')
    NumNodes = 82;
    LeftCortex = 1:34;
    LeftSubcortex = 35:41;
    RightCortex = 42:75;
    RightSubcortex = 76:82;
elseif strcmp(parcellation, 'cust100')
    NumNodes = 220;
    LeftCortex = 1:100;
    LeftSubcortex = 101:110;
    RightCortex = 111:210;
    RightSubcortex = 211:220;
elseif strcmp(parcellation, 'cust250')
    NumNodes = 530;
    LeftCortex = 1:250;
    LeftSubcortex = 251:265;
    RightCortex = 266:515;
    RightSubcortex = 516:530;
    
elseif strcmp(parcellation, 'HCP')
    NumNodes = 360;
    LeftCortex = 1:180;
    %LeftSubcortex = 110;
    RightCortex = 181:360;
    %RightSubcortex = NumNodes;
end

switch normaliseWhat
    case 'Lcortex'
        subjects = 1:6;
        nROIs = LeftCortex;
    case 'LcortexSubcortex'
        subjects = 1:6;
        nROIs = [LeftCortex,LeftSubcortex];
    case 'wholeBrain'
        subjects = 1:2;
        nROIs = 1:NumNodes;
    case 'LRcortex'
        subjects = 1:2;
        nROIs = [LeftCortex,RightCortex];
end

if useCUSTprobes
    startFileName = 'MicroarrayDataWITHcust';
else
    startFileName = 'MicroarrayData';
end

cd ('data/genes/processedData');
load(sprintf('%s%s%dDistThresh%d_CoordsAssigned.mat', startFileName, probeSelection, NumNodes, distanceThreshold));


expressionSubjROI = cell(6,1);
coordinatesSubjROI = cell(6,1);
coordSample = cell(6,1);
expSampNorm = cell(6,1);
expSample = cell(6,1);

entrezIDs = probeInformation.EntrezID; 
load('IDgenes3plus.mat'); 
[~, keep] = intersect(entrezIDs, IDgene); 
%----------------------------------------------------------------------------------
% Normalise data for each subject separately
% Do differential stability calculation:
%1. average expression values for each ROI for differential stability calculation
%2. get DS values for each gene
%3. take top 5% of DS genes
%----------------------------------------------------------------------------------

for sub=subjects
    % normalise data for each subject separately using samples
    expSingleSubj = DataExpression{sub,1};
    coordSingle = DataCoordinatesMNI{sub,1};
    
    switch normaliseWhat
        
        case 'Lcortex'
            
            ind = find(ismember(expSingleSubj(:,2), LeftCortex));
            
        case 'LcortexSubcortex'
            
            LcortexSubcortex = horzcat(LeftCortex,LeftSubcortex);
            ind = find(ismember(expSingleSubj(:,2), LcortexSubcortex));
            
        case 'wholeBrain'
            
            WholeBrain = horzcat(LeftCortex,LeftSubcortex,RightCortex, RightSubcortex);
            ind = find(ismember(expSingleSubj(:,2), WholeBrain));
            
        case 'LRcortex'
            
            LcortexRcortex = horzcat(LeftCortex,RightCortex);
            ind = find(ismember(expSingleSubj(:,2), LcortexRcortex));
            
    end
    expSubj = expSingleSubj(ind,:);
    coord = coordSingle(ind,3:5);
    data = expSubj(:,3:size(expSubj,2));
    data = 2.^(data); 
    if multipleProbes
    data = data(:,keep); 
    end
    %coordSample{sub} = coord;
    ROI = expSubj(:,2);
    % normalise sample x gene data for each subject separately
    %% commented - not normalise
    switch normMethod
        case 'hampel'
            dataNorm = Norm_hampel(data);
            fprintf('Normalising gene expression data\n')
        otherwise
            dataNorm = BF_NormalizeMatrix(data, normMethod);
            fprintf('Normalising gene expression data\n')
    end
    
    expSampNorm{sub} = [ROI, dataNorm];
    coordSample{sub} = [ROI, coord];
    expSample{sub} = [ROI, data];
    ROIs = unique(expSubj(:,2));
    
    % average expression values for each ROI for differential stability calculation
    
    expressionROI = zeros(length(ROIs),size(data,2));
    coordinatesROI = zeros(length(ROIs),3);
    
    for j=1:length(ROIs)
        indROI = find(expSubj(:,2)==(ROIs(j)));
        noProbes = length(indROI);
        fprintf(1,'%u samples for %u ROI found \n', noProbes, ROIs(j))
        % take expression values for a selected entrezID
        expressionRepInt = dataNorm(indROI,:);
        coordinatesRepInt = coord(indROI,:);
        
        % calculate the mean for expression data for a selected entrezID
        expressionROI(j,:)= mean(expressionRepInt,1);
        coordinatesROI(j,:)= mean(coordinatesRepInt,1);
        
    end
    
    S = zeros(length(ROIs),1);
    S(:) = sub;
    expressionSubjROI{sub} = [S, ROIs, expressionROI];
    coordinatesSubjROI{sub} = [S, ROIs, coordinatesROI];
end
%----------------------------------------------------------------------------------
% Combine noramlised data for all subjects.
%----------------------------------------------------------------------------------

expSampNormalisedAll = vertcat(expSampNorm{1}, expSampNorm{2},expSampNorm{3},expSampNorm{4},expSampNorm{5},expSampNorm{6});
combinedCoord = cat(1,coordSample{1}, coordSample{2}, coordSample{3},...
    coordSample{4}, coordSample{5}, coordSample{6});
%% 2. calculate gene x gene coexpression matrix using spearman correlation for all columns starting from second. 
geneCoexpression = corr(expSampNormalisedAll(:,2:end), 'type', 'Spearman'); % calculate sample-sample coexpression
%% 3. Save the matrix and additional information about probes. 
%save('geneXgeneCoexpression.mat', 'geneCoexpression', 'probeInformation');