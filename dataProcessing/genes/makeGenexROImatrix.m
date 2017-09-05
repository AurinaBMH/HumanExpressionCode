% make ROI x gene matrix
useCUSTprobes = true;
% choose what type of probe selection to use, hemisphere, subject list, parcellations, threshols.
probeSelection = 'Variance';% (Variance', LessNoise', 'Mean')
parcellation = {'cust250'};%, aparcaseg, 'cust100', 'cust250'};
distanceThreshold = 2; % first run 30, then with the final threshold 2
subjects = 1:6;
normaliseWhat = 'Lcortex'; % 'Lcortex'; % 'LcortexSubcortex'
normMethod = 'scaledRobustSigmoid';
percentDS = 5;


if useCUSTprobes
    fprintf(1,'Using the data with CUST probes: %s \n', probeSelection)
    startFileName = 'MicroarrayDataWITHcust';
else
    fprintf(1,'Using the data without CUST probes: %s \n', probeSelection)
    startFileName = 'MicroarrayData';
end

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
    
end

cd ('data/genes/processedData')
load(sprintf('%s%s%dDistThresh%d_CoordsAssigned.mat', startFileName, probeSelection, NumNodes, distanceThreshold));

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
    % transform from log2 scale to normal scale and then normalise using
    % scaled robust sigmoid.
    data = 2.^(data);
    
    coordSample{sub} = coord;
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
    expSample{sub} = [ROI, data];
    coordSampROI{sub} = [ROI, coord];
    
end

expSampROI = vertcat(expSampNorm{1}, expSampNorm{2}, expSampNorm{3}, expSampNorm{4}, expSampNorm{5}, expSampNorm{6});
coordSamp = vertcat(coordSampROI{1}, coordSampROI{2}, coordSampROI{3}, coordSampROI{4}, coordSampROI{5}, coordSampROI{6});
ROIs = unique(expSampROI(:,1));
% average expression values for each ROI for differential stability calculation

expressionROI = zeros(length(ROIs),size(expSampROI,2)-1);
coordinatesROI = zeros(length(ROIs),3);

for j=1:length(ROIs)
    indROI = find(expSampROI(:,1)==(ROIs(j)));
    noProbes = length(indROI);
    fprintf(1,'%u samples for %u ROI found \n', noProbes, ROIs(j))
    % take expression values for a selected entrezID
    expressionRepInt = expSampROI(indROI,2:end);
    coordinatesRepInt = coordSamp(indROI,2:end);
    
    % calculate the mean for expression data for a selected entrezID
    expressionROI(j,:)= mean(expressionRepInt,1);
    coordinatesROI(j,:)= mean(coordinatesRepInt,1);
    
end

geneROI = [ROIs, expressionROI];
coordinatesROI = [ROIs, coordinatesROI];

% calculate DS
if strcmp(parcellation, 'aparcaseg')
    DS = calculateDS(DataExpression,DataCoordinatesMNI,parcellation, normaliseWhat, normMethod);
    DS = DS';
    
    
    %----------------------------------------------------------------------------------
    % Take top % of DS genes
    %----------------------------------------------------------------------------------
    fprintf('Selecting genes with highest differential stability \n')
    nrGenes = round(length(DS)*percentDS/100);
    
    [ b, ix ] = sort( DS(:), 'descend' );
    
    DSvalues = zeros(nrGenes, 2);
    for ii=1:nrGenes
        DSvalues(ii,2) = b(ii);
        DSvalues(ii,1) = ix(ii);
    end
    
    %----------------------------------------------------------------------------------
    % Get probeIDs for selected DS genes
    %----------------------------------------------------------------------------------
    
    probes = probeInformation.ProbeName(DSvalues(:,1));
    entrezID = probeInformation.EntrezID(DSvalues(:,1));
    geneID = probeInformation.GeneID(DSvalues(:,1));
    geneName = probeInformation.GeneName(DSvalues(:,1));
    geneSymbol = probeInformation.GeneSymbol(DSvalues(:,1));
    geneInd = DSvalues(:,1);
    geneDS = DSvalues(:,2);
    DSTable = table(geneDS, geneInd, entrezID, geneID, geneName, geneSymbol);
    probeInformation.DS = DS;
end
cd 'forAlex'
save(sprintf('%dparcellation%s.mat', NumNodes, normaliseWhat), 'geneROI',  'coordinatesROI');
if strcmp(parcellation, 'aparcaseg')
save(sprintf('DSgenes%s.mat', normaliseWhat), 'DSTable', 'probeInformation');
end
cd ../../../..

