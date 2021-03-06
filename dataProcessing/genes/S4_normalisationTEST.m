function S4_normalisationTEST(options)
%% Author: Aurina

%close all;
%Last modiffied: 2017-07-31
%Last modiffied: 2017-08-01
%close all;
%clear all;
%------------------------------------------------------------------------------
% Choose options
%------------------------------------------------------------------------------
useCUSTprobes = options.useCUSTprobes;  % choose if you want to use data with CUST probes
probeSelection = options.probeSelections;
parcellation = options.parcellations;
distanceThreshold = options.distanceThreshold;
correctDistance = options.correctDistance;
resolution = options.resolution;
calculateDS = options.calculateDS;
percentDS = options.percentDS;
distanceCorrection = options.distanceCorrection;
coexpressionFor = options.coexpressionFor;
Fit = options.Fit;
doNormalise = options.doNormalise;
normMethod = options.normMethod;
normaliseWhat = options.normaliseWhat;
% choose Lcortex if want tnormalise samples assigned to left cortex separately;
% choose LcortexSubcortex if want to normalise LEFT cortex + left subcortex together
% choose wholeBrain if you want to normalise the whole brain.
% choose LRcortex if you want to normalise left cortex + right cortex.
optionsSave = options;

for p=probeSelection
    
    %options.probeSelection = p;
    %------------------------------------------------------------------------------
    % Define number of subjects and parcellation details based on choises
    %------------------------------------------------------------------------------
    if strcmp(parcellation, 'aparcaseg')
        numNodes = 82;
        LeftCortex = 1:34;
        LeftSubcortex = 35:41;
        RightCortex = 42:75;
        RightSubcortex = 76:82;
    elseif strcmp(parcellation, 'cust100')
        numNodes = 220;
        LeftCortex = 1:100;
        LeftSubcortex = 101:110;
        RightCortex = 111:210;
        RightSubcortex = 211:220;
    elseif strcmp(parcellation, 'cust250')
        numNodes = 530;
        LeftCortex = 1:250;
        LeftSubcortex = 251:265;
        RightCortex = 266:515;
        RightSubcortex = 516:530;
        
    elseif strcmp(parcellation, 'HCP')
        numNodes = 360;
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
            nROIs = 1:numNodes;
        case 'LRcortex'
            subjects = 1:2;
            nROIs = [LeftCortex,RightCortex];
        case 'LcortexSubcortexSEPARATE'
            subjects = 1:6;
            nROIs = [LeftCortex,LeftSubcortex];
            
    end
    %options.subjects = subjects;
    if useCUSTprobes
        startFileName = 'MicroarrayDataWITHcustProbesUpdatedXXX';
    else
        startFileName = 'MicroarrayDataProbesUpdatedXXX';
    end
    
    cd ('data/genes/processedData');
    load(sprintf('%s%s%dDistThresh%d.mat', startFileName, p{1}, numNodes, distanceThreshold));
    
    
    expressionSubjROI = cell(6,1);
    coordinatesSubjROI = cell(6,1);
    coordSample = cell(6,1);
    expSampNorm = cell(6,1);
    expSample = cell(6,1);
    
    entrezIDs = probeInformation.EntrezID;
    %load('IDgenes2plus.mat');
    %[~, keep] = intersect(entrezIDs, IDgene);
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
                
            case 'LcortexSubcortexSEPARATE'
                ind1{1} = find(ismember(expSingleSubj(:,2), LeftCortex));
                ind1{2} = find(ismember(expSingleSubj(:,2), LeftSubcortex));
                
        end
        
        switch normaliseWhat
            case 'LcortexSubcortexSEPARATE'
                % normalise cortical and subcortical samples separately and
                % then combine them
                for kk=1:2
                    expSubj1{kk} = expSingleSubj(ind1{kk},:);
                    coord1{kk} = coordSingle(ind1{kk},3:5);
                    data = expSubj1{kk}(:,3:size(expSubj1{kk},2));
                    ROI1{kk} = expSubj1{kk}(:,2);
                    % normalise sample x gene data for each subject separately
                    %% commented - not normalise
                    switch normMethod
                        case 'hampel'
                            dataNorm1 = Norm_hampel(data);
                            fprintf('Normalising gene expression data\n')
                        otherwise
                            dataNorm1 = BF_NormalizeMatrix(data, normMethod);
                            %dataNorm2 = BF_NormalizeMatrix(dataNorm, normMethod);
                            fprintf('Normalising gene expression data\n')
                    end
                    DATA{kk} = dataNorm1;
                end
                dataNorm = vertcat(DATA{1}, DATA{2});
                ROI = vertcat(ROI1{1}, ROI1{2}); 
                coord = vertcat(coord1{1}, coord1{2}); 
                expSubj = vertcat(expSubj1{1}, expSubj1{2});
            otherwise
                expSubj = expSingleSubj(ind,:);
                coord = coordSingle(ind,3:5);
                data = expSubj(:,3:size(expSubj,2));
                %data = 2.^(data);
                %     if onlyMultipleProbes
                %         data = data(:,keep);
                %     end
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
                        %dataNorm2 = BF_NormalizeMatrix(dataNorm, normMethod);
                        fprintf('Normalising gene expression data\n')
                end
        end
        if doNormalise
            expSampNorm{sub} = [ROI, dataNorm];
        else
            expSampNorm{sub} = [ROI, data];
        end
        % SAVE NON-NORMALISED DATA FOR TEST
        coordSample{sub} = [ROI, coord];
        %expSample{sub} = [ROI, data];
        ROIs = unique(expSubj(:,2));
        
        
        % average expression values for each ROI for differential stability calculation
        
        expressionROI = zeros(length(ROIs),size(data,2));
        coordinatesROI = zeros(length(ROIs),3);
        
        for j=1:length(ROIs)
            indROI = find(expSubj(:,2)==(ROIs(j)));
            noProbes = length(indROI);
            % fprintf(1,'%u samples for %u ROI found \n', noProbes, ROIs(j))
            % take expression values for a selected entrezID
            if doNormalise
                expressionRepInt = dataNorm(indROI,:);
            else
                expressionRepInt = data(indROI,:);
            end
            % try not normalised data for DS
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
    
    if calculateDS
        
        %DS = calculateDS(DataExpression,DataCoordinatesMNI,parcellation, normaliseWhat, normMethod);
        
        %----------------------------------------------------------------------------------
        % Pre-define variables for DS calculation
        %----------------------------------------------------------------------------------
        
        numSubjects = max(subjects);
        inter = cell(numSubjects,numSubjects);
        indexj = cell(numSubjects,numSubjects);
        indexk = cell(numSubjects,numSubjects);
        indexjp = cell(numSubjects,numSubjects);
        corellations = cell(numSubjects,numSubjects);
        numGenes = size(expressionSubjROI{1,1},2)-2;
        
        %----------------------------------------------------------------------------------
        % Get ROIs that are in all subjects
        %----------------------------------------------------------------------------------
        fprintf('Selecting ROIs for DS calculation\n')
        R = cell(1,numSubjects);
        for o=1:numSubjects
            R{:,o} = expressionSubjROI{o}(:,2);
        end
        intersectROIs = mintersect(R{1}, R{2}, R{3}, R{4}, R{5}, R{6});
        
        % use a set list of ROIS that are present in all subjects
        ROIsindex = zeros(length(intersectROIs),numSubjects);
        for j=1:numSubjects
            
            for w=1:length(intersectROIs)
                ROIsindex(w,j) = find(R{j}==intersectROIs(w));
            end
            
        end
        
        %----------------------------------------------------------------------------------
        % For each pair of cubjects, calculate correlation (Spearman) of regional expression
        % for each gene to select genes that have consistent expression patterns
        % through regions between pairs of subjects.
        %----------------------------------------------------------------------------------
        fprintf('Calculating differential stability\n')
        for j=1:numSubjects
            for k=j+1:numSubjects
                
                expSubone = expressionSubjROI{j}(ROIsindex(:,j),3:numGenes+2);
                expSubtwo = expressionSubjROI{k}(ROIsindex(:,k),3:numGenes+2);
                genes = zeros(1,numGenes);
                
                for g=1:numGenes
                    genes(g) = corr(expSubone(:,g),expSubtwo(:,g),'type','Spearman');
                end
                
                corellations{j,k} = genes;
            end
        end
        
        %----------------------------------------------------------------------------------
        % Combina data for all pairs of subjects
        %----------------------------------------------------------------------------------
        fprintf('Combining data for all subjects\n')
        if numSubjects == 2
            c = vertcat(corellations{1,2});
        elseif numSubjects == 6
            c = vertcat(corellations{1,2}, corellations{1,3}, corellations{1,4}, corellations{1,5}, corellations{1,6}, ...
                corellations{2,3}, corellations{2,4}, corellations{2,5}, corellations{2,6}, corellations{3,4}, corellations{3,5}, ...
                corellations{3,6}, corellations{4,5}, corellations{4,6}, corellations{5,6});
        end
        %----------------------------------------------------------------------------------
        % Take the mean for each gene - this is DS score for a gene
        % gene that have most consistent expression pattern through regions will
        % get highest scores
        %----------------------------------------------------------------------------------
        DS = mean(c,1);
        
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
        
        %    probes = probeInformation.ProbeName(DSvalues(:,1));
        % DSProbeTable = table(probes, DSvalues(:,2));
        %----------------------------------------------------------------------------------
        % Take selected genes and calculate sample - sample coexpression
        %----------------------------------------------------------------------------------
        fprintf('Calculating coexpression between samples, performing coexpression-distance correctio and averaging coexpression to ROIs\n')
        %load('DistancesONsurface.mat');
        switch coexpressionFor
            case 'all'
                
                selectedGenes = expSampNormalisedAll(:,2:end);
                switch distanceCorrection
                    case 'Euclidean'
                        % calculate euclidean distance on MNI coordinates
                        sampleDistances = pdist2(combinedCoord(:,2:end), combinedCoord(:,2:end));%
                    case 'GMvolume'
                        % load pre-calculated distances within GM volume
                        load('distancesGM_MNIXXX.mat');
                        sampleDistances = distSamples;
                    case 'Surface'
                        % load pre-calculated distances on surface
                        load('DistancesONsurfaceXXX.mat');
                        sampleDistances = distSamples;
                end
                fprintf(sprintf('%s distance correction is chosen\n', distanceCorrection))
                W = unique(expSampNormalisedAll(:,1));
                ROIs = expSampNormalisedAll(:,1);
                
                [averageCoexpression, parcelCoexpression, correctedCoexpression, Residuals, distExpVect, averageDistance] = calculateCoexpression(sampleDistances, selectedGenes, DSvalues, W, ROIs,nROIs, Fit, correctDistance, resolution);
                
            case 'separate'
                
                expPlotALL = zeros(max(nROIs),max(nROIs),max(subjects));
                %expPlotALL2 = cell(6,1);
                correctedCoexpressionALL = cell(max(subjects),1);
                parcelCoexpressionALL = cell(max(subjects),1);
                indSub1 = 1;
                for sub=subjects
                    
                    
                    selectedGenes = expSampNorm{sub}(:,2:end);
                    indSub = size(expSampNorm{sub},1);
                    sampleDistances = pdist2(coordSample{sub}(:,2:end), coordSample{sub}(:,2:end)); %distancesMNI(indSub1:indSub, indSub1:indSub); %
                    indSub1 = indSub+1;
                    W = unique(expSampNorm{sub}(:,1));
                    ROIs = expSampNorm{sub}(:,1);
                    
                    [expPlot, parcelCoexpression, correctedCoexpression, Residuals, distExpVect, averageDistance] = calculateCoexpression(sampleDistances, selectedGenes, DSvalues, W, ROIs,nROIs, Fit, correctDistance, resolution);
                    expPlotALL(:,:,sub) = expPlot;
                    %expPlotALL2{sub} = expPlot;
                    correctedCoexpressionALL{sub} = correctedCoexpression;
                    parcelCoexpressionALL{sub} = parcelCoexpression;
                    
                end
                averageCoexpression = nanmean(expPlotALL,3);
        end
        
        
        probeInformation.DS = DS';
    end
    SampleCoordinates = sortrows(combinedCoord,1);
    SampleGeneExpression = sortrows(expSampNormalisedAll,1);
    %save(sprintf('DSnew%s', p{1}), 'DS', 'averageCoexpression', 'DSProbeTable', 'expSampNormalisedAll', 'probeInformation');
    if correctDistance
        save(sprintf('%dDS%d%s%s%d%s_%s_distCorrTEST.mat', percentDS, numNodes, normMethod, p{1}, doNormalise, normaliseWhat, resolution), 'SampleCoordinates', 'SampleGeneExpression', 'probeInformation', 'optionsSave', 'averageCoexpression', 'averageDistance');
        fprintf('saving to file1')
    elseif ~correctDistance
        save(sprintf('%dDS%d%s%s%d%s_%s_NOdistCorrTEST.mat', percentDS, numNodes, normMethod, p{1}, doNormalise, normaliseWhat, resolution), 'SampleCoordinates', 'SampleGeneExpression', 'probeInformation', 'optionsSave', 'averageCoexpression', 'averageDistance');
        fprintf('saving to file2')
    end
    %save(sprintf('DSnew%s%s%dBEN', normMethod, p{1}, doNormalise), 'SampleCoordinates', 'SampleGeneExpression', 'probeInformation', 'options');
    cd ../../..
end
%figure; imagesc(expPlotMNI); caxis([-1 1]); colormap([flipud(BF_getcmap('blues',9));[1 1 1];BF_getcmap('reds',9)]); title('Average coexpression')

% A = [averageCoexpressionSeparateMasMin(:),averageCoexpressionSeparateZscore(:)];
% A = A(~any(isnan(A),2),:);
% figure; scatter(A(:,1), A(:,2));
% [r,p] = corr(A(:,1), A(:,2), 'type', 'Spearman')


end