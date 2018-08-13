%  CGE vs distance

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
        coordinates{m} = COG{m}(keepNodes,:);
        A(:,:,m) = matrices{m};
    else
        matrices{m} = ADJS{m};
        A(:,:,m) = matrices{m};
    end
end

nanA = A;
nanA(nanA==0) = NaN;
avWeight = nanmean(nanA,3);
avWeight(isnan(avWeight)) = 0;

distCorr = 'NO';
DS = 100;
correctConnected = false;

for densThreshold = [0.1]
    
    [rgb_colorMatrix,labels] = GiveMeColors('richFeederPeripheral2');
    if densThreshold==0.10
        hubThresh = 20;
    elseif densThreshold==0.15
        hubThresh = 30;
    elseif densThreshold==0.20
        hubThresh = 50;
    elseif densThreshold==0.25
        hubThresh = 60;
    elseif densThreshold==0.30
        hubThresh = 80;
    end
    
    
    Gr = giveMeGroupAdj_variance(matrices, densThreshold);
    Gr = logical(Gr);
    G = Gr.*avWeight;
    
    load(sprintf('%dDS%dscaledRobustSigmoidNGRNAseqQC1Lcortex_ROI_%sdistCorrEuclidean.mat', DS, nROIs, distCorr))
    
    % see CGE-distance relationship for all connected regions
    CGEcon = averageCoexpression.*Gr;
    DISTcon = averageDistance;
    CGEcon(CGEcon==0) = NaN; DISTcon(DISTcon==0) = NaN;
    distExpVect = zeros(length(CGEcon(:)),2);
    distExpVect(:,1) = DISTcon(:);
    distExpVect(:,2) = CGEcon(:);
    % try to correct distance effect for connected links only
    if correctConnected
        indNAN = find(isnan(distExpVect(:,2)));
        Fit{1} = 'linear';
        [~,~,c] = GiveMeFit(distExpVect(:,1),distExpVect(:,2),Fit{1});
        switch Fit{1}
            case 'linear'
                FitCurve = c.p1*distExpVect(:,1) + c.p2;
        end
        % get residuals
        Residuals = distExpVect(:,2) - FitCurve;
        %figure; scatter(distExpVect(:,1), Residuals);
        distExpVect(:,2) = Residuals;
        % reshape into matrix for further use
        CGEcon = reshape(Residuals,[length(G) length(G)]);
    end
    
    distExpVect(any(isnan(distExpVect), 2),:)=[];
    figure; set(gcf,'Position',[300 300 500 1300])
    set(gcf,'color','w');
    subplot(4,1,1); scatter(distExpVect(:,1),distExpVect(:,2),50,'MarkerEdgeColor',[0.35 0.35 0.35],'MarkerFaceColor',[0.65 0.65 0.65],'LineWidth',1.5)
    [r,p] = corr(distExpVect(:,1),distExpVect(:,2));
    ylabel('CGE')
    xlabel('Euclidean distance between regions')
    title(sprintf('Between all connected regions, %d density', densThreshold))
    xlim([0 160])
    ylim([-1 1])
    text(120, 0.7, sprintf('corr = %d',r))
    text(120, 0.5, sprintf('p = %d', p))
    
    % see CGE-distance relationship for hubs
    nodeDeg = degrees_und(Gr);
    nodeStrength = strengths_und(G);
    isHub = nodeDeg>hubThresh;
    maskHub = zeros(numNodes, numNodes);
    maskHub(isHub, isHub) = 1; % rich
    maskHub(~isHub, isHub) = 2; % feeder
    maskHub(isHub, ~isHub) = 2; % feeder
    maskHub(~isHub, ~isHub) = 3; % peripheral
    maskHub(maskHub==0) =  NaN;
    %indHub = maskHub(:);
    for i=1:3
        DISTconHub = DISTcon.*(maskHub==i);
        CGEconHub = CGEcon.*(maskHub==i);
        DISTconHub(DISTconHub==0) = NaN;
        CGEconHub(CGEconHub==0) = NaN;
        
        distExpVectHub = nan(length(CGEcon(:)),2);
        distExpVectHub(:,1) = DISTconHub(:);
        distExpVectHub(:,2) = CGEconHub(:);
        distExpVectHub(any(isnan(distExpVectHub), 2),:)=[];
        
        subplot(4,1,i+1); scatter(distExpVectHub(:,1),distExpVectHub(:,2),50,'MarkerEdgeColor',[0.35 0.35 0.35],'MarkerFaceColor',rgb_colorMatrix(i+1,:),'LineWidth',1.5)
        [r,p] = corr(distExpVectHub(:,1),distExpVectHub(:,2));
        hline = refline;
        ylabel('CGE')
        xlabel('Euclidean distance between regions')
        if i==1
            title(sprintf('Between connected hubs, %d density', densThreshold))
        elseif i==2
            title(sprintf('Between connected feeders, %d density', densThreshold))
        elseif i==3
            title(sprintf('Between connected peripheral, %d density', densThreshold))
        end
        xlim([0 160])
        ylim([-1 1])
        text(120, 0.7, sprintf('corr = %d',r))
        text(120, 0.5, sprintf('p = %d', p))
    end
    
    % % using degree as x value
    RichClubHuman(Gr,DISTcon);
    title(sprintf('Density %d, %d DS, %s distance correction - distance', densThreshold, DS, distCorr))
    ylim([50 65])
    
    RichClubHuman(Gr,averageCoexpression);
    title(sprintf('Density %d, %d DS, %s distance correction - CGEe', densThreshold, DS, distCorr))
    % % using strength as x value
    % RichClubHuman(G,averageCoexpression);
    % title(sprintf('Density %d, %d DS, %s distance correction', densThreshold, DS, distCorr))
    
    
    
    % check if there are genes that are correlated with degree
    
    corrDeg = zeros(size(parcelExpression,2)-1,2);
    corrStrength = zeros(size(parcelExpression,2)-1,2);
    
    for g=1:size(parcelExpression,2)-1
        [corrDeg(g,1),corrDeg(g,2)] = corr(parcelExpression(:,g+1), nodeDeg(parcelExpression(:,1))', 'rows','complete', 'type', 'Spearman');
        [corrStrength(g,1),corrStrength(g,2)] = corr(parcelExpression(:,g+1), nodeStrength(parcelExpression(:,1))', 'rows','complete','type', 'Spearman');
    end
    % find those genes that show significant correlations and make CGE matrix
    % with them
    geneINDdeg = find(corrDeg(:,2)<0.01);
    geneINDstr = find(corrStrength(:,2)<0.000001);
    parcelCoexpressionDeg = corr(parcelExpression(:,geneINDdeg+1)', 'type', 'Spearman'); % calculate sample-sample coexpression
    parcelCoexpressionDeg(logical(eye(size(parcelCoexpressionDeg)))) = NaN; % replace diagonal with NaN
    
    parcelCoexpressionStr = corr(parcelExpression(:,geneINDstr+1)', 'type', 'Spearman'); % calculate sample-sample coexpression
    parcelCoexpressionStr(logical(eye(size(parcelCoexpressionStr)))) = NaN; % replace diagonal with NaN
    
    p = setdiff(1:length(G), parcelExpression(:,1));
    
    if ~isempty(p)
        parcelCoexpressionDeg = insertrows(parcelCoexpressionDeg,NaN,p);
        parcelCoexpressionDeg = insertrows(parcelCoexpressionDeg.', NaN,p).' ; % insert columns
        
        parcelCoexpressionStr = insertrows(parcelCoexpressionStr,NaN,p);
        parcelCoexpressionStr = insertrows(parcelCoexpressionStr.', NaN,p).' ; % insert columns
        
    end
    
    figure; imagesc(parcelCoexpressionDeg); caxis([-1 1]); title(sprintf('Significant genes (deg) %d density, %d genes', densThreshold, length(geneINDdeg)))
    figure; imagesc(parcelCoexpressionStr); caxis([-1 1]); title(sprintf('Significant genes (str) %d density, %d genes', densThreshold, length(geneINDstr)))
    
    RichClubHuman(Gr,parcelCoexpressionDeg);
    RichClubHuman(Gr,parcelCoexpressionStr);
    
    % make a list of genes for enrichment
    geneListDeg = probeInformation.EntrezID; 
    geneListDeg(geneINDdeg,2) = 1; %abs(corrDeg(geneINDdeg,1)); 
    
    geneListStr = probeInformation.EntrezID; 
    geneListStr(geneINDstr,2) = 1; %abs(corrDeg(geneINDstr,1)); 
    
    
    
end













