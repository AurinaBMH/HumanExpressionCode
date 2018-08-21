%  CGE vs distance
clear all; close all;

parc = 'HCP';
tract = 'iFOD2';
weight = 'standard';
brainPart = 'wholeBrain';
groupConn = 'variance'; %  lengthCV, consistency

[coexpData, A, matrices, coordinates, avWeight] = giveConnCoexp(parc,tract,weight,brainPart);
% if strcmp(weight,'standard')
%     % for standard weight (where streamline count is used), divide wegith
%     % by the length of the streamline
%     [~, Alength, matricesLength] = giveConnCoexp(parc,tract,'length',brainPart);
%     for s=1:length(matrices)
%         matrices{s} = matrices{s}./matricesLength{s};
%         matrices{s}(isnan(matrices{s})) = 0; 
%         A(:,:,s) = matrices{s};
%     end
%     
% end
    
% [coexpData, A, matrices, coordinates, avWeight] = giveConnCoexp(parc,tract,weight,'Lcortex');

correctWhat = 'connected'; %'connected'; 'all'; 'none'; %false;
whatGeneSet = 'ALL'; % DS, HSE, DEG, SCZ
giveRC = false;

if strcmp(whatGeneSet, 'DS')
    DSthresh = 0.8;
end

[rgb_colorMatrix,labels] = GiveMeColors('RFPU');
for densThreshold = 0.1

    [Gr] = giveMeRichClub(matrices, coordinates, groupConn ,densThreshold, giveRC);
    %title(sprintf('RC and distance curves, %.2f density', densThreshold))
    Gr = logical(Gr);
    nodeDeg = degrees_und(Gr);
    hubThresh = quantile(nodeDeg,0.5) + quantile(nodeDeg,0.25);
    
    if strcmp(brainPart, 'wholeBrain')
             nodeDeg = nodeDeg(1:180);
             Gr = Gr(1:180,1:180); 
    end

    % see CGE-distance relationship for all connected regions
    parcelExpression = coexpData.parcelExpression;
    probeInformation = coexpData.probeInformation;
    switch whatGeneSet
        case 'ALL'
            selectGenes = find(probeInformation.DS);
            nrGenes = size(coexpData.parcelExpression,2)-1;
        case 'MET'
            geneNames = importMouseGenes('mouse_metabolicGenes.txt');
            mgInd = cell(length(geneNames),1);
            for mg=1:length(geneNames)
                mgInd{mg} = find(strcmp(geneNames{mg},probeInformation.GeneSymbol));
            end
            selectGenes = cell2mat(mgInd);
            nrGenes = length(selectGenes);
        case 'DS'
            selectGenes = find(probeInformation.DS>DSthresh);
            nrGenes = length(selectGenes);
        case 'HSE'
            geneNames = {'BEND5', 'C1QL2', 'CACNA1E', 'COL24A1', 'COL6A1', 'CRYM', 'KCNC3', 'KCNH4', 'LGALS1', ...
                'MFGE8', 'NEFH', 'PRSS12', 'SCN3B', 'SCN4B', 'SNCG', 'SV2C', 'SYT2', 'TPBG', 'VAMP1'};
            mgInd = cell(length(geneNames),1);
            for mg=1:length(geneNames)
                mgInd{mg} = find(strcmp(geneNames{mg},probeInformation.GeneSymbol));
            end
            selectGenes = cell2mat(mgInd);
            nrGenes = length(selectGenes);
        case 'DEG'
            corrDeg = zeros(size(parcelExpression,2)-1,2);
            for g=1:size(parcelExpression,2)-1
                [corrDeg(g,1),corrDeg(g,2)] = corr(parcelExpression(:,g+1), nodeDeg(parcelExpression(:,1))', 'rows','complete', 'type', 'Pearson');
            end
            % find those genes that show significant correlations and make CGE matrix
            % with them
            geneINDdegsign = find(corrDeg(:,2)<0.05); geneINDdegpos = find(corrDeg(:,1)>0.3);
            selectGenes = intersect(geneINDdegsign, geneINDdegpos);
            nrGenes = length(selectGenes);
            % make a list of genes for enrichment
            geneListDeg = probeInformation.EntrezID;
            geneListDeg(selectGenes,2) = 1; %abs(corrDeg(geneINDdeg,1));
        case 'SCZ'
            % make a list of SCZ genes
    end
    for doRand=1
        if doRand==1
            doRandomSet = false;
        else
            doRandomSet = true;
        end
        
        if doRandomSet
            % compare to a random set of genes
            CGErandall = zeros(100, length(parcelExpression(:,1)), length(parcelExpression(:,1)));
            for it=1:100
                ran = randi([2, size(parcelExpression,2)-1],1,nrGenes);
                exprand = parcelExpression(:,ran);
                % calculate coexpression of this set of genes
                
                CGErand = corr(exprand', 'type', 'Spearman');
                CGErand(logical(eye(size(CGErand)))) = NaN; % replace diagonal with NaN
                CGErandall(it,:,:) = CGErand;
            end
            CGEmatrix = squeeze(nanmean(CGErandall,1));
        else
            CGEmatrix = corr(parcelExpression(:,selectGenes+1)', 'type', 'Spearman'); % calculate sample-sample coexpression
            CGEmatrix(logical(eye(size(CGEmatrix)))) = NaN; % replace diagonal with NaN
        end
        
        DISTcon = coexpData.averageDistance;
        DISTcon(logical(eye(size(DISTcon)))) = NaN;
        distExpVect = zeros(length(CGEmatrix(:)),2);
        distExpVect(:,1) = DISTcon(:);
        
        % try to correct distance effect for connected/all or none
        Gr = double(Gr);
        Gr(Gr==0) = NaN;
        
        switch correctWhat
            case 'connected'
                
                % make mask for unconnected links, so we can map them back
                % after correction for only connected is done
                CGEmatrixUncon = zeros(size(CGEmatrix)); 
                CGEmatrixUncon(isnan(Gr)) = 1;
                CGEmatrixUncon(~isnan(Gr)) = 0; 
                
                CGEmatrix = CGEmatrix.*Gr;
                %CGEmatrix(CGEmatrix==0) = NaN;
                distExpVect(:,2) = CGEmatrix(:);
                
                %indNAN = find(isnan(distExpVect(:,2)));
                Fit{1} = 'linear';
                [~,~,c] = GiveMeFit(distExpVect(:,1),distExpVect(:,2),Fit{1});
                [FitCurve] = getFitCurve(Fit,distExpVect,c);
                Residuals = distExpVect(:,2) - FitCurve;
                % reshape into matrix for further use
                CGEmatrix = reshape(Residuals,[length(Gr) length(Gr)]);
                
                CGE1 = CGEmatrixUncon.*coexpData.averageCoexpression; 
                CGE2 = CGEmatrix; 
                CGE2(isnan(CGE2)) = 0; 
                CGEmatrixFull = CGE2+CGE1; 
                
            case 'all'
                
                distExpVect(:,2) = CGEmatrix(:);
                Fit{1} = 'exp';
                [~,~,c] = GiveMeFit(distExpVect(:,1),distExpVect(:,2),Fit{1});
                [FitCurve] = getFitCurve(Fit,distExpVect,c);
                Residuals = distExpVect(:,2) - FitCurve;
                % reshape into matrix for further use
                CGEmatrix = reshape(Residuals,[length(Gr) length(Gr)]);
                CGEmatrixFull = CGEmatrix;
                CGEmatrix = CGEmatrix.*Gr;
                
            case 'none'
                CGEmatrixFull = CGEmatrix;
                CGEmatrix = CGEmatrix.*Gr;
        end
        distExpVect(:,2) = CGEmatrix(:);
        distExpVect(any(isnan(distExpVect), 2),:)=[];
        figure; set(gcf,'Position',[300 300 500 1300])
        set(gcf,'color','w');
        subplot(4,1,1); scatter(distExpVect(:,1),distExpVect(:,2),50,'MarkerEdgeColor',[0.35 0.35 0.35],'MarkerFaceColor',[0.65 0.65 0.65],'LineWidth',1.5)
        [r,p] = corr(distExpVect(:,1),distExpVect(:,2));
        hline = refline;
        ylabel('CGE')
        xlabel('Euclidean distance between regions')
        if doRandomSet
            title(sprintf('Between all connected regions (%s corrected), %.2f density, random (%d genes)', correctWhat, densThreshold, nrGenes))
        else
            title(sprintf('Between all connected regions (%s corrected), %.2f density, %s (%d genes)', correctWhat, densThreshold, whatGeneSet, nrGenes))
        end
        xlim([0 160])
        ylim([-1 1])
        text(120, 0.7, sprintf('corr = %.2f',r))
        text(120, 0.5, sprintf('p = %d', p))
        
        % see CGE-distance relationship for hubs
         Grdeg = Gr;
         Grdeg(isnan(Grdeg)) = 0;
         numNodes = size(Gr,1);
%         nodeDeg = degrees_und(Grdeg);
        
        isHub = nodeDeg>hubThresh;
        maskHub = zeros(numNodes, numNodes);
        maskHub(isHub, isHub) = 1; % rich
        maskHub(~isHub, isHub) = 2; % feeder
        maskHub(isHub, ~isHub) = 2; % feeder
        maskHub(~isHub, ~isHub) = 3; % peripheral
        maskHub(maskHub==0) =  NaN;
        
        for i=1:3
            DISTconHub = DISTcon.*(maskHub==i);
            CGEconHub = CGEmatrix.*(maskHub==i);
            DISTconHub(DISTconHub==0) = NaN;
            CGEconHub(CGEconHub==0) = NaN;
            
            distExpVectHub = nan(length(CGEmatrix(:)),2);
            distExpVectHub(:,1) = DISTconHub(:);
            distExpVectHub(:,2) = CGEconHub(:);
            distExpVectHub(any(isnan(distExpVectHub), 2),:)=[];
            
            subplot(4,1,i+1); scatter(distExpVectHub(:,1),distExpVectHub(:,2),50,'MarkerEdgeColor',[0.35 0.35 0.35],'MarkerFaceColor',rgb_colorMatrix(i,:),'LineWidth',1.5)
            [r,p] = corr(distExpVectHub(:,1),distExpVectHub(:,2));
            hline = refline;
            ylabel('CGE')
            xlabel('Euclidean distance between regions')
            if i==1
                title(sprintf('Between connected hubs (k>%d), %.2f density', hubThresh, densThreshold))
            elseif i==2
                title(sprintf('Between connected feeders, %.2f density', densThreshold))
            elseif i==3
                title(sprintf('Between connected peripheral, %.2f density', densThreshold))
            end
            xlim([0 160])
            ylim([-1 1])
            text(120, 0.7, sprintf('corr = %.2f',r))
            text(120, 0.5, sprintf('p = %d', p))
        end
        
        % % using degree as x value
        if doRand~=1
            RichClubHuman(Grdeg,DISTcon, nodeDeg);
            title(sprintf('Density %.2f - distance', densThreshold))
        end
        
        RichClubHuman(Grdeg,CGEmatrix,nodeDeg)
        if doRandomSet
            title(sprintf('Density %.2f, random (%d genes) - CGE (%s corrected)', densThreshold, nrGenes, correctWhat))
        else
            title(sprintf('Density %.2f, %s (%d genes) - CGE (%s corrected)', densThreshold, whatGeneSet, nrGenes, correctWhat))
        end
        
        % calculate distributions for rich/feeder/peripheral/unconnected
        % links within each distance bin
        
        % define distance bins based on connected edges (unconnected will 
%         A = CGEmatrix; 
%         A(maskHub~=1) = NaN; %.*(maskHub==1); 
        [xThresholds,yMeans, yMedians] = BF_PlotQuantiles(DISTcon(:),CGEmatrix(:),11,0,1); 
%         [xThresholds,yMeans, yMedians] = BF_PlotQuantiles(DISTcon(:),A(:),11,0,1); 
        [Y] = discretize(DISTcon,xThresholds); 
        
        f = figure('color', 'w');
        f.Position = [1000,200,2000,1000];
        for b=1:max(Y(:))
            
            maskBin = Y==b; % all links in this distance bin
            AdjCon = Gr.*Y==b; % connected links in this distance bin
            AdjUnc = maskBin & ~AdjCon; % unconnected links in this distance bin

            rgb_colorMatrix = GiveMeColors('RFPU');
            
            subplot(2,5,b);
            [dataCell, p, stats] = AA_CGE_RFPU(CGEmatrixFull, AdjCon, AdjUnc, nodeDeg, hubThresh, rgb_colorMatrix);
            conUncon(b,2) = p; 
            conUncon(b,1) = stats.tstat; 
            title(sprintf('Connections %.2f - %.2f mm', xThresholds(b), xThresholds(b+1)))
            hold on; 

        end

    end
end

% figure; subplot(2,1,1); histogram(degrees_und(FACT), 50); xlim([0, max([degrees_und(iFOD), degrees_und(FACT)])]); title('FACT'); subplot(2,1,2); ...
%     histogram(degrees_und(iFOD), 50); title('iFOD'); xlim([0, max([degrees_und(iFOD), degrees_und(FACT)])]); 







