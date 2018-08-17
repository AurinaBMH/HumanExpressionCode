%  CGE vs distance
clear all; close all;

parc = 'HCP';
tract = 'iFOD2';
weight = 'FA';
brainPart = 'Lcortex';
groupConn = 'consistency'; %  lengthCV, consistency

[coexpData, A, matrices, coordinates, avWeight] = giveConnCoexp(parc,tract,weight,brainPart);

correctWhat = 'all'; %'connected'; 'all'; 'none'; %false;
whatGeneSet = 'MET'; % DS, HSE, DEG, SCZ
giveRC = true;

if strcmp(whatGeneSet, 'DS')
    DSthresh = 0.8;
end

for densThreshold = [0.1 0.15 0.25]
    
    
    [rgb_colorMatrix,labels] = GiveMeColors('richFeederPeripheral2');
    [Gr] = giveMeRichClub(matrices, coordinates, groupConn ,densThreshold, giveRC);
    title(sprintf('RC and distance curves, %.2f density', densThreshold))
    %figure; imagesc(log(Gr));
    Gr = logical(Gr);
    nodeDeg = degrees_und(Gr);
    G = Gr.*avWeight;
    
    if densThreshold==0.10
        switch tract
            case 'FACT'
                hubThresh = 20;
            case 'iFOD2'
                hubThresh = 25;
        end
    elseif densThreshold==0.15
        switch tract
            case 'FACT'
                hubThresh = 30;
            case 'iFOD2'
                hubThresh = 45;
        end
    elseif densThreshold==0.20
        switch tract
            case 'FACT'
                hubThresh = 50;
            case 'iFOD2'
                hubThresh = 55;
        end
    elseif densThreshold==0.25
        switch tract
            case 'FACT'
                hubThresh = 60;
            case 'iFOD2'
                hubThresh = 65;
        end
    elseif densThreshold==0.30
        switch tract
            case 'FACT'
                hubThresh = 80;
            case 'iFOD2'
                hubThresh = 70;
        end
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
    for crand=1:2
        if crand==1
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
                CGErand = corr(exprand', 'rows','complete', 'type', 'Pearson');
                CGErandall(it,:,:) = CGErand;
            end
            CGEmatrix = squeeze(nanmean(CGErandall,1));
        else
            CGEmatrix = corr(parcelExpression(:,selectGenes+1)', 'type', 'Pearson'); % calculate sample-sample coexpression
            CGEmatrix(logical(eye(size(CGEmatrix)))) = NaN; % replace diagonal with NaN
        end
        
        p = setdiff(1:length(Gr), parcelExpression(:,1));
        
        if ~isempty(p)
            CGEmatrix = insertrows(CGEmatrix,NaN,p);% insert rows
            CGEmatrix = insertrows(CGEmatrix.', NaN,p).' ; % insert columns
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
                
            case 'all'
                
                distExpVect(:,2) = CGEmatrix(:);
                Fit{1} = 'exp';
                [~,~,c] = GiveMeFit(distExpVect(:,1),distExpVect(:,2),Fit{1});
                [FitCurve] = getFitCurve(Fit,distExpVect,c);
                Residuals = distExpVect(:,2) - FitCurve;
                % reshape into matrix for further use
                CGEmatrix = reshape(Residuals,[length(Gr) length(Gr)]);
                CGEmatrix = CGEmatrix.*Gr;
            case 'none'
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
        nodeDeg = degrees_und(Grdeg);
        
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
            
            subplot(4,1,i+1); scatter(distExpVectHub(:,1),distExpVectHub(:,2),50,'MarkerEdgeColor',[0.35 0.35 0.35],'MarkerFaceColor',rgb_colorMatrix(i+1,:),'LineWidth',1.5)
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
        if crand~=1
            RichClubHuman(Grdeg,DISTcon);
            title(sprintf('Density %.2f - distance', densThreshold))
        end
        
        RichClubHuman(Grdeg,CGEmatrix);
        if doRandomSet
            title(sprintf('Density %.2f, random (%d genes) - CGE (%s corrected)', densThreshold, nrGenes, correctWhat))
        else
            title(sprintf('Density %.2f, %s (%d genes) - CGE (%s corrected)', densThreshold, whatGeneSet, nrGenes, correctWhat))
        end
    end
end










