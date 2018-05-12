% correct distance effect within cortex/subcortex and between them
% separately
clear all; close all; 

cd /Users/Aurina/GoogleDrive/Genetics_connectome/HumanExpression/data/genes/processedData
normalisedSeparately = true; 
if normalisedSeparately
    load('100DS82scaledRobustSigmoidRNAseq1LcortexSubcortexSEPARATE_ROI_NOdistCorrTEST.mat')
else
    load('100DS82scaledRobustSigmoidRNAseq1LcortexSubcortex_ROI_NOdistCorr.mat')
end

numNodes = 41;
numCort = 34;
%1. take one hal of the matrix;
avCoexp = maskuHalf(averageCoexpression);
avdist = maskuHalf(averageDistance); %.*double(logical(averageCoexpression)));
%2. Separate groups: C-C, S-S, C-S in coexpression and distance;
isCC = false(numNodes,1);
isCC(1:numCort) = 1;
% get indexes for CC, SS, CS in the matrix lower half to use later when we
% put values back to the matrix
indMAT = zeros(numNodes,numNodes);

indCC = indMAT;
indCC(isCC, isCC) = 1;
indCC = maskuHalf(indCC);
indCC = find(indCC==1);


indSS = indMAT;
indSS(~isCC, ~isCC) = 1;
indSS = maskuHalf(indSS);
indSS = find(indSS==1);


indCS = indMAT;
indCS(~isCC, isCC) = 1;
indCS = maskuHalf(indCS);
indCS = find(indCS==1);

% C-C
CCexp = avCoexp(isCC, isCC);
expSEP{1} = CCexp(~isnan(CCexp));
CCdist = avdist(isCC, isCC);
distSEP{1} = CCdist(~isnan(CCdist));
% S-S
SSexp = avCoexp(~isCC, ~isCC);
expSEP{2} = SSexp(~isnan(SSexp));
SSdist = avdist(~isCC, ~isCC);
distSEP{2} = SSdist(~isnan(SSdist));
% C-S
CSexp = avCoexp(~isCC, isCC);
expSEP{3} = CSexp(:);
CSdist = avdist(~isCC, isCC);
distSEP{3} = CSdist(:);
%3. correct separately for each subgroup
for i=1:3
    if i==2 % for SS
        numThr = 4;
    else
        numThr = 15;
    end
    [xThresholds,yMeans] = BF_PlotQuantiles(distSEP{i},expSEP{i},numThr,1,1); xlabel('Euclidean distance (mm)'); ylabel('Correlated gene expression');ylim([-1 1]); set(gca,'fontsize',15)
    Y = discretize(distSEP{i},xThresholds);
    Residuals = zeros(length(Y),1);
    for val=1:length(Y)
        if ~isnan(expSEP{i}(val))
            Residuals(val) = expSEP{i}(val) - yMeans(Y(val));
        else
            Residuals(val) = NaN;
        end
    end
    corrected{i} = Residuals;
    BF_PlotQuantiles(distSEP{i},corrected{i},numThr,1,1);   xlabel('Euclidean distance (mm)'); ylabel('Correlated gene expression');ylim([-1 1]); set(gca,'fontsize',15)
end
%4. plot it as a function of distance with different colours.
expCorr = zeros(41,41);
expCorr(indCC) = corrected{1};
expCorr(indSS) = corrected{2};
expCorr(indCS) = corrected{3};
expCorr = expCorr'+expCorr;

figure; imagesc(expCorr); caxis([-1,1]);title('Corrected coexpression');
colormap([flipud(BF_getcmap('blues',9));[1 1 1];BF_getcmap('reds',9)]);


figure; imagesc(averageCoexpression); caxis([-1,1]);title('Coexpression');
colormap([flipud(BF_getcmap('blues',9));[1 1 1];BF_getcmap('reds',9)]);
set(gcf,'color','w');

% plot distance relationships separately for groups on corrected data
figure; scatter(distSEP{1}, corrected{1}, 100, 'MarkerEdgeColor',[.45 .45 .45],...
    'MarkerFaceColor',[.64 .87 .93]);
hold on; scatter(distSEP{3}, corrected{3},100, 'MarkerEdgeColor',[.45 .45 .45],...
    'MarkerFaceColor',[1 .65 0]);
hold on; scatter(distSEP{2}, corrected{2},100, 'MarkerEdgeColor',[.45 .45 .45],...
    'MarkerFaceColor',[1 .11 .18]);
set(gcf,'color','w');set(gca,'fontsize',15)
legend({'Cortex to cortex', 'Cortex to subcortex', 'Subcortex to subcortex'})
xlabel('Euclidean distance (mm)')
ylabel('Correlated gene expression')
ylim([-1 1])


% plot distance relationships separately for groups on non-corrected data
figure; scatter(distSEP{1},expSEP{1}, 100, 'MarkerEdgeColor',[.45 .45 .45],...
    'MarkerFaceColor',[.64 .87 .93]);
hold on; scatter(distSEP{3}, expSEP{3},100, 'MarkerEdgeColor',[.45 .45 .45],...
    'MarkerFaceColor',[1 .65 0]);
hold on; scatter(distSEP{2}, expSEP{2},100, 'MarkerEdgeColor',[.45 .45 .45],...
    'MarkerFaceColor',[1 .11 .18]);
set(gcf,'color','w');set(gca,'fontsize',15)
legend({'Cortex to cortex', 'Cortex to subcortex', 'Subcortex to subcortex'})
xlabel('Euclidean distance (mm)')
ylabel('Correlated gene expression')
ylim([-1 1])
