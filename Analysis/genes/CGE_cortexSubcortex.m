% Cortex/subcortex
cd('data/genes/processedData')
load('100DS220scaledRobustSigmoidRNAseq1LcortexSubcortex_ROI_NOdistCorr.mat')

allNodes = zeros(110,1);
allNodes(1:100) = 1;
allNodes = logical(allNodes);

coexp = maskHalf(averageCoexpression);
dist = maskHalf(averageDistance);

CCmask = nan(110, 110);
CCmask(allNodes,allNodes) = 1;

CSmask = nan(110, 110);
CSmask(~allNodes, allNodes) = 1;

SSmask = nan(110, 110);
SSmask(~allNodes, ~allNodes) = 1;

CCexp = coexp.*CCmask;  CCexp = CCexp(:); CCexp(isnan(CCexp)) = [];
CSexp = coexp.*CSmask;  CSexp = CSexp(:); CSexp(isnan(CSexp)) = [];
SSexp = coexp.*SSmask;  SSexp = SSexp(:); SSexp(isnan(SSexp)) = [];


CCdist = dist.*CCmask; CCdist = CCdist(:); CCdist(isnan(CCdist)) = [];
CSdist = dist.*CSmask;  CSdist = CSdist(:); CSdist(isnan(CSdist)) = [];
SSdist = dist.*SSmask; SSdist = SSdist(:); SSdist(isnan(SSdist)) = [];


figure; scatter(CCdist, CCexp, 60, 'MarkerEdgeColor',[.45 .45 .45],...
    'MarkerFaceColor',[.64 .87 .93]);
hold on; scatter(CSdist, CSexp,60, 'MarkerEdgeColor',[.45 .45 .45],...
    'MarkerFaceColor',[1 .65 0]);
hold on; scatter(SSdist, SSexp,60, 'MarkerEdgeColor',[.45 .45 .45],...
    'MarkerFaceColor',[1 .11 .18]);
set(gcf,'color','w');
legend({'Cortex to cortex', 'Cortex to subcortex', 'Subcortex to subcortex'}, 'FontSize', 12)
xlabel('Euclidean distance (mm)', 'FontSize', 15)
ylabel('Correlated gene expression', 'FontSize', 15)
ylim([-1 1])
