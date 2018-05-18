% Cortex/subcortex
cd('data/genes/processedData')
load('100DS82scaledRobustSigmoidRNAseq1LcortexSubcortex_ROI_NOdistCorr.mat')
numCort = 34;
numNodes = 41;

allNodes = zeros(numNodes,1);
allNodes(1:numCort) = 1;
allNodes = logical(allNodes);

coexp = maskHalf(averageCoexpression);
dist = maskHalf(averageDistance);

CCmask = nan(numNodes, numNodes);
CCmask(allNodes,allNodes) = 1;

CSmask = nan(numNodes, numNodes);
CSmask(~allNodes, allNodes) = 1;

SSmask = nan(numNodes, numNodes);
SSmask(~allNodes, ~allNodes) = 1;

CCexp = coexp.*CCmask;  CCexp = CCexp(:); CCexp(isnan(CCexp)) = [];
CSexp = coexp.*CSmask;  CSexp = CSexp(:); CSexp(isnan(CSexp)) = [];
SSexp = coexp.*SSmask;  SSexp = SSexp(:); SSexp(isnan(SSexp)) = [];


CCdist = dist.*CCmask; CCdist = CCdist(:); CCdist(isnan(CCdist)) = [];
CSdist = dist.*CSmask;  CSdist = CSdist(:); CSdist(isnan(CSdist)) = [];
SSdist = dist.*SSmask; SSdist = SSdist(:); SSdist(isnan(SSdist)) = [];


figure; scatter(CCdist, CCexp, 100, 'MarkerEdgeColor',[.55 .55 .55],...
    'MarkerFaceColor',[.64 .87 .93]);
hold on; scatter(CSdist, CSexp,100, 'MarkerEdgeColor',[.55 .55 .55],...
    'MarkerFaceColor',[1 .65 0]);
hold on; scatter(SSdist, SSexp,100, 'MarkerEdgeColor',[.55 .55 .55],...
    'MarkerFaceColor',[1 .11 .18]);
set(gcf,'color','w');set(gca,'fontsize',18)
%legend({'Cortex to cortex', 'Cortex to subcortex', 'Subcortex to subcortex'}, 'FontSize', 18); 
xlabel('Euclidean distance (mm)')
ylabel('Correlated gene expression')
ylim([-1 1])
xticks([20 40 60 80 100 120 140 160])
xlim([0 160])


% Cortex/subcortex
load('100DS82scaledRobustSigmoidRNAseq1LcortexSubcortex_ROI_distCorr.mat')
numCort = 34;
numNodes = 41;

allNodes = zeros(numNodes,1);
allNodes(1:numCort) = 1;
allNodes = logical(allNodes);

coexp = maskHalf(averageCoexpression);
dist = maskHalf(averageDistance);

CCmask = nan(numNodes, numNodes);
CCmask(allNodes,allNodes) = 1;

CSmask = nan(numNodes, numNodes);
CSmask(~allNodes, allNodes) = 1;

SSmask = nan(numNodes, numNodes);
SSmask(~allNodes, ~allNodes) = 1;

CCexp = coexp.*CCmask;  CCexp = CCexp(:); CCexp(isnan(CCexp)) = [];
CSexp = coexp.*CSmask;  CSexp = CSexp(:); CSexp(isnan(CSexp)) = [];
SSexp = coexp.*SSmask;  SSexp = SSexp(:); SSexp(isnan(SSexp)) = [];


CCdist = dist.*CCmask; CCdist = CCdist(:); CCdist(isnan(CCdist)) = [];
CSdist = dist.*CSmask;  CSdist = CSdist(:); CSdist(isnan(CSdist)) = [];
SSdist = dist.*SSmask; SSdist = SSdist(:); SSdist(isnan(SSdist)) = [];


figure; scatter(CCdist, CCexp, 100, 'MarkerEdgeColor',[.55 .55 .55],...
    'MarkerFaceColor',[.64 .87 .93]);
hold on; scatter(CSdist, CSexp,100, 'MarkerEdgeColor',[.55 .55 .55],...
    'MarkerFaceColor',[1 .65 0]);
hold on; scatter(SSdist, SSexp,100, 'MarkerEdgeColor',[.55 .55 .55],...
    'MarkerFaceColor',[1 .11 .18]);
set(gcf,'color','w');set(gca,'fontsize',18)
%legend({'Cortex to cortex', 'Cortex to subcortex', 'Subcortex to subcortex'},'FontSize', 18);
xlabel('Euclidean distance (mm)')
ylabel('Correlated gene expression');
ylim([-1 1])
xticks([20 40 60 80 100 120 140 160])
xlim([0 160])



