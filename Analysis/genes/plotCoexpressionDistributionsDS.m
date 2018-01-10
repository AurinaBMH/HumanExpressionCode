% Plot distributions of DS values
percentDS = [5 20 100];
numNodes = 360; 
normMethod = 'scaledRobustSigmoid'; % '' for sigmoid; zscore for zscore;
probeSelection = {'LessNoise'};
doNormalise = 1; 

cd ('data/genes/processedData');
dataCell = cell(length(percentDS),1); 
variance = zeros(length(percentDS),1); 
for p=1:length(percentDS)
     load(sprintf('%dDS%d%s%s%d.mat', percentDS(p), numNodes, normMethod, probeSelection{1}, doNormalise)) % - generated using S5 script (probes chosen based on variance)
     dataCell{p} = averageCoexpression(:); 
     variance(p) = nanvar(averageCoexpression(:)); 
end
JitteredParallelScatter(dataCell, true, true, true);
xLabels = {'top 5 perc DS', 'top 20 perc DS', '100 perc DS'};
set(gca,'Xtick', [1 2 3], 'XTickLabel',xLabels, 'FontSize', 15);
set(gca,'Ytick', [-0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8], 'YTickLabel',[-0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8], 'FontSize', 15);
set(gca,'box','off');
ylabel('Correlated gene expression','FontSize', 15);