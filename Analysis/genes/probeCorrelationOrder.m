clear all;

cd ('data/genes/processedData');
numProbes = 2;
normMethod = 'scaledRobustSigmoid'; % '' for sigmoid; zscore for zscore;
coexpressionOn = 'sample'; % sample of gene
probeSelections = {'Variance', 'LessNoise', 'Mean', 'PC', 'Random'}; %, 'Random'};
%DSthreshold = -1; % -1 will include all.
doNormalise = true; 
numNodes = 360; 

DSscoresAll = cell(5,1);
allData = cell(5,1);

for i=1:length(probeSelections)
    load(sprintf('DS%d%s%s%d', numNodes, normMethod, probeSelections{i}, doNormalise))
    fprintf('Loading data with probes selected based on %s\n',probeSelections{i})
    %fprintf('Reordering probes\n')% - generated using S5 script (probes chosen based on variance)
    DSscoresAll{i} = probeInformation.DS; probeSelections{i} = probeSelections{i};
    allData{i} = SampleGeneExpression(:,2:end);
end

% entrezIDs for genes that have more than X probes.
fprintf('Loading gene entrezIDs that have %d or more probes\n', numProbes)
load(sprintf('IDgenes%dplus.mat', numProbes));

%DSscores = DS>DSthreshold;
entrezIDs = probeInformation.EntrezID; %probeInformation.EntrezID;
% select those genes that have more than X probes (and maybe also high DS)
[~, keep] = intersect(entrezIDs, IDgene);

allDatafinal = cell(length(probeSelections),1);
coexpValues = cell(length(probeSelections),1);
for k=1:length(probeSelections)
    
    % try to normaise data in columns before calculating coexpression
    selectedData = allData{k}(:,keep);
    allDatafinal{k} = selectedData; %BF_NormalizeMatrix(selectedData, 'zscore');
    switch coexpressionOn
        case 'sample'
            mat = corr(allDatafinal{k}', 'type', 'Spearman');
        case 'gene'
            mat = corr(allDatafinal{k}, 'type', 'Spearman');
    end

    coexpValues{k} = mat(:);
    
end

numGenes = size(allDatafinal{1},2);
avCorr = zeros(length(probeSelections),length(probeSelections));
coexpCorr = zeros(length(probeSelections),length(probeSelections));
% calculate correlation between each way of choosing a probe
for i=1:length(probeSelections)
    for j=i+1:length(probeSelections)
        
        correlation = zeros(numGenes,1);
        for g=1:numGenes
            correlation(g) = corr(allDatafinal{i}(:,g), allDatafinal{j}(:,g), 'type', 'Spearman');
        end
        
        avCorr(j,i) = mean(correlation);
        coexpCorr(j,i) = corr(coexpValues{i}, coexpValues{j}, 'type', 'Spearman');
        
    end
end

figure; imagesc(avCorr); colormap([[1 1 1];BF_getcmap('reds',9)]); caxis([0 1]); colorbar;
xticks(1:length(probeSelections)); yticks(1:length(probeSelections));
title(sprintf('Correlation between genes with %d and more probes (%d genes)', numProbes, numGenes));
xticklabels(probeSelections);
yticklabels(probeSelections);

figure; imagesc(coexpCorr); colormap([[1 1 1];BF_getcmap('reds',9)]); colorbar;
xticks(1:length(probeSelections)); yticks(1:length(probeSelections));
title(sprintf('Correlation between %s-%s coexpression matrices for genes with %d and more probes (%d genes)', coexpressionOn,coexpressionOn,numProbes, numGenes));
xticklabels(probeSelections);
yticklabels(probeSelections);
