clear all;

cd ('data/genes/processedData');
numProbes = 2;
normMethod = ''; % '' for sigmoid; zscore for zscore;
coexpressionOn = 'sample'; % sample of gene
probeSelection = {'Variance', 'PC', 'LessNoise', 'Mean', 'Random'};
DSthreshold = -1; % -1 will include all.

DSscoresAll = cell(5,1);
allData = cell(5,1);

for i=1:length(probeSelection)
    load(sprintf('DSnew%s%s.mat', normMethod, probeSelection{i})) % - generated using S5 script (probes chosen based on variance)
    [a,order] = sort(probeInformation.EntrezID);
    DSscoresAll{i} = DS; probeSelection{i} = probeSelection{i};
    data = expSampNormalisedAll(:,2:end);
    allData{i} = data(:, order);
end

% entrezIDs for genes that have more than X probes.
load(sprintf('IDgenes%dplus.mat', numProbes));

DSscores = DS>DSthreshold;
entrezIDs = a; %probeInformation.EntrezID;
% select those genes that have more than X probes (and maybe also high DS)
%[~, keep] = intersect(entrezIDs, IDgene);
[~, keep] = intersect(entrezIDs(DSscores==1), IDgene);

allDatafinal = cell(5,1);
coexpValues = cell(5,1);
for k=1:5
    
    % try to normaise data in columns before calculating coexpression
    selectedData = allData{k}(:,keep);
    allDatafinal{k} = selectedData; %BF_NormalizeMatrix(selectedData, 'zscore');
    switch coexpressionOn
        case 'sample'
            %mat = corr(allDatafinal{k}, 'type', 'Spearman');
            mat = corr(allDatafinal{k}', 'type', 'Spearman');
        case 'gene'
            mat = corr(allDatafinal{k}, 'type', 'Spearman');
    end
    %coexpression.
    coexpValues{k} = mat(:);
    
end

numGenes = size(allDatafinal{1},2);
avCorr = zeros(5,5);
coexpCorr = zeros(5,5);
% calculate correlation between each way of choosing a probe
for i=1:5
    for j=i+1:5
        
        correlation = zeros(numGenes,1);
        for g=1:numGenes
            correlation(g) = corr(allDatafinal{i}(:,g), allDatafinal{j}(:,g), 'type', 'Spearman');
        end
        
        avCorr(j,i) = mean(correlation);
        coexpCorr(j,i) = corr(coexpValues{i}, coexpValues{j}, 'type', 'Spearman');
        
    end
end

figure; imagesc(avCorr); colormap([[1 1 1];BF_getcmap('reds',9)]); caxis([0 1]); colorbar;
xticks([1 2 3 4 5]); yticks([1 2 3 4 5]);
title(sprintf('Correlation between genes with %d and more probes (%d genes)', numProbes, numGenes));
xticklabels({'Variance','PC','noise','mean', 'random'});
yticklabels({'Variance','PC','noise','mean', 'random'});

figure; imagesc(coexpCorr); colormap([[1 1 1];BF_getcmap('reds',9)]); colorbar;
xticks([1 2 3 4 5]); yticks([1 2 3 4 5]);
title(sprintf('Correlation between %s-%s coexpression matrices for genes with %d and more probes (%d genes)', coexpressionOn,coexpressionOn,numProbes, numGenes));
xticklabels({'Variance','PC','noise','mean', 'random'});
yticklabels({'Variance','PC','noise','mean', 'random'});
