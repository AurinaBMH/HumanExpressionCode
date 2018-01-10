clear all;
doEqual = false;
useCUSTprobes = true;
coexpressionOn = 'gene';
signalThreshold = 0.5; % percentage of samples that a selected probe has expression levels that are higher than background
doOriginal = false; %false;
probeSelection = {'Variance','PC','LessNoise', 'Mean', 'random'};
%numProbes = 3;
%------------------------------------------------------------------------------
% Load the data
%------------------------------------------------------------------------------
cd ('data/genes/processedData');
if useCUSTprobes
    fprintf(1,'Loading the data with CUST probes and assigning variables\n')
    startFileName = 'MicroarrayDataWITHcust';
else
    fprintf(1,'Loading the data without CUST probes and assigning variables\n')
    startFileName = 'MicroarrayData';
    
end
load(sprintf('%s.mat', startFileName));

ProbeID = DataTableProbe.ProbeID{1,1};
% % ------------------------------------------------------------------------------
% % First, find probes that have very noisy data and remove them from consideration
% % Threshold for removing those proges is defined as the percentage of
% % samples a probe has expression higher that background
% % ------------------------------------------------------------------------------
noiseALL = noiseall';
%vertcat(noiseSUBJ{1}, noiseSUBJ{2}, noiseSUBJ{3}, noiseSUBJ{4}, noiseSUBJ{5}, noiseSUBJ{6});
% calculate the percentage of samples that each probe has expression value
% higher than a selected number
signalLevel = sum(noiseALL,1)./size(noiseALL,1);
indKeepProbes = find(signalLevel>signalThreshold);

%------------------------------------------------------------------------------
% ORIGINAL DATA OR FILTERED DATA
%------------------------------------------------------------------------------

if doOriginal
    genes = unique(DataTableProbe.EntrezID{1},'stable');
    listGenes = DataTableProbe.EntrezID{1};
    expression = Expressionall;
else
    genes = unique(DataTableProbe.EntrezID{1}(indKeepProbes),'stable');
    listGenes = DataTableProbe.EntrezID{1}(indKeepProbes);
    signalLevel = signalLevel(indKeepProbes);
    expression = Expressionall(indKeepProbes,:);
end

signalLevelProbes = cell(length(genes),1);
numberProbes = zeros(length(genes),1);

  for gene=1:length(genes)
        
        indGene = find(listGenes==genes(gene));
        signalLevelProbes{gene} = signalLevel(indGene);
        numberProbes(gene) = length(indGene);

  end

for numProbes=2
    if doEqual
        IDgene = genes(numberProbes==numProbes);
    else
        IDgene = genes(numberProbes>=numProbes);
    end
allData = cell(5,1);

for k=1:5
    type = probeSelection{k};
    data = cell(6,1);
    for sub=1:6
        load(sprintf('%s%s.mat', startFileName, type));
        data{sub} = expressionAll{sub};
    end
    
    % combine not-normalised data from all subjects and all samples
    allDataone = vertcat(data{1}, data{2}, data{3}, data{4}, data{5}, data{6});
    % load entrez IDs for genes that have more than 1 probe
    %load(sprintf('IDgenes%dplus.mat', numProbes));
    entrezIDs = probeInformation.EntrezID;
    % select those genes
    
    [~, keep] = intersect(entrezIDs, IDgene);
    
    allData{k} = allDataone(:,keep);
    
    
end

numGenes = size(allData{1},2);
avCorr = zeros(5,5);
coexpCorr = zeros(length(probeSelection),length(probeSelection));

for k=1:length(probeSelection)
    
    switch coexpressionOn
        case 'sample'
            mat = corr(allData{k}', 'type', 'Spearman');
        case 'gene'
            mat = corr(allData{k}, 'type', 'Spearman');
    end
    
    coexpValues{k} = mat(:);
    
end


% calculate correlation between each way of choosing a probe
for i=1:length(probeSelection)
    for j=i+1:length(probeSelection)
        
        correlation = zeros(numGenes,1);
        for g=1:numGenes
            correlation(g) = corr(allData{i}(:,g), allData{j}(:,g), 'type', 'Spearman');
        end
        
        
        avCorr(j,i) = median(correlation);
        coexpCorr(j,i) = corr(coexpValues{i}, coexpValues{j}, 'type', 'Spearman');
        
        
    end
end
figure; imagesc(avCorr);
set(gcf,'color','w');
nice_cmap = [make_cmap('steelblue',50,30,0);flipud(make_cmap('orangered',50,30,0))];
colormap(nice_cmap)
caxis([0 1])
%colormap([[1 1 1];BF_getcmap('reds',9)]); caxis([0 1]); colorbar;
xticks([1 2 3 4 5]); yticks([1 2 3 4 5]);
title(sprintf('Correlation between genes with %d probes (%d genes)', numProbes, numGenes));
xticklabels(probeSelection);
yticklabels(probeSelection);

figure; imagesc(coexpCorr);
set(gcf,'color','w');
nice_cmap = [make_cmap('steelblue',50,30,0);flipud(make_cmap('orangered',50,30,0))];
colormap(nice_cmap)
caxis([0 1])
%colormap([[1 1 1];BF_getcmap('reds',9)]); caxis([0 1]); colorbar;
xticks([1 2 3 4 5]); yticks([1 2 3 4 5]);
title(sprintf('Correlation between samples with %d probes (%d genes)', numProbes, numGenes));
xticklabels(probeSelection);
yticklabels(probeSelection);
end


