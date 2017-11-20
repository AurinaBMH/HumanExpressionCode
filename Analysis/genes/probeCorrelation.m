useCUSTprobes = true;
probeSelection = {'Variance','PC','LessNoise', 'Mean', 'random'};
numProbes = 2;
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
    load(sprintf('IDgenes%d.mat', numProbes));
    entrezIDs = probeInformation.EntrezID;
    % select those genes
    %[a, indord] = sort(ProbeInformation.EntrezID);
   % allData{k} = allDataone{k}(:,indord);
    
    [~, keep] = intersect(entrezIDs, IDgene);
   
    allData{k} = allDataone(:,keep);
    
    
end

numGenes = size(allData{1},2);
avCorr = zeros(5,5);
% calculate correlation between each way of choosing a probe
for i=1:5
    for j=i+1:5
        
        correlation = zeros(numGenes,1);
        for g=1:numGenes
            correlation(g) = corr(allData{i}(:,g), allData{j}(:,g), 'type', 'Spearman');
        end
        
        avCorr(j,i) = mean(correlation);
        
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
xticklabels({'Variance','Less noise','Mean','PC', 'Random'});
yticklabels({'Variance','Less noise','Mean','PC', 'Random'});


