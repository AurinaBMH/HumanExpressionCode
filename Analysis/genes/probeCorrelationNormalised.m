numProbes = 2; 
load(sprintf('IDgenes%dplus.mat', numProbes)); 

load('MicroarrayDataWITHcustMean360DistThresh2_CoordsAssigned.mat')
load('compareProbes.mat')
cd 'forBen'
load('geneSample_aparcasegUpdated.mat'); 
cd ..

DS = DSscores>0.5; 

entrezIDs = probeInformation.EntrezID; 
% select those genes
[~, keep] = intersect(entrezIDs(DS==1), IDgene); 

% allData{1} = variance(:,2:end); 
% allData{2} = noise(:,2:end);
% allData{3} = mean(:,2:end);
% allData{4} = pc(:,2:end);

allDatafinal = cell(4,1);
coexpValues = cell(4,1);
for k=1:4
    
allDatafinal{k} = allData{k}(:,keep); 
mat = corr(allDatafinal{k}, 'type', 'Spearman');
coexpValues{k} = mat(:); 

end

numGenes = size(allDatafinal{1},2); 
avCorr = zeros(4,4); 
coexpCorr = zeros(4,4); 
% calculate correlation between each way of choosing a probe
for i=1:4
    for j=i+1:4
        
       correlation = zeros(numGenes,1); 
        for g=1:numGenes
            correlation(g) = corr(allDatafinal{i}(:,g), allDatafinal{j}(:,g), 'type', 'Spearman'); 
        end
        
        avCorr(j,i) = mean(correlation); 
        coexpCorr(j,i) = corr(coexpValues{i}, coexpValues{j}, 'type', 'Spearman'); 
        
    end
end

figure; imagesc(avCorr); colormap([[1 1 1];BF_getcmap('reds',9)]); caxis([0 1]); colorbar; 
xticks([1 2 3 4]); yticks([1 2 3 4]); 
title(sprintf('Correlation between genes with %d and more probes (%d genes)', numProbes, numGenes)); 
xticklabels({'Variance','Less noise','Mean','maxPC'}); 
yticklabels({'Variance','Less noise','Mean','maxPC'}); 

figure; imagesc(coexpCorr); colormap([[1 1 1];BF_getcmap('reds',9)]); caxis([0 0.05]); colorbar; 
xticks([1 2 3 4]); yticks([1 2 3 4]); 
title(sprintf('Correlation between gene-gene coexpression matrices for genes with %d and more probes (%d genes)', numProbes, numGenes)); 
xticklabels({'Variance','Less noise','Mean','maxPC'}); 
yticklabels({'Variance','Less noise','Mean','maxPC'}); 
