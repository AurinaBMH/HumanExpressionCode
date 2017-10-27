load('DSnewVariance.mat') % - generated using S5 script (probes chosen based on variance)
[a1,order] = sort(probeInformation.EntrezID); 
DSscores{1} = DS; probeSelection{1} = 'variance'; 
allData{1} = expSampNormalisedAll(:,order); 

load('DSnewPC.mat') % - generated using S5 script (probes chosen based on PC)
[a2,order] = sort(probeInformation.EntrezID); 
DSscores{2} = DS; probeSelection{2} = 'PC'; 
allData{2} = expSampNormalisedAll(:,order); 


load('DSnewLessNoise.mat')
[a3,order] = sort(probeInformation.EntrezID); 
DSscores{3} = DS; probeSelection{3} = 'noise'; 
allData{3} = expSampNormalisedAll(:,order); 


load('DSnewMean.mat');
[a4,order] = sort(probeInformation.EntrezID); 
DSscores{4} = DS; probeSelection{1} = 'Mean'; 
allData{4} = expSampNormalisedAll(:,order); 

load('DSnewrandom.mat');
[a5,order] = sort(probeInformation.EntrezID); 
DSscores{5} = DS; probeSelection{1} = 'random'; 
allData{5} = expSampNormalisedAll(:,order); 


numProbes = 2; 
coexpressionOn = 'sample'; 
DSthreshold = -1; 

%cd ('data/genes/processedData');
% entrezIDs for genes that have more than X probes. 
load(sprintf('IDgenes%dplus.mat', numProbes)); 
%load('compareProbesNormalised.mat')% - calculated using only cortex+subcortex to vilter noise (not 3702 samples) 
%load('DSvariance.mat'); % - DS values and probe information calculated on probes selected according to highest variance
%load('MicroarrayDataWITHcustMean360DistThresh2_CoordsAssigned.mat')
% compare probes contains normalised expression data in left cortex when
% different criteria are used to select a probe
% 1 - variance. 
% 2 - less noise
% 3 - mean of all probes
% 4 - highest loading on 1st PC
% it's generated by running all processing scripts and selecting normalised
% expression data combined from all subjects at the sample level. 


DSscores = DS>DSthreshold;
entrezIDs = a1; %probeInformation.EntrezID; 
% select those genes that have more than X probes (and maybe also high DS)
[~, keep] = intersect(entrezIDs, IDgene); 
%[~, keep] = intersect(entrezIDs(DSscores==1), IDgene); 
%  allData{1} = variance(:,2:end); 
%  allData{2} = noise(:,2:end);
%  allData{3} = mean(:,2:end);
%  allData{4} = pc(:,2:end);
%  allData{5} = random(:,2:end);

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