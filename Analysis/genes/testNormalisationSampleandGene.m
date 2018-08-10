
%% load data normalised per sample and per gene
load('100DS360scaledRobustSigmoidNormSampleGenesRNAseq1Lcortex_ROI_NOdistCorrEuclidean.mat')
distMatrix = pdist2(SampleCoordinates(:,2:end),SampleCoordinates(:,2:end)); 
data = SampleGeneExpression(:,2:end);

ordSamples = BF_ClusterReorder(data,distMatrix);
dataORD = data(ordSamples,:); 

ordc1 = BF_ClusterReorder(dataORD', 'corr');
dataORDORDnsg = dataORD(:,ordc1); 

figure('color','w'); imagesc(dataORDORDnsg); title('Normalised samples and genes'); 
caxis([0,1]); colormap([flipud(BF_getcmap('blues',9));[1 1 1]; BF_getcmap('reds',9)]);
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
          
genes = probeInformation.EntrezID; 

%% load normalised data per-gene across samples
load('100DS360scaledRobustSigmoidNormGenesRNAseq1Lcortex_ROI_NOdistCorrEuclidean.mat')
% plot it reordered by inter-sample distance
[~,indGenes] = intersect(probeInformation.EntrezID, genes); 
data = SampleGeneExpression(:,2:end);
data = data(:, indGenes); 

dataORD = data(ordSamples,:); 

% use the same ordering as before to compare
dataORDORDng = dataORD(:,ordc1); 

figure('color','w'); imagesc(dataORDORDng); title('Normalised genes'); 
caxis([0,1]); colormap([flipud(BF_getcmap('blues',9));[1 1 1]; BF_getcmap('reds',9)]);
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])

%% calculate the correlation between samples

for s = 1:size(dataORDORDng,1)
     rs(s) = corr(dataORDORDng(s,:)', dataORDORDnsg(s,:)', 'rows', 'complete'); 
end

figure; histogram(rs, 100); title('Correlations between samples')
fprintf('Median correlation between samples is %d, IQR %d\n', median(rs), iqr(rs))
%% calculate the correlation between genes

for g=1:size(dataORDORDng,2)
    rg(g) = corr(dataORDORDng(:,g), dataORDORDnsg(:,g)); 
end

figure; histogram(rg, 100); title('Correlations between genes')
fprintf('Median correlation between genes is %d, IQR %d\n', median(rg), iqr(rg))











% load('100DS360scaledRobustSigmoidRNAseq1Lcortex_ROI_NOdistCorr.mat')
% % compare expression values (correlation) for only gene normalisation
%  
% 
% % and additional sample normalisation
% % normalise per sample
% dataNormGene = SampleGeneExpression(:,2:end);
% % reorder data based on distance between samples and see if samples close
% % to each other have very different expression
% distMatrix = pdist2(SampleCoordinates(:,2:end),SampleCoordinates(:,2:end)); 
% [ordSamples,R,keepers] = BF_ClusterReorder(dataNormGene,distMatrix);
% 
% dataNormSampleORD = dataNormGene(ordSamples,:); 
% [ordc1,R,keepers] = BF_ClusterReorder(dataNormSampleORD', 'corr');
% dataNormGeneORDORD = dataNormSampleORD(:,ordc1); 
% 
% figure('color','w'); imagesc(dataNormGeneORDORD);
% caxis([0,1]); colormap([flipud(BF_getcmap('blues',9));[1 1 1]; BF_getcmap('reds',9)]);
%         set(gca,'xtick',[])
%         set(gca,'xticklabel',[])
%         set(gca,'ytick',[])
%         set(gca,'yticklabel',[])
%         
%         
% % now plot the same thing after per-sample normalisation
% 
% dataNormSample = (BF_NormalizeMatrix(dataNormGene', 'scaledRobustSigmoid'))'; 
% dataNormSample2 = (BF_NormalizeMatrix(dataNormSample, 'scaledRobustSigmoid')); 
% dataNormSampleORD = dataNormSample2(ordSamples,:); 
% 
% [ordc2,R,keepers] = BF_ClusterReorder(dataNormSampleORD', 'corr');
% dataNormSampleORDORD = dataNormSampleORD(:,ordc1); 
% 
% figure('color','w'); imagesc(dataNormSampleORDORD);
% caxis([0,1]); colormap([flipud(BF_getcmap('blues',9));[1 1 1]; BF_getcmap('reds',9)]);
%         set(gca,'xtick',[])
%         set(gca,'xticklabel',[])
%         set(gca,'ytick',[])
%         set(gca,'yticklabel',[])
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % correlated expression between norm and non-binrm samples
% figure; scatter(dataNormGene(:,1), dataNormSample(:,1)); 
% xlabel('Norm per gene'); ylabel('Norm per sample'); 
% for g=1:10
%     r(g) = corr(dataNormGene(:,g), dataNormSample(:,g)); 
%     figure; scatter(dataNormGene(:,g), dataNormSample(:,g)); 
%     xlabel('Norm per gene'); ylabel('Norm per sample'); 
% end
% 
% figure; histogram(r,100); 
% 
% for i=1:100
%     ind = find(SampleGeneExpression(:,1)==i); 
%     if ~isempty(ind) && length(ind)>1
% 
%     meanExp(i,:) = mean(SampleGeneExpression(ind,:)); 
%     meanCord(i,:) = mean(SampleCoordinates(ind,:)); 
%     elseif ~isempty(ind) && length(ind)==1
%     meanExp(i,:) = SampleGeneExpression(ind,:); 
%     meanCord(i,:) = SampleCoordinates(ind,:); 
%         
%     else
%         meanExp(i,:) = NaN; 
%         meanCord(i,:) = NaN; 
%     end
%         
% end
% 
% meanExp1 = mean(meanExp,2); 
% 
% 
% % calculate distance
% dist = pdist2(meanCord(:,2:end), meanCord(:,2:end)); 
% for g=2:10
% figure; 
% scatter3(meanCord(:,2),meanCord(:,3), meanCord(:,4), meanExp(:,g)*100); 
% end
% 
% figure; 
% scatter3(meanCord(:,2),meanCord(:,3), meanCord(:,4), meanExp1*100); 
% 
% 
% clear all; 
% load('MicroarrayDataWITHcustProbesUpdatedXXXRNAseq220DistThresh2.mat')
% 
% for s=1:6
%     dataSamp{s} = BF_NormalizeMatrix(DataExpression{s}(:,3:end)', 'scaledRobustSigmoid')'; 
%     dataGene{s} = BF_NormalizeMatrix(DataExpression{s}(:,3:end), 'scaledRobustSigmoid'); 
%     dataSampGene{s} = BF_NormalizeMatrix(dataSamp{s}, 'scaledRobustSigmoid'); 
% end
% 
% dataNormSamp = vertcat(dataSamp{1}, dataSamp{2}, dataSamp{3}, ...
%     dataSamp{4}, dataSamp{5}, dataSamp{6}); 
% 
% % like always
% dataNormGene = vertcat(dataGene{1}, dataGene{2}, dataGene{3}, ...
%     dataGene{4}, dataGene{5}, dataGene{6}); 
% 
% % sample and then gene
% dataNormSampleGene = vertcat(dataSampGene{1}, dataSampGene{2}, dataSampGene{3}, ...
%     dataSampGene{4}, dataSampGene{5}, dataSampGene{6}); 
% 
% % correlations between genes
% 
% for g=1:size(dataNormGene,2)
%     rg(g) = corr(dataNormGene(:,g), dataNormSampleGene(:,g)); 
% end
% 
% % correlation between samples
% for s=1:size(dataNormGene,1)
%     rs(s) = corr(dataNormGene(s,:)', dataNormSampleGene(s,:)', 'rows', 'complete'); 
% end
% 
% figure; histogram(rg, 100); title('Correlation between genes'); 
% figure; histogram(rs, 100); title('Correlation between samples'); 
% 


