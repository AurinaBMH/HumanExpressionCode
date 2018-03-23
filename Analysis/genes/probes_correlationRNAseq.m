
% script to evaluate correlations between RNAseq and microarray data and
% make list for enrichment
load('MicroarrayDataWITHcustProbesUpdatedXXXRNAseqnoQC.mat')
ind = cell2mat(cellfun(@(x)any(isnan(x)),avgCorr,'UniformOutput',false)); 
fprintf('%d genes are overlapping between RNA-seq and microarray datasets\n', length(find(ind==0)))

avgCorr = avgCorr(ind==0); 

for i=1:length(find(ind==0))
    maxCor(i) = max(avgCorr{i}); 
end

maxCor = maxCor'; 

figure; colors = [.96 .63 .55; 1 .46 .22]; 
histogram(maxCor, 100,'EdgeColor',[.6 .6 .6],...
    'FaceColor',colors(1,:)); 
xlabel('Average maximum correlation between microarray and RNA-seq expression','FontSize', 14)
ylabel('Number of genes','FontSize', 14')
set(gcf,'color','w'); hold on; 
legendText = sprintf('%d genes', length(find(ind==0))); 
legend(legendText); 