% make picture of ROI-gene matrix
load('100DS360scaledRobustSigmoidRNAseq1Lcortex_ROI_distCorr.mat')

numROIs = length(unique(SampleGeneExpression(:,1)));
ROIs = SampleGeneExpression(:,1);
uROIs = unique(ROIs);
ROIGeneExpression = zeros(numROIs, size(SampleGeneExpression,2));

for r=1:numROIs
    ro = uROIs(r);
    ROIind = find(ROIs==ro);
    if length(ROIind)>1
        ROIGeneExpression(r,:) = mean(SampleGeneExpression(ROIind,:));
    else
        ROIGeneExpression(r,:) = SampleGeneExpression(ROIind,:);
    end
end
% exclude the first column, so the indexes for genes match exactly
ROIGeneExpression = ROIGeneExpression(:,2:end);

[ordr,R,keepers] = BF_ClusterReorder(ROIGeneExpression, 'euclidean');
[ordc,R,keepers] = BF_ClusterReorder(ROIGeneExpression', 'euclidean');

[ordr1,R,keepers] = BF_ClusterReorder(ROIGeneExpression, 'corr');
[ordc1,R,keepers] = BF_ClusterReorder(ROIGeneExpression', 'corr');

A = ROIGeneExpression(ordr1,ordc1);

figure('color','w'); imagesc(ROIGeneExpression);
caxis([0,1]); colormap([flipud(BF_getcmap('blues',9));[1 1 1]; BF_getcmap('reds',9)]);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

figure('color','w'); imagesc(A);
caxis([0,1]); colormap([flipud(BF_getcmap('blues',9));[1 1 1]; BF_getcmap('reds',9)]);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

% make a plot for last supp figure
B = A(1:50, 5000:5100);
figure('color','w'); imagesc(B);
caxis([0,1]); colormap([flipud(BF_getcmap('blues',9));[1 1 1]; BF_getcmap('reds',9)]);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

