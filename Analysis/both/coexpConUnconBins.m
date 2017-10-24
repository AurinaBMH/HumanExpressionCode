% compare coexpression for connected vs unconnected links

clear all; close all; 

cd ('data/connectomes')
load('HCPgroupConnectomes.mat')
cd ..
cd ('genes/processedData')
load('5DSuncorrected.mat')
load('HCPdistances.mat')

Adj = double(logical(AdjCount(1:180, 1:180)));
numThresholds = 16; 

%% ---------------------------------------------------------------------------------
% filter data in bins
[xThresholds,yMeans] = BF_PlotQuantiles(distPlot(:),averageCoexpression(:),numThresholds,false,true);

nanMat = nan(180,180);
nanMat = tril(nanMat);
nanMat(nanMat==0) = 1;

distPlot = distPlot.*nanMat;
averageCoexpression = averageCoexpression.*nanMat;

% bin data into quantiles
%% ---------------------------------------------------------------------------------
figure; 

for i=2:length(xThresholds)
 
    ind = distPlot<=xThresholds(i) & distPlot>xThresholds(i-1);
    
    % map onto connectivity data
    maskCon = Adj & ind;
    maskUncon = ~Adj & ind;
    
    coexpCon = nonzeros(averageCoexpression.*maskCon);
    coexpCon(isnan(coexpCon)) = [];
    coexpUncon = nonzeros(averageCoexpression.*maskUncon);
    coexpUncon(isnan(coexpUncon)) = [];
    
    dataCell{1} = coexpCon;
    dataCell{2} = coexpUncon;
    [h,p(i),ci,stats] = ttest2(dataCell{1}, dataCell{2});
    subplot(5,3,i-1); 
    JitteredParallelScatter(dataCell, true, true, false); ylim([-0.5 0.5]); 
    title(sprintf('Connected VS unconnected in bin %d', i-1));
    ylabel('Coexpression');
    set(gca,'Xtick', [1 2], 'XTickLabel',{sprintf('Connected (%d pairs)', length(coexpCon)), ...
    sprintf('Unconnected (%d pairs)',length(coexpUncon))}); 

end  
cd ../../..
        