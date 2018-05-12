% case study
% test connected VS unconnected coexpression (on corrected da
close all; 
cd ('/Users/Aurina/GoogleDrive/Genetics_connectome/Heritability/data/connectomes')
load('HCPMMP1ANDfslatlas20_acpc_FACT_SIFT2_standard_structnets.mat')
cd ../../..
cd ('HumanExpression/data/genes/processedData')
load('30DS360scaledRobustSigmoidRNAseq1Lcortex_ROI_NOdistCorr_PC.mat'); 
load('30DS360scaledRobustSigmoidRNAseq1Lcortex_ROI_NOdistCorr_MEAN.mat'); 

[groupAdj, consist] = giveMeGroupAdj_variance(ADJS);
adj = double(logical(groupAdj(1:180, 1:180))); 
adjcon = adj; 
adjuncon = adj; 

adjcon(adjcon==0) = NaN; 
halfAdjcon = maskuHalf(adjcon); 

adjuncon(adjuncon~=0) = NaN; 
halfAdjuncon = maskuHalf(adjuncon); 
halfAdjuncon(halfAdjuncon==0) = 1; 


Con = averageCoexpression.*halfAdjcon; Con = Con(:); Con(isnan(Con)) = []; 
unCon = averageCoexpression.*halfAdjuncon; unCon = unCon(:); unCon(isnan(unCon)) = []; 

[h,p,~, stats] = ttest2(Con,unCon,'Vartype','unequal') %(Con,unCon, 'tail', 'right')

data{1} = Con; 
data{2} = unCon; 
JitteredParallelScatter(data)

RichClubHuman(logical(groupAdj(1:180, 1:180)),averageCoexpression);

load('100DS360scaledRobustSigmoidRNAseq1Lcortex_ROI_NOdistCorr_PC.mat'); 

%[groupAdj, consist] = giveMeGroupAdj_variance(ADJS);
adj = double(logical(groupAdj(1:180, 1:180))); 
adjcon = adj; 
adjuncon = adj; 

adjcon(adjcon==0) = NaN; 
halfAdjcon = maskuHalf(adjcon); 

adjuncon(adjuncon~=0) = NaN; 
halfAdjuncon = maskuHalf(adjuncon); 
halfAdjuncon(halfAdjuncon==0) = 1; 


Con = averageCoexpression.*halfAdjcon; Con = Con(:); Con(isnan(Con)) = []; 
unCon = averageCoexpression.*halfAdjuncon; unCon = unCon(:); unCon(isnan(unCon)) = []; 

[h,p, ~, stats] = ttest2(Con,unCon,'Vartype','unequal')%; [p,h, stats] = ranksum(Con,unCon, 'tail', 'right')

data{1} = Con; 
data{2} = unCon; 
JitteredParallelScatter(data)

RichClubHuman(logical(groupAdj(1:180, 1:180)),averageCoexpression);