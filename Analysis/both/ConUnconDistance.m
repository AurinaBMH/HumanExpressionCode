% case study
% test connected VS unconnected coexpression (on corrected da
close all; 
clear all; 
parcellation = 'aparc+aseg';% 'aparc+aseg' , custom200ANDfslatlas20 'custom500ANDfslatlas20. HCPMMP1ANDfslatlas20
numCort = 250; 
weights = 'standard'; 
groupDens = [0.05 0.15 0.3 0.50]; 
tract = 'FACT'; % 'FACT, 'iFOD2'
distanceCorrection = 'NOdistCorr'; % 'distCorr', NOdistCorr
percDS = 100; 
groupConnectome = 'variance'; % 'length variance'

cd ('/Users/Aurina/GoogleDrive/Genetics_connectome/HumanExpression/data/connectomes')
connectomes = load(sprintf('%s_acpc_%s_SIFT2_%s_structnets.mat', parcellation, tract, weights));
%distances = load(sprintf('%s_acpc_%s_SIFT2_length_structnets.mat', parcellation, tract)); 
numSubjects = length(connectomes.SUBS); 
numNodes = size(connectomes.ADJS{1},1); 
cd ../../..
cd ('HumanExpression/data/genes/processedData')
load(sprintf('%dDS%dscaledRobustSigmoidRNAseq1Lcortex_ROI_%s.mat', percDS, numNodes+10, distanceCorrection)); 
% make A matrix and distances
conMatr = zeros(numNodes, numNodes, numSubjects); 
disMatr = zeros(numNodes, numNodes, numSubjects); 
tractLength = zeros(numNodes, numNodes, numSubjects); 
for s=1:numSubjects
    conMatr(:,:,s) = connectomes.ADJS{s};
    %tractLength(:,:,s) = distances.ADJS{s}; 
    disMatr(:,:,s) = pdist2(connectomes.COG{s},connectomes.COG{s}); 
end

conProb = mean(logical(conMatr),3); 
avdistances = mean(disMatr,3); 

%tractLength(tractLength==0) = NaN;
%avtractLength = nanmean(tractLength,3);

avWeight = conMatr; 
avWeight(avWeight==0)=NaN; 
avWeight = nanmean(avWeight,3); 

hemiid = zeros(numNodes,1); 
hemiid(1:numNodes/2) = 1; 
hemiid(numNodes/2+1:numNodes) = 2; 

for gd = 1:length(groupDens)
switch groupConnectome
    case 'variance'
        G = giveMeGroupAdj_variance(connectomes.ADJS, groupDens(gd)); 
    case 'lengthVariance'
        G = fcn_group_average(conMatr,avdistances,hemiid); 
    case 'consistency'
        G = giveMeGroupAdj_consistency(connectomes.ADJS, distances.ADJS ,0.95); 
end
% select only left cortex)
Glc = double(logical(G(1:numCort,1:numCort))); 

% adjcon = Glc; 
% adjcon(Glc==0) = NaN; 
% halfAdjcon = maskuHalf(adjcon); 
% 
% adjuncon = Glc; 
% adjuncon(adjuncon~=0) = NaN; 
% halfAdjuncon = maskuHalf(adjuncon); 
% halfAdjuncon(halfAdjuncon==0) = 1; 
% 
% Con = averageCoexpression.*halfAdjcon; Con = Con(:); Con(isnan(Con)) = []; 
% unCon = averageCoexpression.*halfAdjuncon; unCon = unCon(:); unCon(isnan(unCon)) = []; 
% 
% [h,p,~, stats] = ttest2(Con,unCon,'Vartype','unequal') %(Con,unCon, 'tail', 'right')
% 
% data{1} = Con; 
% data{2} = unCon; 
% JitteredParallelScatter(data)
% ylabel('Correlated gene expression'); 
% xticks([1 2]); xticklabels({'Connected', 'Unconnected'})
% set(gca, 'FontSize', 15); 

RichClubHuman(Glc,averageCoexpression);
end
% TEST(double(logical(G)),averageCoexpression);
% RichClubHuman(Glc,avWeight(1:180,1:180));
% RichClubHuman(Glc,avdistances(1:180,1:180));
% RichClubHuman(Glc,conProb(1:180,1:180));
% RichClubHuman(Glc,avtractLength(1:180,1:180));
