% Rich club in the whole brain
type = 'HCP';

if strcmp(type, 'HCP')
    load('HCPMMP1ANDfslatlas20_acpc_connectome_data.mat')
elseif strcmp(type, 'GenCog')
    load('HCPMMP1ANDfslatlas20_GenCOG_connectome_data.mat');
end

[Adj, consist] = giveMeGroupAdj(standard);
%deg = degrees_und(Adj);

numIter = 5;
numRepeats = 10;
WhatTypeNetwork = 'bu'; 
whatNullModel = 'randmio_und'; %'randmio_und'; %'strength'; %

if strcmp(WhatTypeNetwork, 'bu')
    Adj = logical(Adj); 
end
onlyLeft = false; 
%Adj = logical(groupAdj(1:180,1:180));

if onlyLeft
    Adj = Adj(1:180,1:180); 
end

if strcmp(WhatTypeNetwork, 'wuStrength') || strcmp(WhatTypeNetwork, 'wuStrengthBins')
    kmax = 100; 
else
    kmax = max(sum(logical(Adj)));
end 
strength = sum(Adj); 
deg = degrees_und(Adj);
pThreshold = 0.05;
%[~, E] = BF_PlotQuantiles(strength,strength,100); 
[PhiNorm,PhiTrue,PhiRand] = RichClubPhiNorm(Adj,kmax, numIter,numRepeats,WhatTypeNetwork,whatNullModel); %, doBins);
figure;

% Compute p-values
pValues = zeros(kmax,1);
for i = 1:kmax
    pValues(i) = mean(PhiTrue(i) <= PhiRand(:,i));
end
% Significant phi
isSig = (pValues <= pThreshold);
PhiNormMean = zeros(size(PhiTrue));
for i = 1:length(PhiTrue)
    PhiNormMean(i) = PhiTrue(i)/mean(PhiRand(:,i));
end
% plot the graphs
if strcmp(WhatTypeNetwork, 'wuStrength') || strcmp(WhatTypeNetwork, 'wuStrengthBins')

        subplot(2,1,1); plot(PhiNormMean, '-','Color','r','LineWidth',2);  xlim([0 100]); ylim([min(PhiNormMean)-0.1 max(PhiNormMean)+0.05]); legend('Phi normalised', 'Location','NorthWest'); ...
            hold on;
        plot(find(isSig),PhiNormMean(isSig),'o','Color','r','LineWidth',3); xlim([0 100]); ylim([min(PhiNormMean)-0.1 max(PhiNormMean)+0.05]);
        title (sprintf('Normalised rich club \n%s', type)); xlabel('Strength'); ylabel('Phi (norm)'); set(gca,'FontSize',12,'fontWeight','bold');
        subplot(2,1,2); histogram(strength, 50);  title ('Strength distribution');  xlabel('Strength bins'); set(gca,'FontSize',12,'fontWeight','bold');
%     else
%histogram(strength, 50);
%         subplot(2,1,1); plot(PhiNormMean, '-','Color','r','LineWidth',2);  xlim([0 100]); 
%         ylim([min(PhiNormMean)-0.1 max(PhiNormMean)+0.05]); legend('Phi normalised', 'Location','NorthWest'); ...
%         hold on;
%         plot(find(isSig),PhiNormMean(isSig),'o','Color','r','LineWidth',3); xlim([min(E) max(E)]); ylim([min(PhiNormMean)-0.1 max(PhiNormMean)+0.05]);
%         title (sprintf('Normalised rich club \n%s', type)); xlabel('Strength'); ylabel('Phi (norm)'); set(gca,'FontSize',12,'fontWeight','bold');
        %subplot(2,1,2); histogram(deg, 50); title ('Degree distribution'); xlim([min(deg) max(deg)+2]); xlabel('Degree'); set(gca,'FontSize',12,'fontWeight','bold');

else

subplot(2,1,1); plot(PhiNormMean, '-','Color','r','LineWidth',2);  xlim([min(deg) max(deg)+2]); ylim([0.8 max(PhiNormMean)+0.05]); legend('Phi normalised', 'Location','NorthWest'); ...
hold on;
plot(find(isSig),PhiNormMean(isSig),'o','Color','r','LineWidth',3); xlim([min(deg) max(deg)+2]); ylim([0.8 max(PhiNormMean)+0.05]);
title (sprintf('Normalised rich club \n%s', type)); xlabel('Degree'); ylabel('Phi (norm)'); set(gca,'FontSize',12,'fontWeight','bold');
subplot(2,1,2); histogram(deg, 50); title ('Degree distribution'); xlim([min(deg) max(deg)+2]); xlabel('Degree'); set(gca,'FontSize',12,'fontWeight','bold');
end