function [dMiddleNorm, dMiddle, PhiNormMean, PhiTrue, PhiRand]=PlotRichClub(Adj,DISTmean, WhatTypeNetwork,whatNullModel,numIter,numRepeats, scaling)
%-------------------------------------------------------------------------------
% Loads in data from a computation using networkRC,
% Plots RC curve and degree distribution
%-------------------------------------------------------------------------------

doRankSum = false; % Welch's t-test or ranksum test?
meanOrMedian = 'mean'; % mean distance

%-------------------------------------------------------------------------------
% Load in the data
%-------------------------------------------------------------------------------


if strcmp(WhatTypeNetwork, 'bu') || strcmp(WhatTypeNetwork, 'bd')
    Adj = logical(Adj);
end

if strcmp(WhatTypeNetwork, 'bd')
    [~,~,deg] = degrees_dir(Adj);
else
    deg = degrees_und(Adj);
end

if strcmp(WhatTypeNetwork, 'wuStrength') || strcmp(WhatTypeNetwork, 'wuStrengthBins')
    kmax = 100;
else
    kmax = max(deg);
end

pThreshold = 0.05;
%[~, E] = BF_PlotQuantiles(strength,strength,100);
[PhiNorm,PhiTrue,PhiRand] = RichClubPhiNorm(Adj,kmax, numIter,numRepeats,WhatTypeNetwork,whatNullModel); %, doBins);

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
%-------------------------------------------------------------------------------

f = figure('color','w');
axMain = subplot(4,1,2:3); hold on;
axMain.Position = [0.13 0.313 0.775 0.377];

%===============================================================================
% Distance curve
%===============================================================================
kRange = min(deg):max(deg);
% Add a baseline (1)
plot([min(deg),max(deg)],ones(2,1),':k');
dMean = zeros(length(kRange),2);
dAll = cell(length(kRange),2);
pVals = zeros(length(kRange),1);
for i = 1:length(kRange)
    isHub = double(deg > kRange(i));
    isRich = isHub'*isHub;
    isRich(eye(logical(size(isRich)))) = 0;
    notRich = ~isRich;
    notRich(eye(logical(size(isRich)))) = 0;
    dAll{i,1} = DISTmean(Adj & isRich);
    dAll{i,2} = DISTmean(Adj & notRich);
    
    % Difference
    if all(isnan(dAll{i,1})) || all(isnan(dAll{i,2}))
        pVals(i) = NaN;
    else
        if doRankSum
            pVals(i) = ranksum(dAll{i,1},dAll{i,2},'tail','right');
        else
            [~,pVals(i)] = ttest2(dAll{i,1},dAll{i,2},'Tail','right','Vartype','unequal');
        end
    end
end

switch meanOrMedian
    case 'mean'
        dMiddle = cellfun(@nanmean,dAll);
    case 'median'
        dMiddle = cellfun(@nanmedian,dAll);
end
dMiddleNorm = 1+scaling.*(dMiddle(:,1)-dMiddle(1,1))/(max(dMiddle(:,1))-dMiddle(1,1));

h_dCurve = plot(kRange,dMiddleNorm,'color',[49/255 130/255 189/255],'LineWidth',2.5);
sigD = pVals < 0.05;
plot(kRange(sigD),dMiddleNorm(sigD),'o','MarkerEdgeColor',[49/255 130/255 189/255],...
    'MarkerFaceColor',brighten([49/255 130/255 189/255],0.2),'LineWidth',1,'MarkerSize',9)

%===============================================================================
% Add rich-club curve:
%===============================================================================
h_RCcurve = plot(PhiNormMean, '-','Color','r','LineWidth',2.5);
% Significance at p = 0.05 as circles
plot(find(isSig),PhiNormMean(isSig),'o','MarkerEdgeColor','r',...
    'MarkerFaceColor',[1 .41 .38],'LineWidth',1,'MarkerSize',9)

% Add shaded rectangle
ylimNow = [axMain.YLim(1),axMain.YLim(2)];
%h_rect = rectangle('Position',[30,ylimNow(1),60,ylimNow(2)-ylimNow(1)+.1],...
% 'EdgeColor','none','FaceColor',ones(3,1)*0.90);
%uistack(h_rect,'bottom');
axMain.YLim = ylimNow;

% Set axis limits, legend, labels:
axMain.YLim = [0.95 max([max(PhiNormMean(isSig)),max(dMiddleNorm(sigD))])+0.05];
xlabel('Degree, k','fontsize',16);
ylabel('\Phi_{norm}','fontsize',16);
get(gca, 'XTick');
set(gca, 'FontSize', 16)
box off;
legend([h_RCcurve,h_dCurve],{'\Phi_{norm}','d'},'Location','NorthWest','fontsize',12);
legend('boxoff')
box('off')
axMain.XLim = [min(deg)-0.5,max(deg)+2];
axDD_type = subplot(4,1,1);

N = arrayfun(@(x)sum(deg==x),kRange);
a = bar(kRange,N,'EdgeColor',[0.35 0.35 0.35],'FaceColor',[0.35 0.35 0.35]); 

xticks([]); box off;
xlim([min(deg)-0.5 max(deg)+2]);
ylim([0 max(a.YData)+1]);

ylabel('Frequency','fontsize',16);
get(gca, 'YTick');
set(gca, 'FontSize', 16)
axDD_type.Position = [0.1300    0.776    0.775    0.149];
f.Position = [700,700,560,353];
box off;
end


