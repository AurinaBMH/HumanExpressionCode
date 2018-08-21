function RichClubHuman(Adj,averageCoexpression,nodeData)
% ------------------------------------------------------------------------------
% Function plots coexpression for rich/feeder/peripheral lins as a function
% of degree using mean to summarise coexpression at each threshold
%-------------------------------------------------------------------------------
% ------------------------------------------------------------------------------

pThreshold = 0.05;
% ------------------------------------------------------------------------------
%% INPUTS:
% ------------------------------------------------------------------------------

% ------------- (1) analyzeWhat: what to plot the rich club curves for:
realLinkData = averageCoexpression;

% -------------- (3) labelNodesHow: how to define "rich"
labelNodesHow = 'hub-kth'; extraParam = {'degree',0};

numBins = 'all'; % Range of k to plot the rich club curve across

plotDist = false; % plot full distributions (at each k)

% ------------------------------------------------------------------------------
%% Assign data measured at each link in the network
% ------------------------------------------------------------------------------
% All directed connections:
linkedAdj = Adj;
% Add a mask to only include particular types of connections
numNeurons = size(linkedAdj,1);
nn = linspace(1,numNeurons,numNeurons);
allLinkData = realLinkData;
allLinkData(Adj==0) = NaN;

numNodes = length(Adj);

% ------------------------------------------------------------------------------
% nodeData (~degree) should not change with different nulls, which preserve the in/out degree of all nodes
% ------------------------------------------------------------------------------
if nargin<3
[~,~,nodeData] = AdjLabelNodes(labelNodesHow,Adj,extraParam,'bu');
else
    nodeData = nodeData;
end
%nodeData = degrees_und(Adj);
% ------------------------------------------------------------------------------
% Get groups of links based on their degree (or use all k)
% ------------------------------------------------------------------------------
sortK = sort(nodeData,'descend');
maxK = sortK(2); % Up to the second-highest k
if strcmp(numBins,'all')
    % All k are a bin
    kr = min(nodeData):maxK;
else % numeric number of bins
    kr = linspace(min(nodeData),maxK,numBins);
    kr = round(kr);
end
krAll = min(nodeData):max(nodeData);

% Proportion of nodes called a 'hub' through kr
propisHub = arrayfun(@(x) mean(nodeData > x),kr);
% ------------------------------------------------------------------------------
% Go through each class of links and compute statistics on the set of link
% data compared to the nulls:

whatLinks = {'rich','feeder','local'};
% Ben Fulcher, 2015-01-06
% Computes a t-test between special and non-special links for each k, and each link-type:
tStats = zeros(length(kr),length(whatLinks));
allHubHub = cell(1,1); % make a 1-component cell for consistency with null version
allHubHub{1} = cell(length(kr),length(whatLinks));
for i = 1:length(kr)
    for j = 1:length(whatLinks) % loop across rich, feedin, feedout, and local connections:
        
        % Keep links between high degree nodes that have linkData (i.e., links that exist):
        % (assumption that nodeData will be the same for all null networks)
        % Ben Fulcher, 2015-01-06 -- changed from >= to > to match actual rich-club coeff definition
        r = (nodeData > kr(i));
        
        keepMe = logical(zeros(numNodes,numNodes));
        % Add ones for a given type of link:
        switch j
            case 1 % 'rich'
                keepMe(r,r) = 1;
            case 2 % 'feeder'
                keepMe(~r,r) = 1;
                keepMe(r,~r) = 1;
            case 3 % 'local'
                keepMe(~r,~r) = 1;
                
        end
        % Remove missing data (NaNs encode no link):
        keepMe(isnan(allLinkData)) = 0;
        
        % Keep values assigned to each remaining link as this element of hubhubData
        linkDataSpecial = allLinkData(keepMe);
        allHubHub{1}{i,j} = linkDataSpecial;
        
        notSpecial = (~isnan(allLinkData) & ~keepMe);
        linkDataNotSpecial = allLinkData(notSpecial);
        
        % 2-sample t-test for special links greater than non-special links:

        [~,p] = ttest2(linkDataSpecial,linkDataNotSpecial,'Vartype','unequal', 'Tail','right');
        
        tStats(i,j) = p;
    end
end

% ------------------------------------------------------------------------------
%% Plot as rich plots
% ------------------------------------------------------------------------------
myColors = GiveMeColors('RFPU'); % [BF_getcmap('spectral',4,1),BF_getcmap('set2',4,1)];
plotOnOne = true; % plot all on one figure
includeHist = true;
plotJustRich = true;
includeStd = false;
plotNulls = false; % plot each null trajectory in the figure
sigThresh = 0.05;

for j = 1:length(whatLinks)
    if plotJustRich && (j < 1)
        break
    end
    if plotDist
        JitteredParallelScatter(allHubHub{1}(:,j))
    else
        if includeHist
            if (~plotOnOne || (plotOnOne==1 && j==1))
                figure('color','w');
                sp=subplot(5,3,1:6);
                if islogical(Adj)
                    
                    N = arrayfun(@(x)sum(nodeData==x),krAll);
                    bar(krAll,N,'EdgeColor',[0.35 0.35 0.35],'FaceColor',[0.35 0.35 0.35])
                else
                    histogram(nodeData,70,'EdgeColor',[0.35 0.35 0.35],'FaceColor',[0.35 0.35 0.35]);
                end
                xlim([min(nodeData)-0.8,max(nodeData)+0.8]);
                %xlim([min(nodeData)-0.5,max(nodeData)+0.5]);
                xticks([]); box off;
                ylabel('Frequency', 'FontSize', 18)
                get(gca, 'YTick');
                set(gca, 'FontSize', 16)
                sp=subplot(5,3,7:15);
                hold on;
                
            end
        else
            if (~plotOnOne || (plotOnOne==1 && j==1))
                figure('color','w'); hold on
            end
        end
    end

    
    % Plot flat line for rich
    if j==1
        
        plot([kr(1),krAll(end)],ones(2,1)*nanmean(allHubHub{1}{1,j}),':','color','k','LineWidth',3)
        
    end
    
    realTrajectory = cellfun(@nanmean,allHubHub{1}(:,j));
    getYlim(:,j) = cellfun(@nanmean,allHubHub{1}(:,j)); 
    xlim([min(nodeData)-0.5,max(nodeData)+0.5]);
    %ylim([0, max(getYlim(:))]);
    xlabel(extraParam{1})
    
    
    realStd = cellfun(@std,allHubHub{1}(:,j));
    
    % p-values from 2-sample t-test with unequal variances:
    pvalues = tStats(:,j);
    
    isSig = (pvalues < sigThresh); % significantly higher than null
    
    % mean (real data trajectory):
    lineStyle = '-'; markerStyle = 'o';
    
    if any(isSig)
        plot(kr(isSig),realTrajectory(isSig),markerStyle,'MarkerEdgeColor',myColors(j,:),...
            'MarkerFaceColor',brighten(myColors(j,:),+0.5),'LineWidth',1,'MarkerSize',9)
    end
    
    % mean trajectory:
    plot(kr,realTrajectory,lineStyle,'color',myColors(j,:),'LineWidth',3)
    
    % +/- std:
    if includeStd
        plot(kr,realTrajectory+realStd,lineStyle,'color',myColors(j,:),'LineWidth',1)
        plot(kr,realTrajectory-realStd,lineStyle,'color',myColors(j,:),'LineWidth',1)
    end
    
    xLimits = get(gca,'xlim'); yLimits = get(gca,'ylim');
    
    if ~plotJustRich
        text(xLimits(1)+0.1*diff(xLimits),yLimits(1)+0.9*diff(yLimits)-j/20,whatLinks{j},'color',myColors(j,:),'FontSize',18)
    end

    divisionText = '';
    
    
end

ax = gca;
if plotJustRich
    ax.FontSize = 18;
else
    set(gcf, 'Position', [500 500 750 500])
end


% Add light gray rectangle
%ylimNow = [ax.YLim(1),ax.YLim(2)];
% if ax.YLim(1)<0
%     h_rect = rectangle('Position',[42,ylimNow(1),12,ylimNow(2)+(ylimNow(1).*(-1))],'EdgeColor','none','FaceColor',ones(3,1)*0.90);
% else
%     h_rect = rectangle('Position',[42,ylimNow(1),12,ylimNow(2)-ylimNow(1)],'EdgeColor','none','FaceColor',ones(3,1)*0.90);
% end
%uistack(h_rect,'bottom');
%set(gca,'ylim',ylimNow);
set(gcf, 'Position', [500 500 750 500])

axisName = {'Mean correlated', 'gene expression'};
ylabel(axisName, 'FontSize', 18)
xlabel('Node degree, k','FontSize', 18);

sp=subplot(5,3,7:15);
pos=get(sp,'Position');
deg = degrees_und(Adj);
kRange = min(deg):max(deg);

end


