clear all; 

cd ('data/genes/processedData');

normMethod = 'scaledRobustSigmoid'; % '' for sigmoid; zscore for zscore;
probeSelection = {'Variance', 'LessNoise'};
onlyMultipleProbes = true; 
numNodes = 360; 
percentDS = 100; 
doNormalise = true; 

load('MicroarrayDataProbesUpdated.mat')
% select genes that have multiple probes, so thay will be sub-selected for
% comparison
[v, ind] = unique(DataTableProbe.EntrezID{1});

duplicate_ind = setdiff(1:size(DataTableProbe.EntrezID{1}, 1), ind);
duplicate_value = unique(DataTableProbe.EntrezID{1}(duplicate_ind));
% numProbes = 3; %(2 or 3)
% 
% load(sprintf('IDgenes%dplus.mat', numProbes));
% DSscoresAll = cell(length(probeSelection),1);
%ind = intersect(probeInformation.EntrezID, duplicate_value); 

for op=1:length(probeSelection)
    load(sprintf('%dDS%d%s%s%d.mat', percentDS, numNodes, normMethod, probeSelection{op}, doNormalise)) % - generated using S5 script (probes chosen based on variance)
    % reorder entrez IDs for each probe selection
    %[entrezIDs,order] = sort(probeInformation.EntrezID);
   % DSone = DS(order); 
    DS = probeInformation.DS; 
    if onlyMultipleProbes
        % find genes with multiple probes in the reordered
        % version and sub-select them. 
    [~, keep] = intersect(probeInformation.EntrezID, duplicate_value);
    DSscoresAll{op} = DS(keep); 
    
    else
        % id include all genes
    DSscoresAll{op} = DS; 
    
    end

end

% make a plot
sz=50; 
figure; title ('Correlation between DS scores'); 
set(gcf,'color','w');
scatter(DSscoresAll{1}, DSscoresAll{2}, sz, 'MarkerEdgeColor',[0.45 0.45 0.45],...
              'MarkerFaceColor',[.98 .85 .37]);
        [r,p] = corr(DSscoresAll{1}, DSscoresAll{2}, 'type', 'Spearman'); 
        xlabel(sprintf('DS for probes selected based on %s', probeSelection{1})); 
        ylabel(sprintf('DS for probes selected based on %s', probeSelection{2}));



% r = zeros(length(probeSelection),length(probeSelection)); 
% p = zeros(length(probeSelection),length(probeSelection)); 
% f=1; 
% for i=1:length(probeSelection)
%     for j=i+1:length(probeSelection)
% 
%         subplot(2,length(probeSelection),f); scatter(DSscoresAll{i}, DSscoresAll{j}, sz, 'filled');
%         [r(i,j),p(i,j)] = corr(DSscoresAll{i}, DSscoresAll{j}, 'type', 'Spearman'); 
%         xlabel(sprintf('DS for probes selected based on %s', probeSelection{i})); 
%         ylabel(sprintf('DS for probes selected based on %s', probeSelection{j}));
%         f=f+1; 
%     end
% end
cd ../../..

