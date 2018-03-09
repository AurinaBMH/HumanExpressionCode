clear all; 

cd ('data/genes/processedData');


%numProbes = 3; %(2 or 3)

options.ExcludeCBandBS =  true;
options.useCUSTprobes = true; 
options.probeSelections = {'RNAseq'};
options.parcellations = {'cust100'};
options.onlyMultipleProbes = false; 
options.distanceThreshold = 2; % first run 30, then with the final threshold 2
options.normaliseWhat = 'Lcortex';
options.normMethod = 'scaledRobustSigmoid'; 
options.percentDS =  10;
options.doNormalise = true;
options.resolution = 'ROI'; 

normMethod = options.normMethod;  % '' for sigmoid; zscore for zscore;
probeSelection = options.probeSelections; 
onlyMultipleProbes = options.onlyMultipleProbes; 
percentDS = options.percentDS; 
numNodes = [82, 360]; 
doNormalise = options.doNormalise; 

load('MicroarrayDataWITHcustProbesUpdatedXXX.mat')
[v, ind] = unique(DataTableProbe.EntrezID{1});

duplicate_ind = setdiff(1:size(DataTableProbe.EntrezID{1}, 1), ind);
duplicate_value = unique(DataTableProbe.EntrezID{1}(duplicate_ind));

%load(sprintf('IDgenes%dplus.mat', numProbes));
DSscoresAll = cell(length(probeSelection),1);


for op=1:length(numNodes)
    
    load(sprintf('%dDS%d%s%s%d%s', percentDS, numNodes(op), normMethod, probeSelection{1}, doNormalise, normaliseWhat)) %  % - generated using S5 script (probes chosen based on variance)
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
              'MarkerFaceColor',[.53 .83 .97]); 
        [r,p] = corr(DSscoresAll{1}, DSscoresAll{2}, 'type', 'Spearman'); 
        xlabel(sprintf('DS for probes selected based on %d node parcellation', numNodes(1))); 
        ylabel(sprintf('DS for probes selected based on %d node parcellation', numNodes(2)));
        
        
% r = zeros(length(numNodes),length(numNodes)); 
% p = zeros(length(numNodes),length(numNodes)); 
% f=1; 
% for i=1:length(numNodes)
%     for j=i+1:length(numNodes)
% 
%         subplot(2,length(numNodes),f); scatter(DSscoresAll{i}, DSscoresAll{j}, sz, 'filled');
%         [r(i,j),p(i,j)] = corr(DSscoresAll{i}, DSscoresAll{j}, 'type', 'Spearman'); 
%         xlabel(sprintf('DS for probes selected based on %d node parcellation', numNodes(i))); 
%         ylabel(sprintf('DS for probes selected based on %d node parcellation', numNodes(j)));
%         f=f+1; 
%     end
% end
cd ../../..

