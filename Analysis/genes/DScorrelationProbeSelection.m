clear all; 

cd ('data/genes/processedData');
load('MicroarrayDataWITHcustProbesUpdatedXXX.mat')

options.ExcludeCBandBS =  true;
options.useCUSTprobes = true;
options.updateProbes = 'reannotator'; %'Biomart', 'reannotator', 'no'; 
options.probeSelections = {'RNAseq'};
options.parcellations = {'cust100'};
options.signalThreshold = 0.5; 
options.RNAseqThreshold = 0.2; 
options.correctDistance = true; 
options.onlyMultipleProbes = true; 
options.calculateDS = true;
options.distanceCorrection = 'Euclidean'; 
options.coexpressionFor = 'all';
options.Fit = {'removeMean'};
options.distanceThreshold = 2; % first run 30, then with the final threshold 2
options.normaliseWhat = 'Lcortex';
options.normMethod = 'scaledRobustSigmoid'; 
options.percentDS =  10;
options.doNormalise = true;
options.resolution = 'ROI'; 

normMethod = options.normMethod; % '' for sigmoid; zscore for zscore;
probeSelection = {'Variance', 'LessNoise'};
onlyMultipleProbes = true; 
numNodes = 360; 
percentDS = options.percentDS; 
doNormalise = options.doNormalise;
normaliseWhat = options.normaliseWhat; 
% select genes that have multiple probes, so thay will be sub-selected for
% comparison
signalLevel = sum(noiseall,2)./size(noiseall,2);
indKeepProbes = find(signalLevel>=options.signalThreshold);

% sduplicated values calculated only on those probes that pass QC
% threshold;
[v, ind] = unique(DataTableProbe.EntrezID{1}(indKeepProbes));
duplicate_ind = setdiff(1:size(DataTableProbe.EntrezID{1}(indKeepProbes), 1), ind);
entrezSelected = DataTableProbe.EntrezID{1}(indKeepProbes);
duplicate_value = unique(entrezSelected(duplicate_ind));

% numProbes = 3; %(2 or 3)
% 
% load(sprintf('IDgenes%dplus.mat', numProbes));
% DSscoresAll = cell(length(probeSelection),1);
%ind = intersect(probeInformation.EntrezID, duplicate_value); 

for op=1:length(probeSelection)
    load(sprintf('%dDS%d%s%s%d%s', percentDS, numNodes, normMethod, probeSelection{op}, doNormalise, normaliseWhat)) % - generated using S5 script (probes chosen based on variance)
    % reorder entrez IDs for each probe selection
    %[entrezIDs,order] = sort(probeInformation.EntrezID);
   % DSone = DS(order); 
   [genesMultiple, keep] = intersect(probeInformation.EntrezID, duplicate_value);
    DS = probeInformation.DS; 
    if onlyMultipleProbes
        % find genes with multiple probes in the reordered
        % version and sub-select them. 
    %[~, keep] = intersect(probeInformation.EntrezID, duplicate_value);
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

