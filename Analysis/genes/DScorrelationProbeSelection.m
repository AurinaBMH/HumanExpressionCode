clear all; 

cd ('data/genes/processedData');

normMethod = 'scaledRobustSigmoid'; % '' for sigmoid; zscore for zscore;
probeSelection = {'Variance', 'LessNoise', 'Mean'};
onlyMultipleProbes = false; 
numNodes = 220; 
percentDS = 5; 
doNormalise = true; 
numProbes = 3; %(2 or 3)

load(sprintf('IDgenes%dplus.mat', numProbes));
DSscoresAll = cell(length(probeSelection),1);


for op=1:length(probeSelection)
    load(sprintf('%dDS%d%s%s%d.mat', percentDS, numNodes, normMethod, probeSelection{op}, doNormalise)) % - generated using S5 script (probes chosen based on variance)
    % reorder entrez IDs for each probe selection
    %[entrezIDs,order] = sort(probeInformation.EntrezID);
   % DSone = DS(order); 
    DS = probeInformation.DS; 
    if onlyMultipleProbes
        % find genes with multiple probes in the reordered
        % version and sub-select them. 
    [~, keep] = intersect(probeInformation.EntrezID, IDgene);
    DSscoresAll{op} = DS(keep); 
    
    else
        % id include all genes
    DSscoresAll{op} = DS; 
    
    end

end

% make a plot
sz=10; 
figure; title ('Correlation between DS scores'); 
set(gcf,'color','w');
r = zeros(length(probeSelection),length(probeSelection)); 
p = zeros(length(probeSelection),length(probeSelection)); 
f=1; 
for i=1:length(probeSelection)
    for j=i+1:length(probeSelection)

        subplot(2,length(probeSelection),f); scatter(DSscoresAll{i}, DSscoresAll{j}, sz, 'filled');
        [r(i,j),p(i,j)] = corr(DSscoresAll{i}, DSscoresAll{j}, 'type', 'Spearman'); 
        xlabel(sprintf('DS for probes selected based on %s', probeSelection{i})); 
        ylabel(sprintf('DS for probes selected based on %s', probeSelection{j}));
        f=f+1; 
    end
end
cd ../../..

