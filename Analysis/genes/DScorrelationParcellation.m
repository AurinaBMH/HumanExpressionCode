clear all; 

cd ('data/genes/processedData');

normMethod = 'scaledRobustSigmoid'; % '' for sigmoid; zscore for zscore;
probeSelection = {'Mean'};
onlyMultipleProbes = false;
percentDS = 5; 
numNodes = [82, 220, 360]; 
doNormalise = true; 
numProbes = 3; %(2 or 3)

load(sprintf('IDgenes%dplus.mat', numProbes));
DSscoresAll = cell(length(probeSelection),1);


for op=1:length(numNodes)
    
    load(sprintf('%dDS%d%s%s%d.mat', percentDS, numNodes(op), normMethod, probeSelection{1}, doNormalise)) % - generated using S5 script (probes chosen based on variance)
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
r = zeros(length(numNodes),length(numNodes)); 
p = zeros(length(numNodes),length(numNodes)); 
f=1; 
for i=1:length(numNodes)
    for j=i+1:length(numNodes)

        subplot(2,length(numNodes),f); scatter(DSscoresAll{i}, DSscoresAll{j}, sz, 'filled');
        [r(i,j),p(i,j)] = corr(DSscoresAll{i}, DSscoresAll{j}, 'type', 'Spearman'); 
        xlabel(sprintf('DS for probes selected based on %d node parcellation', numNodes(i))); 
        ylabel(sprintf('DS for probes selected based on %d node parcellation', numNodes(j)));
        f=f+1; 
    end
end
cd ../../..

