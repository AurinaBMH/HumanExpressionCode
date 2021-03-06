% How many genes have multiple probes
%------------------------------------------------------------------------------
% Choose options
%------------------------------------------------------------------------------
clear all; 
useCUSTprobes = true;
signalThreshold = 0.5; % percentage of samples that a selected probe has expression levels that are higher than background
doOriginal = false; %false;
%------------------------------------------------------------------------------
% Load the data
%------------------------------------------------------------------------------
cd ('data/genes/processedData');

if useCUSTprobes
    fprintf(1,'Loading the data with CUST probes and assigning variables\n')
    startFileName = 'MicroarrayDataWITHcust';
else
    fprintf(1,'Loading the data without CUST probes and assigning variables\n')
    startFileName = 'MicroarrayData';
end

load(sprintf('%s.mat', startFileName));

% % ------------------------------------------------------------------------------
% % Find best representative in a set of duplicates using maxVar and remove all others:
% % ------------------------------------------------------------------------------

ProbeID = DataTableProbe.ProbeID{1,1};
% % ------------------------------------------------------------------------------
% % First, find probes that have very noisy data and remove them from consideration
% % Threshold for removing those proges is defined as the percentage of
% % samples a probe has expression higher that background
% % ------------------------------------------------------------------------------
noiseALL = noiseall'; 
%vertcat(noiseSUBJ{1}, noiseSUBJ{2}, noiseSUBJ{3}, noiseSUBJ{4}, noiseSUBJ{5}, noiseSUBJ{6});
% calculate the percentage of samples that each probe has expression value
% higher than a selected number
signalLevel = sum(noiseALL,1)./size(noiseALL,1);
indKeepProbes = find(signalLevel>=signalThreshold);

%------------------------------------------------------------------------------
% ORIGINAL DATA OR FILTERED DATA
%------------------------------------------------------------------------------

if doOriginal
    genes = unique(DataTableProbe.EntrezID{1},'stable');
    listGenes = DataTableProbe.EntrezID{1};
    expression = Expressionall; 
else
    genes = unique(DataTableProbe.EntrezID{1}(indKeepProbes),'stable');
    listGenes = DataTableProbe.EntrezID{1}(indKeepProbes);
    signalLevel = signalLevel(indKeepProbes);
    expression = Expressionall(indKeepProbes,:);
end

signalLevelProbes = cell(length(genes),1);
numberProbes = zeros(length(genes),1);
k=1; 
l=1; 
m=1; 
for gene=1:length(genes)
    
    indGene = find(listGenes==genes(gene));
    signalLevelProbes{gene} = signalLevel(indGene);
    numberProbes(gene) = length(indGene);
    if numberProbes(gene)==2
        [r(k),p(k)] = corr(expression(indGene(1),:)', expression(indGene(2),:)', 'type', 'Spearman'); 
        variance12(k) = var(expression(indGene(1),:)'); 
        variance22(k) = var(expression(indGene(2),:)'); 
        G{k} = indGene; 
        k=k+1; 
    end
    if numberProbes(gene)==3
        [r1(l),p1(l)] = corr(expression(indGene(1),:)', expression(indGene(2),:)'); 
        [r2(l),p2(l)] = corr(expression(indGene(1),:)', expression(indGene(3),:)'); 
        [r3(l),p3(l)] = corr(expression(indGene(2),:)', expression(indGene(3),:)'); 
        variance13(l) = var(expression(indGene(1),:)'); 
        variance23(l) = var(expression(indGene(2),:)'); 
        variance33(l) = var(expression(indGene(3),:)'); 
        l=l+1;  
    end
    
    % save entrez ID for genes with multiple probes to be used to check the
    % influence of probe selection
    if numberProbes(gene)>4
        IDgene(m) = listGenes(indGene(1));
         m=m+1; 
    end
   
    
end

% [nrProbesSummary(:,1),nrProbesSummary(:,2)] = hist(numberProbes,unique(numberProbes)); 
% 
% 
% nice_cmap = [make_cmap('steelblue',50,30,0);flipud(make_cmap('orangered',50,30,0))];
% figure; histogram(rall,100, 'Normalization', 'probability', 'facecolor',nice_cmap(98,:),'facealpha',.5,'edgecolor',nice_cmap(100,:)); 
% hold on; ...
%     histogram(r, 100, 'Normalization', 'probability','facecolor',nice_cmap(20,:),'facealpha',.5,'edgecolor',nice_cmap(10,:));
% leg1 = sprintf('Before filtering, %d genes', length(rall)); 
% leg2 = sprintf('After filtering, %d genes', length(r)); 
% 
% legend(leg1, leg2)
% xlabel('Correlation between probes for the same gene'); ylabel('Probability'); 
% set(gcf,'color','w');
% 
% save('IDgenes5plus.mat', 'IDgene'); 

uniqNrProbes = unique(numberProbes);
HowManyProbes = zeros(length(uniqNrProbes),1);
HowManyGenes = zeros(length(uniqNrProbes),1);

for i=1:length(uniqNrProbes)
    
    HowManyProbes(i) = uniqNrProbes(i);
    HowManyGenes(i) = length(find(numberProbes==uniqNrProbes(i)));
    
end
probeSummary = table(HowManyGenes,HowManyProbes);
figure; 
scatter(probeSummary.HowManyProbes, probeSummary.HowManyGenes, 'filled'); ...
    hold on; plot(probeSummary.HowManyProbes, probeSummary.HowManyGenes);
xlabel('Number of probes'); ylabel('Number of genes'); 

%figure; hist(r1,100); figure; hist(r2,100);  figure; hist(r3,100);
cd ../../..
% check the proportion of noise



