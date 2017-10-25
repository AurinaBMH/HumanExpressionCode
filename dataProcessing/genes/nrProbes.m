% How many genes have multiple probes
%------------------------------------------------------------------------------
% Choose options
%------------------------------------------------------------------------------

useCUSTprobes = true;
signalThreshold = 0.5; % percentage of samples that a selected probe has expression levels that are higher than background
doOriginal = false;
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
indKeepProbes = find(signalLevel>signalThreshold);

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
        [r(k),p(k)] = corr(expression(indGene(1),:)', expression(indGene(2),:)'); 
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
    if numberProbes(gene)>1
        IDgene(m) = listGenes(indGene(1));
         m=m+1; 
    end
   
    
end

save('IDgenes2plus.mat', 'IDgene'); 

uniqNrProbes = unique(numberProbes);
HowManyProbes = zeros(length(uniqNrProbes),1);
HowManyGenes = zeros(length(uniqNrProbes),1);

for i=1:length(uniqNrProbes)
    
    HowManyProbes(i) = uniqNrProbes(i);
    HowManyGenes(i) = length(find(numberProbes==uniqNrProbes(i)));
    
end
probeSummary = table(HowManyGenes,HowManyProbes);
figure; hist(r1,100); figure; hist(r2,100);  figure; hist(r3,100);
cd ../../..
% check the proportion of noise



