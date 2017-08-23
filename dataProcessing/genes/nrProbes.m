% How many genes have multiple probes
%------------------------------------------------------------------------------
% Choose options
%------------------------------------------------------------------------------

useCUSTprobes = true;
probeSelection = 'PC';% (Variance', LessNoise', 'Mean', 'PC')
signalThreshold = 0.5; % percentage of samples that a selected probe has expression levels that are higher than background
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

cd ..
cd ('rawData');

% % ------------------------------------------------------------------------------
% % Find best representative in a set of duplicates using maxVar and remove all others:
% % ------------------------------------------------------------------------------
expressionSelected = cell(6,1);
noiseSUBJ = cell(6,1);
fileNoise = 'PACall.csv';

ProbeID = DataTableProbe.ProbeID{1,1};
% % ------------------------------------------------------------------------------
% % First, find probes that have very noisy data and remove them from consideration
% % Threshold for removing those proges is defined as the percentage of
% % samples a probe has expression higher that background
% % ------------------------------------------------------------------------------
for subject = 1:6
    folder = sprintf('normalized_microarray_donor0%d', subject);
    cd (folder);
    noise2filter = csvread(fileNoise);
    [~,probeList] = intersect(noise2filter(:,1),ProbeID, 'stable');
    noise2filter = (noise2filter(probeList,2:end))';
    noiseSUBJ{subject} = noise2filter;
    cd ..
end
% combine noise data for all subjects
noiseALL = vertcat(noiseSUBJ{1}, noiseSUBJ{2}, noiseSUBJ{3}, noiseSUBJ{4}, noiseSUBJ{5}, noiseSUBJ{6});
% calculate the percentage of samples that each probe has expression value
% higher than a selected number
signalLevel = sum(noiseALL,1)./size(noiseALL,1);
indKeepProbes = find(signalLevel>signalThreshold);


cd ../../..

parcellation = 'aparcaseg';%, 'cust100', 'cust250'};
distanceThreshold = 2; % first run 30, then with the final threshold 2
normMethod = 'zscore';
normaliseWhat = 'Lcortex'; %(LcortexSubcortex, wholeBrain, LRcortex)

% choose Lcortex if want to normalise samples assigned to left cortex separately;
% choose LcortexSubcortex if want to normalise LEFT cortex + left subcortex together
% choose wholeBrain if you want to normalise the whole brain.
% choose LRcortex if you want to normalise left cortex + right cortex.
%------------------------------------------------------------------------------
% Assign variables
%------------------------------------------------------------------------------
if strcmp(parcellation, 'aparcaseg')
    NumNodes = 82;
elseif strcmp(parcellation, 'cust100')
    NumNodes = 220;
elseif strcmp(parcellation, 'cust250')
    NumNodes = 530;
end
%------------------------------------------------------------------------------
% ORIGINAL DATA OR FILTERED DATA
%------------------------------------------------------------------------------
doOriginal = false;

if doOriginal
    genes = unique(DataTableProbe.EntrezID{1},'stable');
    listGenes = DataTableProbe.EntrezID{1};
else
    genes = unique(DataTableProbe.EntrezID{1}(indKeepProbes),'stable');
    listGenes = DataTableProbe.EntrezID{1}(indKeepProbes);
    signalLevel = signalLevel(indKeepProbes);
end

signalLevelProbes = cell(length(genes),1);
numberProbes = zeros(length(genes),1);
for gene=1:length(genes)
    
    indGene = find(listGenes==genes(gene));
    signalLevelProbes{gene} = signalLevel(indGene);
    numberProbes(gene) = length(indGene);
    
end

uniqNrProbes = unique(numberProbes);
HowManyProbes = zeros(length(uniqNrProbes),1);
HowManyGenes = zeros(length(uniqNrProbes),1);
for i=1:length(uniqNrProbes)
    
    HowManyProbes(i) = uniqNrProbes(i);
    HowManyGenes(i) = length(find(numberProbes==uniqNrProbes(i)));
    
end
probeSummary = table(HowManyGenes,HowManyProbes);

% check the proportion of noise



