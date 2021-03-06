%% Author: Aurina
%% Date modified: 2016-03-10
%% Date modified: 2017-07-14
%% Date modified: 2017-08-02
%% This script:
%   1. Loads all microarray data MicroarrayData.mat file
%   2. Uses 1 of 5 possible probe selection options: by max variance or by
%   highest PC, less noise, mean of all probes and random.
%   3. Finds probes with duplicate entrezIDs
%   4. Chooses one probe from them according to the selected option
%   5. Removes exluded probes from expression data and all other fields.
%   6. Saves MicroarrayDataMaxVAR.mat'or MicroarrayDataPCA.mat with Expressionall, Probe information, Sample information.
%%
% choose if you want to use data with CUST probes or without them:
% UseDataWithCUSTprobes = 1; if UseDataWithCUSTprobes=0, it will load data
% without cust probes.
clear all;
useCUSTprobes = false;
probeSelection = 'RNAseq'; %{'Mean', 'Variance', 'LessNoise', 'Random', 'PC', 'RNAseq', DS, CV};% probeSelection = {'Mean', 'Variance', 'LessNoise', 'Random', 'PC'};
RNAseqThreshold = 0.2;
signalThreshold = 0.5;% percentage of samples that a selected probe has expression levels that are higher than background
rng shuffle % for selecting different seed for random probe selection
%------------------------------------------------------------------------------
% Load the data
%------------------------------------------------------------------------------
cd ('data/genes/processedData');

if useCUSTprobes
    fprintf(1,'Loading the data with CUST probes and assigning variables\n')
    
    startFileName = 'MicroarrayDataWITHcust';
else
    fprintf(1,'Loading the data without CUST probes and assigning variables\n')
    startFileName = 'MicroarrayDataProbesUpdated';
end
fprintf(1,sprintf('Probe selection based on %s is chosen\n', probeSelection))
load(sprintf('%s.mat', startFileName));

cd ..
cd ('rawData');

%------------------------------------------------------------------------------
% Calculate probe selection criteria for each subject separately (non
% normalised data)
% Choose probe with max probe selection criteria on average
%------------------------------------------------------------------------------

%% find repeating entrezIDs and calculate variances for them, then take a probe of max variance.
%------------------------------------------------------------------------------
% Find best representative in a set of duplicates using maxVar and remove all others:
%------------------------------------------------------------------------------
expressionSelected = cell(6,1);
noiseSUBJ = cell(6,1);

ProbeID = DataTableProbe.ProbeID{1,1};
% % ------------------------------------------------------------------------------
% First, find probes that have very noisy data and remove them from consideration
% Threshold for removing those proges is defined as the percentage of
% samples a probe has expression higher that background
% % ------------------------------------------------------------------------------

signalLevel = sum(noiseall,2)./size(noiseall,2);
indKeepProbes = find(signalLevel>=signalThreshold);

% remove selected probes from data and perform other calculations only on
% non-noisy probes
ProbeName = DataTableProbe.ProbeName{1,1}(indKeepProbes);
ProbeID = ProbeID(indKeepProbes);
EntrezID = DataTableProbe.EntrezID{1,1}(indKeepProbes);
%GeneID = DataTableProbe.GeneID{1,1}(indKeepProbes);
GeneSymbol = DataTableProbe.GeneSymbol{1,1}(indKeepProbes);
%GeneName = DataTableProbe.GeneName{1,1}(indKeepProbes);

probeInformationALL.ProbeName = ProbeName;
probeInformationALL.ProbeID = ProbeID;
probeInformationALL.EntrezID = EntrezID;
probeInformationALL.GeneSymbol = GeneSymbol;


% if choosing probes based on RNAseq, then use data only from 1 subject
if strcmp(probeSelection, 'RNAseq')
    
    [correlations, avgCorr, indProbe, genes] = selectProbeRNAseq(DataTable, EntrezID, indKeepProbes, RNAseqThreshold);
    nSub = 1;
elseif strcmp(probeSelection, 'DS')
    [indProbe, avCorr] = selectProbeDS(EntrezID, DataTable, indKeepProbes);
    nSub = 1;
else
    nSub = 6;
end


Uniq = unique(EntrezID);
ProbeList = zeros(length(Uniq),2);
indMsubj = zeros(length(Uniq),nSub);

for subj = 1:nSub
    expression = (DataTable.Expression{subj}(indKeepProbes,:))';
    noise = DataTable.Noise{subj}(indKeepProbes,:)';
    
    for k=1:length(Uniq)
        %fprintf(1,'Processing entrez ID %u\n',Uniq(k))
        % find indexes for repeating entrexIDs
        indRepEntrezIDs = find(EntrezID==(Uniq(k)));
        expRepEntrezIDs = expression(:,indRepEntrezIDs);
        if length(indRepEntrezIDs) >=2
            
            %fprintf(1,'%d duplicates found\n', length(indRepEntrezIDs));
            
            % take expression values for a selected entrezID
            % calculate variances for expression data for a selected entrezID
            %switch probeSelection
            if strcmp(probeSelection, 'Variance')
                %case 'Variance'
                %fprintf(1,'Performing probe selection using Max variance\n');
                measure = var(expRepEntrezIDs,0,1);
                % determine max var value
                [MaxV, indMaxV] = max(measure);
                
                
            elseif strcmp(probeSelection, 'PC')
                %fprintf(1,'Performing probe selection using max PC\n');
                % substract the mean before doing pca
                expRepEntrezIDsNOmean = expRepEntrezIDs-mean(expRepEntrezIDs); 
                measure = pca(expRepEntrezIDsNOmean,'Centered',false);
                % determine max PC loading
                [MaxV, indMaxV] = max(measure(:,1));
                
                
            elseif strcmp(probeSelection, 'CV')
                %fprintf(1,'Performing probe selection using max PC\n');
                measure = std(expRepEntrezIDs)./mean(expRepEntrezIDs);
                % determine max PC loading
                [MaxV, indMaxV] = max(measure);
                
            elseif strcmp(probeSelection, 'LessNoise')
                %fprintf(1,'Performing probe selection using less noise criteria\n');
                noiseRepEntrezIDs = noise(:,indRepEntrezIDs);
                % find probe with most signal in it compared to noise
                measure = sum(noiseRepEntrezIDs,1);
                % determine probe with max signal
                [MaxV, indMaxV] = max(measure);
                
            elseif strcmp(probeSelection, 'Mean')
                expressionSelected{subj}(:,k) = mean(expRepEntrezIDs,2);
                [indMaxV] = randsample(1:size(expRepEntrezIDs,2),1);
                
                % choose one of probes as a place holder (at random)
                % for probeInformation (expression values are mean of
                % al probes)
                
            elseif strcmp(probeSelection, 'Random')
                %fprintf(1,'Performing probe selection using random selection\n');
                
                % determine max var value
                [indMaxV] = randsample(1:size(expRepEntrezIDs,2),1);
            elseif strcmp(probeSelection, 'RNAseq')
                indMaxV = indProbe(k);
            elseif strcmp(probeSelection, 'DS')
                %indMaxV = indProbe(k);
                indMsubj(k,subj) = indProbe(k);
                
            end
            
            if (strcmp(probeSelection, 'Mean') || strcmp(probeSelection, 'Variance') || strcmp(probeSelection, 'Random') ||...
                    strcmp(probeSelection, 'PC') || strcmp(probeSelection, 'LessNoise') || strcmp(probeSelection, 'RNAseq') || strcmp(probeSelection, 'CV'))
                %indMsubj(k,subj) = indProbe(k); %(indMaxV);
                %if NaN, use NaN;
                if isnan(indMaxV)
                    indMsubj(k,subj) = NaN;
                    
                else
                    indMsubj(k,subj) = indRepEntrezIDs(indMaxV);
                end
            end
            
        else
            if strcmp(probeSelection, 'Mean')
                expressionSelected{subj}(:,k) = expRepEntrezIDs;
                indMsubj(k,subj) = indRepEntrezIDs;
                
            elseif strcmp(probeSelection, 'RNAseq')
                if isnan(indProbe(k))
                    indMsubj(k,subj) = NaN;
                else
                    indMsubj(k,subj) = indRepEntrezIDs;
                end
            else
                indMsubj(k,subj) = indRepEntrezIDs;
            end
            
        end
        
    end
end

for j=1:length(Uniq)
    
    indINlist = mode(indMsubj(j,:),2);
    if isnan(indINlist)
        ProbeList(j,1) = NaN;
        ProbeList(j,2) = NaN;
    else
        ProbeList(j,1) = ProbeID(indINlist);
        ProbeList(j,2) = EntrezID(indINlist);
    end
    
end
% %% check to exclude probes with not maximum variance or max PC and repeating entrezIDs (assign NaN value to
% % them)
[reordered, reorder_ind] = sort(EntrezID);
% reorder all values based on sorted entrezIDs, because this is the order
% of selected probes (1...n) from Unique.
EntrezID = EntrezID(reorder_ind);
ProbeID = ProbeID(reorder_ind);
%GeneID = GeneID(reorder_ind);
GeneSymbol = GeneSymbol(reorder_ind);
%GeneName = GeneName(reorder_ind);
ProbeName = ProbeName(reorder_ind);

[a,ind2rem] = setdiff(ProbeID, ProbeList(:,1));
ProbeID(ind2rem) = NaN;
EntrezID(ind2rem) = NaN;

EntrezID(isnan(EntrezID)) = [];
%GeneID(isnan(ProbeID)) = [];
GeneSymbol(isnan(ProbeID)) = [];
%GeneName(isnan(ProbeID)) = [];
ProbeName(isnan(ProbeID)) = [];




% % ------------------------------------------------------------------------------
% Check if all genes that are left have unique gene symbols
% % ------------------------------------------------------------------------------

% [~, ind] = unique(GeneSymbol, 'stable');
% % find duplicate indices
% duplicate_ind = setdiff(1:size(GeneSymbol,1), ind);
%
% toExclude = cell(length(duplicate_ind),1);
% for gene = 1:length(duplicate_ind)
%
%     % if gene names are not matching, exclude both of them
%     symbs = GeneSymbol(duplicate_ind(gene));
%     test_ind = find(strcmp(GeneSymbol, symbs{1}));
%     % check if for all duplicated instances gene name is the same
%
%     names = GeneName(test_ind);
%     doMatch = zeros(length(test_ind));
%     for j=1:length(test_ind)
%         for i=j+1:length(test_ind)
%
%             doMatch(j,i) = strcmp(names{j}, names{i});
%
%         end
%     end
%     % if they're not matching, record indexes to exclude later
%     if sum(doMatch)==0
%         toExclude{gene,:} = test_ind;
%     end
%
% end
excludeGenes = []; %unique(cell2mat(toExclude));
% make a vector if indexes for genes to keep
keepGenes = setdiff(1:size(GeneSymbol,1), excludeGenes, 'stable');% changed from stable

probeInformation.EntrezID = EntrezID(keepGenes); %[reordered, reorder_ind] = sort(probeInformation.EntrezID);
%probeInformation.GeneID = GeneID(keepGenes);
probeInformation.GeneSymbol = GeneSymbol(keepGenes);
%probeInformation.GeneName = GeneName(keepGenes);
probeInformation.ProbeName = ProbeName(keepGenes);
ProbeID2nd = ProbeID; % assign probe ID values to a variable (will be used to filter gene expression values)
ProbeID(isnan(ProbeID)) = []; % remove redundant probes
probeInformation.ProbeID = ProbeID(keepGenes);

cd ..
cd ('processedData')
expressionAll = cell(6,1);
sampleInfo = cell(6,1);
for subject=1:6
    if strcmp(probeSelection, 'Variance') || strcmp(probeSelection, 'PC') ...
            || strcmp(probeSelection, 'LessNoise') || strcmp(probeSelection, 'Random') ...
            || strcmp(probeSelection, 'RNAseq') || strcmp(probeSelection, 'DS') || strcmp(probeSelection, 'CV')
        
        % exclude NaN probes keeping 1 probe for 1 entrezID.
        fprintf(1,'Combining and saving the data for subject %u\n', subject)
        Expression = DataTable.Expression{subject,1}(indKeepProbes,:); % filter noisy probes
        Expression = Expression(reorder_ind,:); % reorder probes according to sorted entrezIDs
        Expression(isnan(ProbeID2nd),:) = []; % exclude probes that were not selected
        Expression = Expression(keepGenes,:); % exclude genes that had duplicated and non-matching names
        expressionAll{subject} = Expression';
        
    elseif strcmp(probeSelection, 'Mean')
        
        % exclude NaN probes keeping 1 probe for 1 entrezID.
        fprintf(1,'Combining and saving the data for subject %u\n', subject)
        expressionAll{subject} = expressionSelected{subject}(:,keepGenes); % exclude genes that had duplicated and non-matching names
    end
    % combine sample information variables to a structure.
    SampleInformation.StructureNames = DataTable.StructureName{subject,1};
    SampleInformation.MMCoordinates = DataTable.MMcoordinates{subject,1};
    SampleInformation.MRIvoxCoordinates = DataTable.MRIvoxCoordinates{subject,1};
    sampleInfo{subject} = SampleInformation;
    
end

if strcmp(probeSelection, 'RNAseq')
    save(sprintf('%s%s%dRNAthr%dnoisethr.mat', startFileName, probeSelection, RNAseqThreshold, signalThreshold), 'expressionAll', 'probeInformation' , 'sampleInfo', 'avgCorr', 'probeInformationALL', 'genes');
else
    save(sprintf('%s%s.mat', startFileName, probeSelection), 'expressionAll', 'probeInformation' , 'sampleInfo');
end
cd ../../..
