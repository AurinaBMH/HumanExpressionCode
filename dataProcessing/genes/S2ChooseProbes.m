%% Author: Aurina
%% Date modified: 2016-03-10
%% Date modified: 2017-07-14
%% This script:
%   1. Loads all microarray data MicroarrayData.mat file
%   2. Uses 1 of 2 possible probe selection options: by max variance or by highest PC
%   3. Finds probes with duplicate entrezIDs
%   4. Chooses one probe from them according to the selected option
%   5. Removes exluded probes from expression data and all other fields.
%   6. Saves MicroarrayDataMaxVAR.mat'or MicroarrayDataPCA.mat with Expressionall, Probe information, Sample information.
%%
%choose one of those
% To use probes with max variance choose UseMaxVarProbes = 1; keeping the other option = 0;
% To use probes with highest PC choose UseHighestPCProbes = 1; keeping the other option = 0;
% choose if you want to use data with CUST probes or without them:
% UseDataWithCUSTprobes = 1; if UseDataWithCUSTprobes=0, it will load data
% without cust probes.

UseDataWithCUSTprobes = false;
probeSelection = 'Variance';% (Variance', LessNoise', 'Mean', 'PC')

%------------------------------------------------------------------------------
% Load the data
%------------------------------------------------------------------------------
cd ('data/genes/processedData');
if UseDataWithCUSTprobes
    fprintf(1,'Loading the data with CUST probes and assigning variables\n')
    load('MicroarrayDataWITHCUST.mat');
else
    fprintf(1,'Loading the data without CUST probes and assigning variables\n')
    load('MicroarrayData.mat');
end
cd ..
cd ('rawData');

ProbeName = DataTableProbe.ProbeName{1,1};
ProbeID = DataTableProbe.ProbeID{1,1};
EntrezID = DataTableProbe.EntrezID{1,1};
GeneID = DataTableProbe.GeneID{1,1};
GeneSymbol = DataTableProbe.GeneSymbol{1,1};
GeneName = DataTableProbe.GeneName{1,1};
%------------------------------------------------------------------------------
% Normalise each subject data separately
% Calculate probe selection criteria for each subject separately (non
% normalised data)
% Choose probe with max probe selection criteria on average
%------------------------------------------------------------------------------
% sigmoid normalise data for MaxVariance calculation
% ExpressionallNorm = BF_NormalizeMatrix(Expressionall');
% ExpressionallNorm = ExpressionallNorm';




%% find repeating entrezIDs and calculate variances for them, then take a probe of max variance.
%
% % ------------------------------------------------------------------------------
% % Find best representative in a set of duplicates using maxVar and remove all others:
% % ------------------------------------------------------------------------------
Uniq = unique(EntrezID(:,1));
N = histc(EntrezID, Uniq);
ProbeList = zeros(length(Uniq),2);
indMsubj = zeros(length(Uniq),6);
expressionSelected = cell(6,1);

fileNoise = 'PACall.csv';

for subj = 1:6
    expression = (DataTable.Expression{subj})';
    if strcmp (probeSelection, 'LessNoise')
        folder = sprintf('normalized_microarray_donor0%d', subj);
        cd (folder);
        noise = csvread(fileNoise);
        [~,probeList] = intersect(noise(:,1),ProbeID, 'stable');
        noise = (noise(probeList,2:end))';
        cd ..
    end
    % load noise level matrix for each subject here
    
    for k=1:length(Uniq)
        fprintf(1,'Processing entrez ID %u\n',Uniq(k))
        % find indexes for repeating entrexIDs
        indRepEntrezIDs = find(EntrezID==(Uniq(k)));
        if length(indRepEntrezIDs) >=2
            fprintf(1,'%d duplicates found\n', length(length(indRepEntrezIDs)));
            
            % take expression values for a selected entrezID
            expRepEntrezIDs = expression(:,indRepEntrezIDs);
            % calculate variances for expression data for a selected entrezID
            switch probeSelection
                
                case 'Variance'
                    fprintf(1,'Performing probe selection using Max variance\n');
                    measure = var(expRepEntrezIDs,0,1);
                    % determine max var value
                    [MaxV, indMaxV] = max(measure);
                    
                case 'PC'
                    fprintf(1,'Performing probe selection using max PC\n');
                    measure = pca(expRepEntrezIDs,'Centered',false);
                    % determine max PC loading
                    [MaxV, indMaxV] = max(measure(:,1));
                    
                case 'LessNoise'
                    fprintf(1,'Performing probe selection using less noise criteria\n');
                    noiseRepEntrezIDs = noise(:,indRepEntrezIDs);
                    % find probe with most signal in it compared to noise
                    measure = sum(noiseRepEntrezIDs,1);
                    % determine probe with max signal
                    [MaxV, indMaxV] = max(measure);
                    
                case 'Mean'
                    expressionSelected{subj}(:,k) = mean(expRepEntrezIDs,2);
                    
            end
        end
        
        if strcmp(probeSelection, 'Variance') || strcmp(probeSelection, 'PC') || strcmp(probeSelection, 'LessNoise')
            if length(indRepEntrezIDs) >=2
                indMsubj(k,subj) = indRepEntrezIDs(indMaxV);
                
                
            else
                indMsubj(k,subj) = indRepEntrezIDs;
            end
        end
    end
end

if strcmp(probeSelection, 'Variance') || strcmp(probeSelection, 'PC') || strcmp(probeSelection, 'LessNoise')
    indProbe = zeros(length(Uniq),1);
    for j=1:length(Uniq)
        [hcount,index] = hist(indMsubj(j,:),unique(indMsubj(j,:)));
        [~, ind] = max(hcount);
        indINlist = index(ind);
        ProbeList(j,1) = ProbeID(indINlist);
        ProbeList(j,2) = EntrezID(indINlist);
    end
    %% check to exclude probes with not maximum variance or max PC and repeating entrezIDs (assign NaN value to
    % them)
    for q=1:length(ProbeID)
        if ProbeID(q)~=(ProbeList(:,1))
            ProbeID(q) = NaN;
            EntrezID(q) = NaN;
        end
    end
    
    EntrezID(isnan(EntrezID)) = [];
    ProbeInformation.EntrezID = EntrezID;
    
    GeneID(isnan(ProbeID)) = [];
    ProbeInformation.GeneID = GeneID;
    
    GeneSymbol(isnan(ProbeID)) = [];
    ProbeInformation.GeneSymbol = GeneSymbol;
    
    GeneName(isnan(ProbeID)) = [];
    ProbeInformation.GeneName = GeneName;
    
    ProbeName(isnan(ProbeID)) = [];
    ProbeInformation.ProbeName = ProbeName;
    
    RemoveProbes = ProbeID;
    
    ProbeID(isnan(ProbeID)) = [];
    ProbeInformation.ProbeID = ProbeID;
    
    for subject=1:6
        
        
        % exclude NaN probes keeping 1 probe for 1 entrezID.
        fprintf(1,'Combining and saving the data for subject %u\n', subject)
        Expression = DataTable.Expression{subject,1};
        Expression(isnan(RemoveProbes),:) = [];
        Expression = Expression';
        
        
        % combine sample information variables to a structure.
        SampleInformation.StructureNames = DataTable.StructureName{subject,1};
        SampleInformation.MMCoordinates = DataTable.MMcoordinates{subject,1};
        SampleInformation.MRIvoxCoordinates = DataTable.MRIvoxCoordinates{subject,1};
        cd ..
        cd ('processedData');
        if ~UseDataWithCUSTprobes
            save(sprintf('MicroarrayData%sS0%d.mat', probeSelection, subject), 'Expression', 'ProbeInformation' , 'SampleInformation');
        else
            save(sprintf('MicroarrayDataWITHcust%sS0%d.mat', probeSelection, subject), 'Expression', 'ProbeInformation' , 'SampleInformation');
        end
        
    end
    
else
    
    for subject=1:6
        
        % exclude NaN probes keeping 1 probe for 1 entrezID.
        fprintf(1,'Combining and saving the data for subject %u\n', subject)
        Expression = expressionSelected{subject};
        
        % combine sample information variables to a structure.
        SampleInformation.StructureNames = DataTable.StructureName{subject,1};
        SampleInformation.MMCoordinates = DataTable.MMcoordinates{subject,1};
        SampleInformation.MRIvoxCoordinates = DataTable.MRIvoxCoordinates{subject,1};
        
        
        ProbeInformation.EntrezID = unique(EntrezID, 'stable');
        ProbeInformation.GeneID = unique(GeneID, 'stable');
        ProbeInformation.GeneSymbol = unique(GeneSymbol, 'stable');
        ProbeInformation.GeneName = unique(GeneName, 'stable');
        ProbeInformation.ProbeName = unique(ProbeName, 'stable');
        cd ..
        cd ('processedData');
        if ~UseDataWithCUSTprobes
            save(sprintf('MicroarrayData%sS0%d.mat', probeSelection, subject), 'Expression', 'ProbeInformation' , 'SampleInformation');
        else
            save(sprintf('MicroarrayDataWITHcust%sS0%d.mat', probeSelection, subject), 'Expression', 'ProbeInformation' , 'SampleInformation');
        end
        
    end
    
end




