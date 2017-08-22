%% Author: Aurina
%% Date modified: 2016-03-10
%% Date modified: 2017-07-14
%% Date modified: 2017-08-02
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

% remove selected probes from data and perform other calculations only on
% non-noisy probes
ProbeName = DataTableProbe.ProbeName{1,1}(indKeepProbes);
ProbeID = ProbeID(indKeepProbes);
EntrezID = DataTableProbe.EntrezID{1,1}(indKeepProbes);
GeneID = DataTableProbe.GeneID{1,1}(indKeepProbes);
GeneSymbol = DataTableProbe.GeneSymbol{1,1}(indKeepProbes);
GeneName = DataTableProbe.GeneName{1,1}(indKeepProbes);

Uniq = unique(EntrezID(:,1));
N = histc(EntrezID, Uniq);
ProbeList = zeros(length(Uniq),2);
indMsubj = zeros(length(Uniq),6);

for subj = 1:6
    expression = (DataTable.Expression{subj}(indKeepProbes,:))';
    %if strcmp (probeSelection, 'LessNoise')
    %folder = sprintf('normalized_microarray_donor0%d', subj);
    %cd (folder);
    noise = noiseSUBJ{subj}; %csvread(fileNoise);
    %[~,probeList] = intersect(noise(:,1),ProbeID, 'stable');
    noise = (noise(:,indKeepProbes));
    %noiseALL{subj} = noise;
    %cd ..
    %end
    % load noise level matrix for each subject here
    
    for k=1:length(Uniq)
        fprintf(1,'Processing entrez ID %u\n',Uniq(k))
        % find indexes for repeating entrexIDs
        indRepEntrezIDs = find(EntrezID==(Uniq(k)));
        expRepEntrezIDs = expression(:,indRepEntrezIDs);
        if length(indRepEntrezIDs) >=2
            fprintf(1,'%d duplicates found\n', length(length(indRepEntrezIDs)));
            
            % take expression values for a selected entrezID
            
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
        else
            switch probeSelection
                case 'Mean'
                    expressionSelected{subj}(:,k) = expRepEntrezIDs;
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
        Expression = DataTable.Expression{subject,1}(indKeepProbes,:);
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
        save(sprintf('%s%sS0%d.mat', startFileName, probeSelection, subject), 'Expression', 'ProbeInformation' , 'SampleInformation');
        
        
    end
    
end




