%% Author: Aurina
%% Date modified: 2016-03-22
%% Date modified: 2017-07-14
%% Date modified: 2017-08-29 - noise level data added
%% Date modified: 2017-11-10 - manual gene symbol naming fixed
%% Date modified: 2018-01-10 - option to update probe--> gene assingment added
%% This script:
%   1. Loads all microarray data from excell files for each subject
%   2. Excludes custom probes;
%   3. Excludes probes with missing entrezIDs or updates probe to gene
%   assignment (depending on options chosen)
%   4. Saves expression data, coordinates, sample structure names for all samples
%   5. Saves data for separate subjects to DataTable
%   6. Saves data for all subjects combined as variables 'MicorarrayData.mat' file
%%
%------------------------------------------------------------------------------
% Choose options
%------------------------------------------------------------------------------
clear all; 
ExcludeCBandBS = true; % will exclude samples from brainstem and cerebellum
useCUSTprobes = false; % use CUST probes
updateProbes = true; % update probe --> gene assignment based on latest data

if ExcludeCBandBS
    fprintf(1,'Brainstem and cerebellum samples will be EXCLUDED\n')
else
    fprintf(1,'Brainstem and cerebellum samples will be INCLUDED\n')
end

if updateProbes
    fprintf(1,'Probe to gene assignment will be UPDATED\n')
else
    fprintf(1,'Probe to gene assignment will NOT be UPDATED\n')
end

if useCUSTprobes
    fprintf(1,'CUST probes will be INCLUDED\n')
    startFileName = 'MicroarrayDataWITHcust';
else
    fprintf(1,'CUST probes will be EXCLUDED\n')
    startFileName = 'MicroarrayData';
end

cd ('data/genes/rawData');
%% load probe information (same for all subjects)
fprintf(1,'Loading Probes.xlsx file\n')
FileProbes = 'Probes.xlsx';
ProbeTable = readtable(FileProbes);
ProbeID = ProbeTable.probe_id;
EntrezID = ProbeTable.entrez_id;
ProbeName =  ProbeTable.probe_name;
GeneID = ProbeTable.gene_id;
GeneSymbol = ProbeTable.gene_symbol;
GeneName = ProbeTable.gene_name;

% do not remove probes with missing
% entrez IDs because some of them can be updated
fprintf(1,'%d probes with missing entrez IDs\n', sum(isnan(EntrezID)));
% creat a Data cell to store the output
headerdata = {'Expression' , 'MMcoordinates', 'StructureName', 'MRIvoxCoordinates', 'Noise', 'SampleID', 'WellID'};
headerprobe = { 'ProbeID', 'EntrezID','ProbeName', 'GeneSymbol'};
Data = cell(6,7);
DataProbe = cell(1,4);

%------------------------------------------------------------------------------
% Go to each subject's directory and take the data
%------------------------------------------------------------------------------

for subj=1:6
    fprintf(1,'Loading data for %u subject\n', subj)
    folder = sprintf('normalized_microarray_donor0%d', subj);
    cd (folder);
    %%load information specific for each subject
    FileMicroarray = 'MicroarrayExpression.csv';
    FileAnnot = 'SampleAnnot.xlsx';
    FileNoise = 'PACall.csv';
    Expression = csvread(FileMicroarray);
    noise = csvread(FileNoise);
    
    Expression(:,1) = [];                         % exclude probe IDs from expression matrix
    [~,~,SlabType] = xlsread(FileAnnot, 'D:D');
    [~,~, StructureName] = xlsread(FileAnnot, 'F:F');
    SlabType(1) = [];                           % remove headline
    StructureName(1) = [];                      % remove headline
    MMcoordinates = xlsread(FileAnnot, 'K:M');
    MRIvoxCoordinates = xlsread(FileAnnot, 'H:J');
    SampleID = xlsread(FileAnnot, 'A:A');
    WellID = xlsread(FileAnnot, 'C:C');
    [~,probeList] = intersect(noise(:,1),ProbeID, 'stable');
    noise = noise(probeList,2:end);
    
    % To exclude expression data and coordinates for braintem (BS) and cerebellum (CB)
    % exclude columns in expression and rows in coordinates if slabtype is CB or BS
    if ExcludeCBandBS
        fprintf('Excluding brainstem and cerebellum data\n')
        
        BS = strfind(SlabType, 'BS'); BSind = find(~cellfun(@isempty,BS));
        CB = strfind(SlabType, 'CB'); CBind = find(~cellfun(@isempty,CB));
        BSandCBind = [BSind;CBind];
        
        fprintf(1,'%d cerebellum and brainstem samples to remove\n', length(BSandCBind))
        
        Expression(:,BSandCBind) = NaN;
        MMcoordinates(BSandCBind,:) = NaN;
        MRIvoxCoordinates(BSandCBind,:) = NaN;
        SampleID(BSandCBind,:) = NaN;
        WellID(BSandCBind,:) = NaN; 
        noise(:,BSandCBind) = NaN;
        StructureName(BSandCBind) = {NaN};
    end
    
    % for nan columns
    % keep only existing expression values
    Expression = Expression(:,all(~isnan(Expression)));
    % keep only existing coordinates
    MMcoordinates = MMcoordinates(all(~isnan(MMcoordinates),2),:); % for nan rows
    MRIvoxCoordinates = MRIvoxCoordinates(all(~isnan(MRIvoxCoordinates),2),:); % for nan rows
    SampleID = SampleID(all(~isnan(SampleID),2),:); % for nan rows
    WellID = WellID(all(~isnan(WellID),2),:); % for nan rows
    noise = noise(:,all(~isnan(noise)));
    % keep only existing structure names
    StructureName(cellfun(@(StructureName) any(isnan(StructureName)),StructureName)) = [];
    
    % assign output to Data cell;
    Data{subj,1} = Expression;
    Data{subj,2} = MMcoordinates;
    Data{subj,3} = StructureName;
    Data{subj,4} = MRIvoxCoordinates;
    Data{subj,5} = noise;
    Data{subj,6} = SampleID; 
    Data{subj,7} = WellID; 
    cd ..
    
end
%
% %------------------------------------------------------------------------------
% % Make a table from all the data
% %------------------------------------------------------------------------------
DataTable = dataset({Data, headerdata{:}});

if ~useCUSTprobes
    %------------------------------------------------------------------------------
    % Remove CUST probes:
    %------------------------------------------------------------------------------
    % Remove all CUST probes (assign NaN values for all custom probes)
    fprintf(1,'Removing CUST probes\n')
    cust = strfind(ProbeName, 'CUST');
    remInd = find(~cellfun(@isempty,cust));
    fprintf(1,'%d CUST probes removed\n', length(remInd))
    ProbeName(remInd) = {NaN};
    ProbeID(remInd) = NaN;
    
ProbeName(isnan(ProbeID)) = [];
EntrezID(isnan(ProbeID)) = [];
GeneID(isnan(ProbeID)) = [];
GeneSymbol(isnan(ProbeID)) = [];
GeneName(isnan(ProbeID)) = [];

    for s=1:6
       DataTable.Expression{s,1}(isnan(ProbeID),:) = []; 
       DataTable.Noise{s,1}(isnan(ProbeID),:) = [];
    end

ProbeID(isnan(ProbeID)) = [];
    
end

if updateProbes
    %------------------------------------------------------------------------------
    % Update probe to gene annotation based on latest data
    %------------------------------------------------------------------------------
    % load re-annotated probes
    updatedProbes = importProbes('mart_export_updatedProbes.txt');
    
    % find probes that are mentioned more than once
    [~, i, j] = unique(updatedProbes.probeNames,'first');
    indexToDupes = find(not(ismember(1:numel(updatedProbes.probeNames),i)));
    
    % exclude those probes
    updatedProbes(indexToDupes,:) = [];
    % now each probe is mapped to one gene
    % It's time to compare allen probes with re-annotated probes.
    allenProbes = table;
    allenProbes.probeNames = ProbeName;
    allenProbes.geneName = GeneSymbol;
    allenProbes.NCBIgeneID = EntrezID;
    
    updatedMatching = comparehg38VSAllen(updatedProbes, allenProbes);
    
    % sumarise the information in numbers
    nm = length(find(updatedMatching.compare==1));
    fprintf(1,'%d probes are matching\n', nm)
    nmm = length(find(updatedMatching.compare==0));
    fprintf(1,'%d probes are mismatching\n', nmm)
    nu = length(find(updatedMatching.compare==2));
    fprintf(1,'%d probes are introduced with IDs\n', nu)
    
    % find probes that are in both lists and replace geneSymbol and entrezIDs
    % with updated.
    [probesSelect, INDold, INDnew] = intersect(ProbeName, updatedMatching.probeNames);
    
    ProbeName = ProbeName(INDold);
    ProbeID = ProbeID(INDold);
    EntrezID = updatedMatching.NCBIgeneID(INDnew);
    GeneSymbol = updatedMatching.geneName(INDnew);
    
    
    for s=1:6
        DataTable.Expression{s,1} = DataTable.Expression{s,1}(INDold,:);
        DataTable.Noise{s,1} = DataTable.Noise{s,1}(INDold,:);
    end
    
else
    % if chosen not to update probes, then remove ones with missing entrezIDs
    fprintf(1,'Removing probes not mapped to genes\n')
    fprintf(1,'Removing %d irrelevant probes\n', sum(isnan(EntrezID)))
    %------------------------------------------------------------------------------
    % Make a table from all the data
    %------------------------------------------------------------------------------
    % delete probes that were excluded during the quality control from expression and noise matrices
    for s=1:6
        DataTable.Expression{s,1}(isnan(EntrezID),:) = [];
        DataTable.Noise{s,1}(isnan(EntrezID),:) = [];
    end
    
    ProbeName(isnan(EntrezID)) = [];
    GeneSymbol(isnan(EntrezID)) = [];
    ProbeID(isnan(EntrezID)) = [];
    EntrezID(isnan(EntrezID)) = [];
end

%------------------------------------------------------------------------------
% Assign ProbeIDs, EntrezIDs and ProbeNames to Data cell.
%------------------------------------------------------------------------------

DataProbe{1,1} = ProbeID;
DataProbe{1,2} = EntrezID;
DataProbe{1,3} = ProbeName;
DataProbe{1,4} = GeneSymbol;
DataTableProbe = dataset({DataProbe, headerprobe{:}});
%------------------------------------------------------------------------------
% Combine data for all subjects
%------------------------------------------------------------------------------
fprintf(1,'Combining data for all subjects\n')
Expressionall = horzcat(DataTable{1,1}, DataTable{2,1}, DataTable{3,1}, DataTable{4,1}, DataTable{5,1}, DataTable{6,1});
Coordinatesall = vertcat(DataTable{1,2}, DataTable{2,2}, DataTable{3,2}, DataTable{4,2}, DataTable{5,2}, DataTable{6,2});
StructureNamesall = vertcat(DataTable{1,3}, DataTable{2,3}, DataTable{3,3}, DataTable{4,3}, DataTable{5,3}, DataTable{6,3});
MRIvoxCoordinatesAll = vertcat(DataTable{1,4}, DataTable{2,4}, DataTable{3,4}, DataTable{4,4}, DataTable{5,4}, DataTable{6,4});
noiseall = horzcat(DataTable{1,5}, DataTable{2,5}, DataTable{3,5}, DataTable{4,5}, DataTable{5,5}, DataTable{6,5});

%------------------------------------------------------------------------------
% Save relevant variables to a MicroarrayData.mat file
%------------------------------------------------------------------------------
cd ..
cd ('processedData');

fprintf(1,'Saving data to the file\n')
save(sprintf('%sProbesUpdated.mat', startFileName), 'DataTable','DataTableProbe', 'Expressionall', 'Coordinatesall', 'StructureNamesall', 'MRIvoxCoordinatesAll', 'noiseall');
cd ../../..
