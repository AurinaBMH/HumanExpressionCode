
%% Author: Aurina
%% Date modified: 2016-03-22
%% Date modified: 2017-07-14
%% Date modified: 2017-08-29 - noise level data added
%% Date modified: 2017-11-10 - manual gene symbol naming fixed
%% This script:
%   1. Loads all microarray data from excell files for each subject
%   2. Excludes custom probes;
%   3. Excludes probes with missing entrezIDs
%   4. Saves expression data, coordinates, sample structure names for all samples
%   5. Saves data for separate subjects to DataTable
%   6. Saves data for all subjects combined as variables 'MicorarrayData.mat' file
%%
%------------------------------------------------------------------------------
% Choose options
%------------------------------------------------------------------------------
% RemoveCUSTProbes = true; will exclude CUST probes
% ExcludeCBandBS = true; will exclude samples from brainstem and cerebellum
useCUSTprobes = true;
ExcludeCBandBS = true;

if useCUSTprobes
    fprintf(1,'CUST probes will be included\n')
    startFileName = 'MicroarrayDataWITHcust';
else
    fprintf(1,'CUST probes will be excluded\n')
    startFileName = 'MicroarrayData';
end

cd ('data/genes/rawData');
%% load probe information (same for all subjects)
fprintf(1,'Loading Probes.xlsx file\n')
FileProbes = 'Probes.xlsx';
% in Probes.xlsx file probes from 1070378 down are deleted as there are no information about gene samples.
ProbeTable = readtable(FileProbes);
ProbeID = ProbeTable.probe_id;
EntrezID = ProbeTable.entrez_id;
ProbeName =  ProbeTable.probe_name;
GeneID = ProbeTable.gene_id;
GeneSymbol = ProbeTable.gene_symbol;
GeneName = ProbeTable.gene_name;


% assign NaN values for all probes with missing entrezIDs
% this is the final list of probes to be used in max var calculations
ProbeID(isnan(EntrezID)) = NaN;
fprintf(1,'%d probes with missing entrez IDs\n', sum(isnan(EntrezID)));
% creat a Data cell to store the output
headerdata = {'Expression' , 'MMcoordinates', 'StructureName', 'MRIvoxCoordinates', 'Noise'};
headerprobe = { 'ProbeID', 'EntrezID','ProbeName', 'GeneID', 'GeneSymbol', 'GeneName'};
Data = cell(6,5);
DataProbe = cell(1,6);
%%
%------------------------------------------------------------------------------
%Go to each subject's directory and take the data
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
    % rows from 57860 are removed corresponding to probes from 1070378 down as there are no information about gene samples.
    % Expression = removerows(Expression,57860:size(Expression,1));
    Expression(:,1) = [];                         % exclude probe IDs from expression matrix
    [~,~,SlabType] = xlsread(FileAnnot, 'D:D');
    [~,~, StructureName] = xlsread(FileAnnot, 'F:F');
    SlabType(1) = [];                           % remove headline
    StructureName(1) = [];                      % remove headline
    MMcoordinates = xlsread(FileAnnot, 'K:M');
    MRIvoxCoordinates = xlsread(FileAnnot, 'H:J');
    [~,probeList] = intersect(noise(:,1),ProbeID, 'stable');
    noise = noise(probeList,2:end);
    % keep only non custom probes with existing entrezIDs
    Expression(isnan(ProbeID),:) = [];
    
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
        noise(:,BSandCBind) = NaN;
        StructureName(BSandCBind) = {NaN};
    end
    
    
    % for nan columns
    % keep only existing expression values
    Expression = Expression(:,all(~isnan(Expression)));
    % keep only existing coordinates
    MMcoordinates = MMcoordinates(all(~isnan(MMcoordinates),2),:); % for nan rows
    MRIvoxCoordinates = MRIvoxCoordinates(all(~isnan(MRIvoxCoordinates),2),:); % for nan rows
    noise = noise(:,all(~isnan(noise)));
    % keep only existing structure names
    StructureName(cellfun(@(StructureName) any(isnan(StructureName)),StructureName)) = [];
    
    % assign output to Data cell;
    Data{subj,1} = Expression;
    Data{subj,2} = MMcoordinates;
    Data{subj,3} = StructureName;
    Data{subj,4} = MRIvoxCoordinates;
    Data{subj,5} = noise;
    cd ..
    
end
%%
%------------------------------------------------------------------------------
% Keep only existing ProbeNames, EntrezIDs and ProbeIDs and other gene related information.
%------------------------------------------------------------------------------
fprintf(1,'Removing irrelevant probes\n')
fprintf(1,'Removing %d irrelevant probes\n', sum(isnan(ProbeID)))

ProbeName(isnan(ProbeID)) = [];
EntrezID(isnan(ProbeID)) = [];
GeneID(isnan(ProbeID)) = [];
GeneSymbol(isnan(ProbeID)) = [];
GeneName(isnan(ProbeID)) = [];
ProbeID(isnan(ProbeID)) = [];

%------------------------------------------------------------------------------
% Check for mistakes in probe naming and fix them plus update expression values for those probes that are excluded
%------------------------------------------------------------------------------
% load ncbi data
load('HomoSampiens_geneInfo20171111.mat'); 

[matches, GeneSymbol, GeneName, nrUpdated, checkEntrezID] = checkGene(Homosapiens, GeneSymbol, GeneName, EntrezID); 

onlyExisting = matches;
E1 = EntrezID; 
onlyExisting(any(isnan(matches), 2), :) = [];
E1(any(isnan(matches), 2), :) = [];

ind = find(sum(onlyExisting,2)==0); 
nrMissing = unique(E1(ind)); 

fprintf('%d genes where entrezID exists but both name and symbol doesnt match\n', length(nrMissing));

% save a list of entrezIDs that were not found 
fileID = fopen('entrezIDs.txt','w');
fprintf(fileID, '%d\n',checkEntrezID);
fclose(fileID);

% check if entrezIDs that were not found were discontinued
% run this from the terminal
awk 'NR==FNR{pats[$0]; next} $3 in pats' entrezIDs.txt gene_history.txt > listOfentrezIDs.txt

%load the list 
load ('listOfentrezIDs.txt')


% remove those that were discontinued and replace those that are changed
% with new ones. 


% write code here
% run the function again to make changes for new entries. 
[matches, GeneSymbol, GeneName, nrUpdated, checkEntrezID] = checkGene(Homosapiens, GeneSymbol, GeneName, EntrezID);




% after that check if each uniwue gene ID has a unique EntrezID and the
% other way around using selectDuplicates(IDone, IDtwo, changetoNumber).



