%% Author: Aurina
%% Date modified: 2016-03-22
%% Date modified: 2017-07-14
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
    startFileName = 'MicroarrayDataWITHcust';
else
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
%------------------------------------------------------------------------------
% Remove CUST probes:
%------------------------------------------------------------------------------
        if ~useCUSTprobes
        % Remove all CUST probes (assign NaN values for all custom probes)
            fprintf(1,'Removing CUST probes\n')
            cust = strfind(ProbeName, 'CUST'); 
            remInd = find(~cellfun(@isempty,cust)); 
            fprintf(1,'%d CUST probes removed\n', length(remInd))
            ProbeName(remInd) = {NaN}; 
            ProbeID(remInd) = NaN;
        end

        % assign NaN values for all probes with missing entrezIDs
        % this is the final list of probes to be used in max var calculations
        ProbeID(isnan(EntrezID)) = NaN;
        fprintf(1,'%d probes with missing entrez IDs\n', sum(isnan(EntrezID))); 
        % creat a Data cell to store the output
        headerdata = {'Expression' , 'MMcoordinates', 'StructureName', 'MRIvoxCoordinates'};
        headerprobe = { 'ProbeID', 'EntrezID','ProbeName', 'GeneID', 'GeneSymbol', 'GeneName'};
        Data = cell(6,4);
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
                Expression = csvread(FileMicroarray);
                % rows from 57860 are removed corresponding to probes from 1070378 down as there are no information about gene samples.   
                % Expression = removerows(Expression,57860:size(Expression,1));
                Expression(:,1) = [];                         % exclude probe IDs from expression matrix
                [~,~,SlabType] = xlsread(FileAnnot, 'D:D');
                [~,~, StructureName] = xlsread(FileAnnot, 'F:F');
                SlabType(1) = [];                           % remove headline
                StructureName(1) = [];                      % remove headline
                MMcoordinates = xlsread(FileAnnot, 'K:M');
                MRIvoxCoordinates = xlsread(FileAnnot, 'H:J');

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
                    StructureName(BSandCBind) = {NaN};
                end
                    
                 % for nan columns 
                 % keep only existing expression values
                 Expression = Expression(:,all(~isnan(Expression))); 
                 % keep only existing coordinates
                 MMcoordinates = MMcoordinates(all(~isnan(MMcoordinates),2),:); % for nan rows
                 MRIvoxCoordinates = MRIvoxCoordinates(all(~isnan(MRIvoxCoordinates),2),:); % for nan rows
                 % keep only existing structure names
                 StructureName(cellfun(@(StructureName) any(isnan(StructureName)),StructureName)) = [];                  
                 
                % assign output to Data cell;
                 Data{subj,1} = Expression;
                 Data{subj,2} = MMcoordinates;
                 Data{subj,3} = StructureName;
                 Data{subj,4} = MRIvoxCoordinates;
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
% Assign ProbeIDs, EntrezIDs and ProbeNames to Data cell.   
%------------------------------------------------------------------------------ 

            DataProbe{1,1} = ProbeID;
            DataProbe{1,2} = EntrezID;
            DataProbe{1,3} = ProbeName;
            DataProbe{1,4} = GeneID;
            DataProbe{1,5} = GeneSymbol;
            DataProbe{1,6} = GeneName;
 
%------------------------------------------------------------------------------
% Make a table from all the data
%------------------------------------------------------------------------------ 
        DataTable = dataset({Data, headerdata{:}});
        DataTableProbe = dataset({DataProbe, headerprobe{:}});

%------------------------------------------------------------------------------
% Combine data for all subjects
%------------------------------------------------------------------------------ 
        fprintf(1,'Combining data for all subjects\n')
        Expressionall = horzcat(DataTable{1,1}, DataTable{2,1}, DataTable{3,1}, DataTable{4,1}, DataTable{5,1}, DataTable{6,1});
        Coordinatesall = vertcat(DataTable{1,2}, DataTable{2,2}, DataTable{3,2}, DataTable{4,2}, DataTable{5,2}, DataTable{6,2});
        StructureNamesall = vertcat(DataTable{1,3}, DataTable{2,3}, DataTable{3,3}, DataTable{4,3}, DataTable{5,3}, DataTable{6,3});
        MRIvoxCoordinatesAll = vertcat(DataTable{1,4}, DataTable{2,4}, DataTable{3,4}, DataTable{4,4}, DataTable{5,4}, DataTable{6,4});
        
%------------------------------------------------------------------------------
% Save relevant variables to a MicroarrayData.mat file
%------------------------------------------------------------------------------ 
        cd ..
        cd ('processedData');

         fprintf(1,'Saving data with CUST probes to the file\n')
         save(sprintf('%s.mat', startFileName), 'DataTable','DataTableProbe', 'Expressionall', 'Coordinatesall', 'StructureNamesall', 'MRIvoxCoordinatesAll');
  