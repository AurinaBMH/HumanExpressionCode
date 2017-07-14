%% Author: Aurina
%% Date modified: 2016-03-22
%% This script:
%   1. Loads all microarray data from excell files for each subject
%   2. Excludes custom probes;
%   3. Excludes probes with missing entrezIDs
%   4. Saves expression data, coordinates, sample structure names for all samples
%   5. Saves data for separate subjects to DataTable
%   6. Saves data for all subjects combined as variables 'MicorarrayData.mat' file
%%
clear all;

%choose options (if you'd like to exclude braintem and cerebellum data
%points before choosing most relevant probes: ExcludeCBandBS = 1; if not ExcludeCBandBS = 0;
RemoveCUSTProbes = 1;
ExcludeCBandBS = 0;


        cd ('/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Microarray/');
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
        %Remove probes:
        %------------------------------------------------------------------------------
        if RemoveCUSTProbes ==1
        % Remove all CUST probes (assign NaN values for all custom probes)
            fprintf(1,'Removing CUST probes\n')
            for prob=1:size(ProbeName)
                Probe = ProbeName{prob}; 
                if strcmp(Probe(1),'C')
                   ProbeName{prob} = NaN; 
                   ProbeID(prob) = NaN;
                end
            end
        end
        % assign NaN values for all probes with missing entrezIDs
        % this is the final list of probes to be used in max var calculations
        ProbeID(isnan(EntrezID)) = NaN;
        % creat a Data cell to store the output
        headerdata = {'Expression' , 'MMcoordinates', 'StructureName', 'MRIvoxCoordinates'};
        headerprobe = { 'ProbeID', 'EntrezID','ProbeName', 'GeneID', 'GeneSymbol', 'GeneName'};
        Data = cell(6,4);
        DataProbe = cell(1,6);

        %% go to each subject's directory and take the data
        for subj=1:6
            fprintf(1,'Loading data for %u subject\n', subj)
            folder = sprintf('normalized_microarray_donor0%d', subj);
            cd (folder);
                %% load information specific for each subject
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
                if ExcludeCBandBS == 1
                    fprintf('Excluding brainstem and cerebellum data\n')
                    for roi = 1:size(Expression,2)
                        if strcmp(SlabType{roi}, 'BS') || strcmp(SlabType{roi}, 'CB')
                            Expression(:,roi) = NaN;
                            MMcoordinates(roi,:) = NaN;
                            MRIvoxCoordinates(roi,:) = NaN;
                            StructureName{roi} = NaN;
                        end
                    end
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
                 cd ('/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Microarray/');
        % 
        end
        %% keep only existing ProbeNames, EntrezIDs and ProbeIDs and other gene related information.
            fprintf(1,'Removing irrelevant probes\n')
            
            ProbeName(isnan(ProbeID)) = [];

            EntrezID(isnan(ProbeID)) = [];
  
            GeneID(isnan(ProbeID)) = [];

            GeneSymbol(isnan(ProbeID)) = [];
  
            GeneName(isnan(ProbeID)) = [];
            
            ProbeID(isnan(ProbeID)) = [];

        %% assign ProbeIDs, EntrezIDs and ProbeNames to Data cell.   

            DataProbe{1,1} = ProbeID;
            DataProbe{1,2} = EntrezID;
            DataProbe{1,3} = ProbeName;
            DataProbe{1,4} = GeneID;
            DataProbe{1,5} = GeneSymbol;
            DataProbe{1,6} = GeneName;
  

        %% make a table from all the data
        DataTable = dataset({Data, headerdata{:}});
        DataTableProbe = dataset({DataProbe, headerprobe{:}});

        %% combine expression and coordinate values for all subjects
        fprintf(1,'Combining data for all subjects\n')
        Expressionall = horzcat(DataTable{1,1}, DataTable{2,1}, DataTable{3,1}, DataTable{4,1}, DataTable{5,1}, DataTable{6,1});
        Coordinatesall = vertcat(DataTable{1,2}, DataTable{2,2}, DataTable{3,2}, DataTable{4,2}, DataTable{5,2}, DataTable{6,2});
        StructureNamesall = vertcat(DataTable{1,3}, DataTable{2,3}, DataTable{3,3}, DataTable{4,3}, DataTable{5,3}, DataTable{6,3});
        MRIvoxCoordinatesAll = vertcat(DataTable{1,4}, DataTable{2,4}, DataTable{3,4}, DataTable{4,4}, DataTable{5,4}, DataTable{6,4});

        %% save relevant variables to a MicroarrayData.mat file
        if ~RemoveCUSTProbes
         fprintf(1,'Saving data with CUST probes to the file\n')
         save('MicroarrayDataWITHCUST.mat', 'DataTable','DataTableProbe', 'Expressionall', 'Coordinatesall', 'StructureNamesall', 'MRIvoxCoordinatesAll');
         clear all; 
        else
         fprintf(1,'Saving data without CUST probes to the file\n')
         save('MicroarrayData.mat', 'DataTable','DataTableProbe', 'Expressionall', 'Coordinatesall', 'StructureNamesall', 'MRIvoxCoordinatesAll');
         clear all;
        end