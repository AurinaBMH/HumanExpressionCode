%% Author: Aurina
%Last modiffied: 2016-04-15
%Last modiffied: 2017-07-20

% This script loads MicroarrayDataPCAS[subject].mat file and
% separates samples according to Allen labels to right/left, cortex/subcortex
% resulting in 1 file for each subject.

%------------------------------------------------------------------------------
% Choose what to separate
%------------------------------------------------------------------------------
% separate according to sides
useCUSTprobes = true;
probeSelection = 'Variance';% (Variance', LessNoise', 'Mean')
sides = {'right', 'left'};
% separate according to brain part: cortex/subcortex
brainParts = {'Cortex', 'Subcortex'};
% loop throgh all 6 subjects
subjects = 1:6;
cd ('data/genes/processedData')
%------------------------------------------------------------------------------
% Separate samples according to Allen institute Structure names.
%------------------------------------------------------------------------------

if useCUSTprobes
    fprintf(1,'Using the data with CUST probes\n')
    startFileName = 'MicroarrayDataWITHcust';
else
    fprintf(1,'Using the data without CUST probes\n')
    startFileName = 'MicroarrayData';
end

for subject = subjects
    % load data file
    
    FileName = sprintf('%s%sS0%d.mat',startFileName, probeSelection, subject);
    
    load(FileName);
    
    % select words that belong to cortex
    cortex = {'cortex', 'gyrus', 'gyri', 'sulcus', 'occipital pole', 'planum temporale', 'cuneus', ...
        'frontal pole', 'operculum', 'planum polare', 'temporal pole' , 'paracentral lobule' };%
    
    % sepect samples according to parts and sides
    index = cell(length(cortex),1);
    for side = sides
        for brainPart = brainParts
            % rename variables to keep the original
            expressionS = Expression;
            sampleInformationS = SampleInformation;
            
            % chack each structure name for side and brain pars separation
            for i=1:length(sampleInformationS.StructureNames)
                Structure = sampleInformationS.StructureNames{i};
                
                %check if structure name has a part related to cortex
                for j = 1:length(cortex)
                    index{j} = strfind(Structure, cortex{j});
                end
                
                %check if cortex related name part was found
                emptyCells = cellfun(@isempty,index);
                
                if strcmp(brainPart{1}, 'Cortex') % will exclude subcortical samples
                    
                    if sum(emptyCells) == length(cortex)% && isempty(index2)
                        expressionS(i,:) = NaN;
                        sampleInformationS.MRIvoxCoordinates(i,:) = NaN;
                        sampleInformationS.MMCoordinates(i,:) = NaN;
                        sampleInformationS.StructureNames{i} = 'remove';
                    end
                    
                elseif strcmp(brainPart{1}, 'Subcortex') % will exclude cortical samples
                    
                    if sum(emptyCells) < length(cortex)% && isempty(index2)
                        expressionS(i,:) = NaN;
                        sampleInformationS.MRIvoxCoordinates(i,:) = NaN;
                        sampleInformationS.MMCoordinates(i,:) = NaN;
                        sampleInformationS.StructureNames{i} = 'remove';
                        
                    end
                end
                % there are several samples that are not labeled 'left',
                % 'right' - like corpus {'corpus callosum'}, so we exclude them.
                index2 = strfind(Structure, side{1});
                if isempty(index2)
                    expressionS(i,:) = NaN;
                    sampleInformationS.MRIvoxCoordinates(i,:) = NaN;
                    sampleInformationS.MMCoordinates(i,:) = NaN;
                    sampleInformationS.StructureNames{i} = 'remove'; %"; %{double.empty(0)};
                end
            end
            
            % exclude nonexisting data
            expressionS(any(isnan(expressionS),2),:) = [];
            expression.(side{1}).(brainPart{1}) = expressionS;
            
            sampleInformationS.MRIvoxCoordinates(any(isnan(sampleInformationS.MRIvoxCoordinates),2),:) = [];
            sampleInformation.(side{1}).(brainPart{1}).MRIvoxCoordinates = sampleInformationS.MRIvoxCoordinates;
            
            sampleInformationS.MMCoordinates(any(isnan(sampleInformationS.MMCoordinates),2),:) = [];
            sampleInformation.(side{1}).(brainPart{1}).MMCoordinates = sampleInformationS.MMCoordinates ;
            
            sampleInformationS.StructureNames(strcmp('remove',sampleInformationS.StructureNames)) = [];
            sampleInformation.(side{1}).(brainPart{1}).StructureNames = sampleInformationS.StructureNames;
            
            probeInformation = ProbeInformation;
        end
    end
    %------------------------------------------------------------------------------
    % Save data for each subject separately
    %------------------------------------------------------------------------------
    
    SaveFileName = sprintf('%s%ssepS0%d.mat',startFileName, probeSelection, subject);
    
    save (SaveFileName, 'expression', 'probeInformation', 'sampleInformation');
end
cd ../../..

