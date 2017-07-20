%% Author: Aurina
%Last modiffied: 2016-04-14
%Last modiffied: 2017-07-20

% It assigns samples to parcellation separately for cortex/subcortex,
% right/left sides using only relevant parcellation voxels to make an assignment.

% closest points are chosen using dsearchn function, which assignes a parcellation coordinate index to each microarray sample coordinate.
% output: original microarray coordinates to be assigned (coordsToAssign)
% and new assigned coordinates (coordsAssigned) and expression data are saved in
% each subject's folder separately for each parcellation and left/right cortex/subcortex.


%------------------------------------------------------------------------------
% Clear grounds
%------------------------------------------------------------------------------


%% choose if you want to use data with CUST probes
useCUSTprobes = false;

%% choose what type of probe selection to use, hemisphere, subject list, parcellations, threshols.
probeSelection = 'PC';% (Variance', LessNoise', 'Mean')

% UsePCAprobes=1;
% UseMaxVarprobes=0;

%LeftHemisphere = 0; % Hemisphere = 1 - left; Hemisphere = 0 - right.
%OnlyCortex = 1; % OnlyCortex = 1 - only cortex; OnlyCortex = 0 - only subcortex;

%if LeftHemisphere == 1
% subjects = [1 2 3 4 5 6]; % for subjects 1,2 run with both hemispheres; for 3:6 run only with left.
%elseif LeftHemisphere == 0
% subjects = [1 2];
%end
parcellations = {'aparcaseg'};%, 'cust100', 'cust250'};
distanceThreshold = 30;
subjects = 1:6;

%% select variables according to side/brain part selections
sides = {'left', 'right'};
brainParts = {'Cortex', 'Subcortex'};
% if LeftHemisphere == 1
%     Side = Sides{1};
% elseif LeftHemisphere == 0
%     Side = Sides{2};
% else
%     error('Incorrect index for choosing hemisphere: must be 1 or 0');
% end


%
% if OnlyCortex == 1
%     BrainPart = brainParts {1};
% elseif OnlyCortex == 0
%     BrainPart = brainParts {2};
% else
%     error('Incorrect index for choosing cortexVS subcortex: must be 1 or 0');
% end


%% define subject numbers, parcellations and assignment distance thresholds

% Subj = cell(6,1);
% Subj{1} = '/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Microarray/S01_H0351_2001';
% Subj{2} = '/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Microarray/S02_H0351_2002';
% Subj{3} = '/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Microarray/S03_H0351_1009';
% Subj{4} = '/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Microarray/S04_H0351_1012';
% Subj{5} = '/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Microarray/S05_H0351_1015';
% Subj{6} = '/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Microarray/S06_H0351_1016';


for subject = subjects
    cd ('data/genes/parcellations')
    subjectDir = sprintf('S0%d_H0351', subject);
    cd (subjectDir)
    for parcellation = parcellations
        
        %for DistanceThreshold = distanceThreshold
        
        fprintf('Subject %u parcellation %s assignment distance threshold %u\n; ', subject, parcellation{1}, distanceThreshold )
        
        
        %------------------------------------------------------------------------------
        % Load parcellations
        %------------------------------------------------------------------------------%%
        %S = Subj{subject};
        % load the parcellation
        if strcmp(parcellation, 'aparcaseg')
            Folder = 'default_NativeAnat';
            %FolderName = strcat(S, Folder);
            cd (Folder);
            [~, data_parcel]=read2('defaultparc_NativeAnat.nii');
            NumNodes = 82;
            LeftCortex = 34;
            LeftSubcortex = 41;
            RightCortex = 75;
            RightSubcortex = NumNodes;
        elseif strcmp(parcellation, 'cust100')
            Folder = 'custom100_NativeAnat';
            %FolderName = strcat(S, Folder);
            cd (FolderName);
            [~, data_parcel]=read2('customparc_NativeAnat.nii');
            NumNodes = 220;
            LeftCortex = 100;
            LeftSubcortex = 110;
            RightCortex = 210;
            RightSubcortex = NumNodes;
        elseif strcmp(parcellation, 'cust250')
            Folder = 'custom250_NativeAnat';
            %FolderName = strcat(S, Folder);
            cd (FolderName);
            [~, data_parcel]=read2('customparc_NativeAnat.nii');
            NumNodes = 530;
            LeftCortex = 250;
            LeftSubcortex = 265;
            RightCortex = 515;
            RightSubcortex = NumNodes;
            
        end
        cd ../../../
        %%
        %------------------------------------------------------------------------------
        % Load microarray data
        %------------------------------------------------------------------------------
        cd ('processedData');
        
        if  useCUSTprobes
            fprintf('Loading MicroarrayDataWITHcust%ssepS0%d.mat\n', probeSelection, subject)
            load(sprintf('MicroarrayDataWITHcust%ssepS0%d.mat', probeSelection, subject));
        else
            fprintf('Loading MicroarrayData%ssepS0%d.mat\n', probeSelection, subject)
            load(sprintf('MicroarrayData%ssepS0%d.mat', probeSelection, subject));
        end
        
        for side = sides
            for brainPart = brainParts
                
                coords2assign = sampleInformation.(side{1}).(brainPart{1}).MRIvoxCoordinates;
                
                if distanceThreshold <30

                    load('CoordsAssignedAll.mat');
                    cd ../../../
                end
                %%
                %------------------------------------------------------------------------------
                % Find coordinates for nonzero elements in parcellation according to side and brain part (will be used to assign microarray samples to)
                %------------------------------------------------------------------------------
                % get coordinate values
                if strcmp(side, 'left') && strcmp(brainPart, 'Cortex')
                    Text = sprintf('Subject %d LEFT cortex\n', subject);
                    fprintf(1,Text)
                    [Coordx,Coordy,Coordz] = ind2sub(size(data_parcel),find(data_parcel > 0 & data_parcel <= LeftCortex));
                elseif strcmp(side, 'left') && strcmp(brainPart, 'Subcortex')
                    Text = sprintf('Subject %d LEFT subcortex\n', subject);
                    fprintf(1,Text)
                    [Coordx,Coordy,Coordz] = ind2sub(size(data_parcel),find(data_parcel > LeftCortex & data_parcel <= LeftSubcortex));
                elseif strcmp(side, 'right') && strcmp(brainPart, 'Cortex')
                    Text = sprintf('Subject %d RIGHT cortex\n', subject);
                    fprintf(1,Text)
                    [Coordx,Coordy,Coordz] = ind2sub(size(data_parcel),find(data_parcel > LeftSubcortex & data_parcel <= RightCortex));
                elseif strcmp(side, 'right') && strcmp(brainPart, 'Subcortex')
                    Text = sprintf('Subject %d RIGHT subcortex\n', subject);
                    fprintf(1,Text)
                    [Coordx,Coordy,Coordz] = ind2sub(size(data_parcel),find(data_parcel > RightCortex & data_parcel <= RightSubcortex));
                end
                
                coordsNonzeroParcel = cat(2,Coordx,Coordy,Coordz);
                %%
                %------------------------------------------------------------------------------
                % For each microarray coordinate find a closest coordinate in parcellation
                %------------------------------------------------------------------------------
                
                coordsAssigned =  zeros(size(coords2assign));
                % finds closest point
                % T = delaunayn(unique(coordsNonzeroParcel));
                k = dsearchn(coordsNonzeroParcel,coords2assign);
                
                for i = 1:length(k)
                    coordsAssigned(i,:) = coordsNonzeroParcel(k(i),:);
                end
                
                %------------------------------------------------------------------------------
                % Salculate the distance between original and reassigned coordinate
                %------------------------------------------------------------------------------
                assignDistance = zeros(length(coordsAssigned),1);
                coordsNONassigned = zeros(length(coordsAssigned),3);
                coordsToAssignAll = coords2assign;
                
                for j=1:length(coordsAssigned)
                    
                    assignDistance(j,1) = pdist2(coordsAssigned(j,:), coords2assign(j,:));
                    
                    if assignDistance(j)>distanceThreshold
                        coordsNONassigned(j,:,:,:)=coords2assign(j,:,:,:);
                        coordsAssigned(j,:)=0;
                        coords2assign(j,:)=0;
                        
                    end
                    
                end
                %remove zero elements
                coords2assign( ~any(coords2assign,2), : ) = NaN;
                coordsAssigned( ~any(coordsAssigned,2), : ) = NaN;
                coordsNONassigned( ~any(coordsNONassigned,2), : ) = [];
                
                %------------------------------------------------------------------------------
                % Extract intensity values for assigned coordinates
                %------------------------------------------------------------------------------
                
                if distanceThreshold == 30
                    
                    %cd ('processedData')
                    coordsAssignedALL.(side{1}).(brainPart{1}) = coordsAssigned;
                    save(sprintf('CoordsAssignedAllS0%d.mat', subject), 'coordsAssignedALL');
                
                
                elseif distanceThreshold <30
                    
                    Int=nonzeros(data_parcel);
                    int = unique(Int);
                    
                    intensity_all=zeros(size(coordsAssignedALL.(side{1}).(brainPart{1}),1), 1);
                    
                    for i=1:size(coordsAssignedALL.(side{1}).(brainPart{1}),1)
                        
                        intensity = data_parcel(coordsAssignedALL.(side{1}).(brainPart{1})(i,1), coordsAssignedALL.(side{1}).(brainPart{1})(i,2), coordsAssignedALL.(side{1}).(brainPart{1})(i,3));
                        intensity_all(i) = intensity;
                        
                    end
                    
                    %% sort nonzero coordinates according to intensity value
                    full_informationALL = [coordsAssigned intensity_all];
                    full_information = full_informationALL;
                    full_informationALL(any(isnan(full_informationALL),2),:) = NaN;
                    full_information(any(isnan(full_informationALL),2),:) = [];
                    
                    full_information_sorted = sortrows(full_information,4);
                    
                    sampleIntensity = full_informationALL(:,4);
                    IntANDExpression = [sampleIntensity Expression];
                    % exclude nonassigned points
                    IntANDExpression(any(isnan(IntANDExpression),2),:) = NaN;
                    IntANDExpression(any(isnan(IntANDExpression),2),:) = [];
                    
                    fprintf(1,'Sorts expression data according to ROIs\n')
                    
                    % sort according to intensity (roi)
                    IntANDExpressionSorted = sortrows(IntANDExpression,1);
                    % assign outputs to structure
                    information.(side{1}).(brainPart{1}).raw = full_information;
                    information.(side{1}).(brainPart{1}).sorted = full_information_sorted;
                    coordinates.(side{1}).(brainPart{1}).assigned = coordsAssigned;
                    coordinates.(side{1}).(brainPart{1}).toassign = coords2assign;
                    coordinates.(side{1}).(brainPart{1}).toassignALL = coordsToAssignAll;
                    coordinates.(side{1}).(brainPart{1}).nonassigned = coordsNONassigned;
                    
                    data.(side{1}).(brainPart{1}).expression = IntANDExpressionSorted;
                    data.(side{1}).(brainPart{1}).information = full_information_sorted;
                    %------------------------------------------------------------------------------
                    % Save output
                    %------------------------------------------------------------------------------
                    
                    %save(sprintf('%d_DistThresh%d_S0%d_MRIvoxCoordsAssigned.mat', NumNodes, distanceThreshold, subject), ...
                      %  'coordinates', 'information');
                    
                    %save (sprintf('%d_DistThresh%d_S0%d_ExpressionProbe%s.mat', NumNodes, distanceThreshold, subject, probeSelection ), 'data' );
                    %savefig (L, sprintf('%d_DistThresh%d_S0%d_ExpressionProbeMaxVar.fig', NumNodes, Threshold, subject));
                    
                    
                    %savefig (L, sprintf('%d_DistThresh%d_S0%d_ExpressionProbePCA_%s%s.fig', NumNodes, Threshold, subject, Side, BrainPart));
                end
                
                %% plot
                %     figure; scatter3(coordsAssigned(:,1), coordsAssigned(:,2), coordsAssigned(:,3), 'bo');...
                %         hold on; scatter3(coordsToAssign(:,1), coordsToAssign(:,2), coordsToAssign(:,3), 'g*');...
                %         hold on; scatter3(coordsNonassigned(:,1), coordsNonassigned(:,2), coordsNonassigned(:,3), 'r.');
                %         hold on;
                %         % plot arrows where the point moved
                %         for j = 1:length(coordsToAssign)
                %             dp = coordsToAssign(j,:) - coordsAssigned(j,:);
                %             q = quiver3(coordsToAssign(j,1),coordsToAssign(j,2),coordsToAssign(j,3),-dp(1), -dp(2), -dp(3), 'LineWidth', 1.5);
                %             hold on;
                %         end
                
                %savefig(sprintf('%d_DistThresh%d_coordMapping_MRI.fig', NumNodes, DistanceThreshold ));
                
                %close all;
            end
            
        end
        
    end
    cd ../../..
end

