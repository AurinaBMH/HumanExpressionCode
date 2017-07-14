%% Script to assign Allen institute microarray sample data to the closest point on a brainmask
% Last edit: Aurina 2016-04-15
% It assigns samples to parcellation separately for cortex/subcortex,
% right/left sides using only relevant parcellation voxels to make an assignment.

% closest points are chosen using dsearchn function, which assignes a parcellation coordinate index to each microarray sample coordinate. 
% output: original microarray coordinates to be assigned (coordsToAssign)
% and new assigned coordinates (coordsAssigned) and expression data are saved in
% each subject's folder separately for each parcellation and left/right cortex/subcortex. 

%% Clear grounds 
clear all; close all;

%% choose if you want to use data with CUST probes 
UseCUSTprobes = 0;
Thr = 0;
%% choose what type of probe selection to use, hemisphere, subject list, parcellations, threshols. 
UsePCAprobes=1;
UseMaxVarprobes=0;

LeftHemisphere = 0; % Hemisphere = 1 - left; Hemisphere = 0 - right.
OnlyCortex = 1; % OnlyCortex = 1 - only cortex; OnlyCortex = 0 - only subcortex; 

if LeftHemisphere == 1
    subjects = [1 2 3 4 5 6]; % for subjects 1,2 run with both hemispheres; for 3:6 run only with left. 
elseif LeftHemisphere == 0
    subjects = [1 2];
end
Parcellations = {'aparcaseg', 'cust100', 'cust250'};%, 'cust100', 'cust250'};
ThresholdList = [2];

%% select variables according to side/brain part selections
Sides = {'left', 'right'};
if LeftHemisphere == 1
    Side = Sides{1};
elseif LeftHemisphere == 0
    Side = Sides{2};
else
    error('Incorrect index for choosing hemisphere: must be 1 or 0');
end

BrainParts = {'Cortex', 'Subcortex'};

if OnlyCortex == 1
    BrainPart = BrainParts {1};
elseif OnlyCortex == 0
    BrainPart = BrainParts {2};
else
    error('Incorrect index for choosing cortexVS subcortex: must be 1 or 0');
end


%% define subject numbers, parcellations and assignment distance thresholds

Subj = cell(6,1);
Subj{1} = '/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Microarray/S01_H0351_2001';
Subj{2} = '/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Microarray/S02_H0351_2002';
Subj{3} = '/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Microarray/S03_H0351_1009';
Subj{4} = '/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Microarray/S04_H0351_1012';
Subj{5} = '/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Microarray/S05_H0351_1015';
Subj{6} = '/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Microarray/S06_H0351_1016';


for subject = subjects

    for Parcellation = Parcellations
        
        for DistanceThreshold = ThresholdList
            
        fprintf('Subject %u Parcellation %s Assignment distance threshold %u\n; ', subject, Parcellation{1}, DistanceThreshold )

        %% load the data for parcellation
            S = Subj{subject};
            % load the parcellation
            if strcmp(Parcellation, 'aparcaseg')
                Folder = '/default_NativeAnat';
                FolderName = strcat(S, Folder);
                cd (FolderName);
                [~, data_parcel]=read('defaultparc_NativeAnat.nii');
                NumNodes = 82;
                LeftCortex = 34;
                LeftSubcortex = 41; 
                RightCortex = 75;
                RightSubcortex = NumNodes;
            elseif strcmp(Parcellation, 'cust100')
                Folder = '/custom100_NativeAnat';
                FolderName = strcat(S, Folder);
                cd (FolderName);
                [~, data_parcel]=read('customparc_NativeAnat.nii');
                NumNodes = 220;
                LeftCortex = 100;
                LeftSubcortex = 110; 
                RightCortex = 210;
                RightSubcortex = NumNodes;
            elseif strcmp(Parcellation, 'cust250')
                Folder = '/custom250_NativeAnat';
                FolderName = strcat(S, Folder);
                cd (FolderName);  
                [~, data_parcel]=read('customparc_NativeAnat.nii');
                NumNodes = 530;
                LeftCortex = 250;
                LeftSubcortex = 265; 
                RightCortex = 515;
                RightSubcortex = NumNodes;

            end
            %% load microarray data
            cd ('/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Microarray/');
            if  UsePCAprobes ==1  && UseCUSTprobes == 0
                fprintf('Loading MicroarrayDataPCAS0%u_GeneThr%d_%s%s.mat\n', subject, round(Thr), Side, BrainPart)
                load(sprintf('MicroarrayDataPCAS0%d_GeneThr%d_%s%s.mat', subject, round(Thr), Side, BrainPart));
            elseif UseMaxVarprobes==1 && UseCUSTprobes == 0
                fprintf('Loading MicroarrayMaxVarS0%u_GeneThr%d_%s%s.mat\n', subject,round(Thr), Side, BrainPart)
                load(sprintf('MicroarrayDataMaxVarS0%d_GeneThr%d_%s%s.mat', subject,round(Thr), Side, BrainPart));
            else 
                error('Choose data file with CUST probes');   
            end
            coordsToAssign = SampleInformation.MRIvoxCoordinates;

            cd (FolderName);
            if DistanceThreshold <30
                load(sprintf('CoordsAssignedAll_%s%s.mat', Side, BrainPart));
            end

         %% find coordinates for nonzero elements in parcellation according to side and brain part (will be used to assign microarray samples to)
            % get coordinate values
            if LeftHemisphere == 1 && OnlyCortex == 1
                Text = sprintf('Subject %d LEFT cortex\n', subject);
                fprintf(1,Text)
                [Coordx,Coordy,Coordz] = ind2sub(size(data_parcel),find(data_parcel > 0 & data_parcel <= LeftCortex));
            elseif LeftHemisphere == 1 && OnlyCortex == 0
                Text = sprintf('Subject %d LEFT subcortex\n', subject);
                fprintf(1,Text)
                [Coordx,Coordy,Coordz] = ind2sub(size(data_parcel),find(data_parcel > LeftCortex & data_parcel <= LeftSubcortex));
            elseif LeftHemisphere == 0 && OnlyCortex == 1
                Text = sprintf('Subject %d RIGHT cortex\n', subject);
                fprintf(1,Text)
                [Coordx,Coordy,Coordz] = ind2sub(size(data_parcel),find(data_parcel > LeftSubcortex & data_parcel <= RightCortex));
            elseif LeftHemisphere == 0 && OnlyCortex == 0
                Text = sprintf('Subject %d RIGHT subcortex\n', subject);
                fprintf(1,Text)
                [Coordx,Coordy,Coordz] = ind2sub(size(data_parcel),find(data_parcel > RightCortex & data_parcel <= RightSubcortex));
            end

            coordsNonzeroParcel = cat(2,Coordx,Coordy,Coordz);

        %% for each microarray coordinate find a closest coordinate in parcellation
            coordsAssigned =  zeros(size(coordsToAssign));
            % finds closest point
            % T = delaunayn(unique(coordsNonzeroParcel));
            k = dsearchn(coordsNonzeroParcel,coordsToAssign);

                for i = 1:length(k)
                    coordsAssigned(i,:) = coordsNonzeroParcel(k(i),:);
                end
                
            %% calculate the distance between original and reassigned coordinate
            assignDistance = zeros(length(coordsAssigned),1);
            coordsNonassigned = zeros(length(coordsAssigned),3);
            coordsToAssignAll = coordsToAssign;
            
                for i=1:length(coordsAssigned)

                    assignDistance(i,1) = pdist2(coordsAssigned(i,:), coordsToAssign(i,:));

                    if assignDistance(i)>DistanceThreshold
                    coordsNonassigned(i,:,:,:)=coordsToAssign(i,:,:,:);
                    coordsAssigned(i,:)=0;
                    coordsToAssign(i,:)=0;

                    end

                end
            %remove zero elements
            coordsToAssign( ~any(coordsToAssign,2), : ) = NaN; 
            coordsAssigned( ~any(coordsAssigned,2), : ) = NaN; 
            coordsNonassigned( ~any(coordsNonassigned,2), : ) = []; 

            %% extract intensity values for assigned coordinates

            cd (FolderName);
            if DistanceThreshold == 30
                coordsAssignedALL = coordsAssigned; 
                save(sprintf('CoordsAssignedAll_%s%s.mat', Side, BrainPart), 'coordsAssignedALL');
            end

            if DistanceThreshold <30
                
                Int=nonzeros(data_parcel);
                int = unique(Int);

                Intensity_all=zeros(size(coordsAssignedALL,1), 1);

                for i=1:size(coordsAssignedALL,1)

                    Intensity = data_parcel(coordsAssignedALL(i,1), coordsAssignedALL(i,2), coordsAssignedALL(i,3));
                    Intensity_all(i) = Intensity;

                end

                %% sort nonzero coordinates according to intensity value 
                Full_informationALL = [coordsAssigned Intensity_all];
                Full_information = Full_informationALL;
                Full_informationALL(any(isnan(Full_informationALL),2),:) = NaN;
                Full_information(any(isnan(Full_informationALL),2),:) = [];
                Full_information_sorted = sortrows(Full_information,4);

                SampleIntensity = Full_informationALL(:,4);
                IntANDExpression = [SampleIntensity Expression];
                % exclude nonassigned points
                IntANDExpression(any(isnan(IntANDExpression),2),:) = NaN;
                IntANDExpression(any(isnan(IntANDExpression),2),:) = [];

                fprintf(1,'Sorts expression data according to ROIs\n')

                % sort according to intensity (roi)
                IntANDExpressionSorted = sortrows(IntANDExpression,1);

                
                %% save the output

                save(sprintf('%d_DistThresh%d_S0%d_MRIvoxCoordsAssigned_%s%s.mat', NumNodes, DistanceThreshold, subject, Side, BrainPart), ... 
                    'coordsAssigned', 'coordsToAssign', 'coordsToAssignAll', 'coordsNonassigned', 'Full_informationALL', 'Full_information_sorted');
                if UseMaxVarprobes == 1;
                        save (sprintf('%d_DistThresh%d_S0%d_ExpressionProbeMaxVar_GeneThr%d_%s%s.mat', NumNodes, DistanceThreshold, subject, round(Thr), Side, BrainPart), 'IntANDExpressionSorted', 'Full_information_sorted' );
                        %savefig (L, sprintf('%d_DistThresh%d_S0%d_ExpressionProbeMaxVar.fig', NumNodes, Threshold, subject));
                elseif UsePCAprobes == 1;
                        save (sprintf('%d_DistThresh%d_S0%d_ExpressionProbePCA_GeneThr%d_%s%s.mat', NumNodes, DistanceThreshold, subject, round(Thr),  Side, BrainPart), 'IntANDExpressionSorted', 'Full_information_sorted' );
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
    
end
