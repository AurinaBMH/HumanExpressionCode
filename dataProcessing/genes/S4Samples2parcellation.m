%% Author: Aurina
%Last modiffied: 2016-04-14
%Last modiffied: 2017-07-20
%Last modiffied: 2017-07-31

% It assigns samples to parcellation separately for cortex/subcortex,
% right/left sides using only relevant parcellation voxels to make an assignment.

% closest points are chosen using dsearchn function, which assignes a parcellation coordinate index to each microarray sample coordinate.
% output: original microarray coordinates to be assigned (coordsToAssign)
% and new assigned coordinates (coordsAssigned) and expression data are saved in
% each subject's folder separately for each parcellation and left/right cortex/subcortex.


%------------------------------------------------------------------------------
% Choose options
%------------------------------------------------------------------------------
% choose if you want to use data with CUST probes
clear all;
useCUSTprobes = false;
% choose what type of probe selection to use, hemisphere, subject list, parcellations, threshols.
probeSelections = {'RNAseq'}; % {'Mean', 'Variance', 'LessNoise'};
parcellations = {'cust100'}; %,'cust100', 'cust250', 'aparcaseg', HCP};
distanceThreshold = 2; % first run with 40 for one probeSelection (will make a file that fits all)
% then run with 2 for all probe selections. 
subjects = 1:6;

%------------------------------------------------------------------------------
% Create variables to save data in
%------------------------------------------------------------------------------
DataExpression = cell(length(subjects),1);
DataCoordinatesMRI = cell(length(subjects),1);
DataCoordinatesMNI = cell(length(subjects),1);
coordsAssignedALL = cell(length(subjects),1);
%------------------------------------------------------------------------------
% Select variables according to side/brain part selections
%------------------------------------------------------------------------------
sides = {'left', 'right'};

%------------------------------------------------------------------------------
% Do assignment for all subjects
%------------------------------------------------------------------------------
for parcellation = parcellations
    for t=1:length(probeSelections)
        assignDistance = cell(length(subjects),1); 
        for subject = subjects
            cd ('data/genes/parcellations')
            subjectDir = sprintf('S0%d_H0351', subject);
            cd (subjectDir)
            
            if strcmp(parcellation, 'HCP')
                brainParts = {'Cortex'}; %, 'Subcortex'};
            else
                brainParts = {'Cortex','Subcortex'};
            end
            fprintf('Subject %u parcellation %s assignment distance threshold %u\n; ', subject, parcellation{1}, distanceThreshold )
            
            %------------------------------------------------------------------------------
            % Load parcellations
            %------------------------------------------------------------------------------
            if strcmp(parcellation, 'aparcaseg')
                parcName = 'default_NativeAnat';
                [~, data_parcel]=read('defaultparc_NativeAnat.nii');
                NumNodes = 82;
                LeftCortex = 1:34;
                LeftSubcortex = 35:41;
                RightCortex = 42:75;
                RightSubcortex = 76:82;
            elseif strcmp(parcellation, 'cust100')
                parcName = 'custom100_NativeAnat';
                [~, data_parcel]=read('customparc100_NativeAnatFixed.nii');
                NumNodes = 220;
                LeftCortex = 1:100;
                LeftSubcortex = 101:110;
                RightCortex = 111:210;
                RightSubcortex = 211:220;
            elseif strcmp(parcellation, 'cust250')
                parcName = 'custom250_NativeAnat';
                [~, data_parcel]=read('customparc250_NativeAnatFixed.nii');
                NumNodes = 530;
                LeftCortex = 1:250;
                LeftSubcortex = 251:265;
                RightCortex = 266:515;
                RightSubcortex = 516:530;
            elseif strcmp(parcellation, 'HCP')
                parcName = 'HCP';
                [~, data_parcel]=read('HCPMMP1_acpc_uncorr.nii');
                NumNodes = 360;
                LeftCortex = 1:180;
                %LeftSubcortex = 35:41;
                RightCortex = 181:360;
                %RightSubcortex = 76:82;
                
            end
            cd ../../
            
            %------------------------------------------------------------------------------
            % Load microarray data
            %------------------------------------------------------------------------------
            cd ('processedData');
            
            if useCUSTprobes
                fprintf(1,'Using the data with CUST probes %s for %d\n', probeSelections{t}, subject)
                startFileName = 'MicroarrayDataWITHcust';
            else
                fprintf(1,'Using the data without CUST probes %s for %d\n', probeSelections{t}, subject)
                startFileName = 'MicroarrayDataProbesUpdated';
            end
            
            
            %load(sprintf('ProbeInformation%s.mat', probeSelection));
            load(sprintf('%s%ssep.mat', startFileName, probeSelections{t}));
            
            
            for side = sides
                for brainPart = brainParts
                    
                    coords2assign = sampleInfoSubjects{subject}.(side{1}).(brainPart{1}).MRIvoxCoordinates;
                    
                    if distanceThreshold < 35
                        
                        load(sprintf('CoordsAssignedAll%d.mat', NumNodes));
                        
                    end
                    %------------------------------------------------------------------------------
                    % Find coordinates for nonzero elements in parcellation according to side and brain part (will be used to assign microarray samples to)
                    %------------------------------------------------------------------------------
                    % get coordinate values
                    if strcmp(side, 'left') && strcmp(brainPart, 'Cortex')
                        Text = sprintf('Subject %d LEFT cortex\n', subject);
                        fprintf(1,Text)
                        [Coordx,Coordy,Coordz] = ind2sub(size(data_parcel),find(ismember(data_parcel, LeftCortex)));
                    elseif strcmp(side, 'left') && strcmp(brainPart, 'Subcortex')
                        Text = sprintf('Subject %d LEFT subcortex\n', subject);
                        fprintf(1,Text)
                        [Coordx,Coordy,Coordz] = ind2sub(size(data_parcel),find(ismember(data_parcel, LeftSubcortex)));
                    elseif strcmp(side, 'right') && strcmp(brainPart, 'Cortex')
                        Text = sprintf('Subject %d RIGHT cortex\n', subject);
                        fprintf(1,Text)
                        [Coordx,Coordy,Coordz] = ind2sub(size(data_parcel),find(ismember(data_parcel, RightCortex)));
                    elseif strcmp(side, 'right') && strcmp(brainPart, 'Subcortex')
                        Text = sprintf('Subject %d RIGHT subcortex\n', subject);
                        fprintf(1,Text)
                        [Coordx,Coordy,Coordz] = ind2sub(size(data_parcel),find(ismember(data_parcel, RightSubcortex)));
                    end
                    
                    coordsNonzeroParcel = cat(2,Coordx,Coordy,Coordz);
                    
                    %------------------------------------------------------------------------------
                    % For each microarray coordinate find a closest coordinate in parcellation
                    %------------------------------------------------------------------------------
                    
                    coordsAssigned =  zeros(size(coords2assign));
                    % find closest point
                    k = dsearchn(coordsNonzeroParcel,coords2assign);
                    
                    for i = 1:size(k,1)
                        coordsAssigned(i,:) = coordsNonzeroParcel(k(i),:);
                    end
                    
                    %------------------------------------------------------------------------------
                    % Salculate the distance between original and reassigned coordinate
                    %------------------------------------------------------------------------------
                    %(size(coordsAssigned,1),1);
                    coordsNONassigned = zeros(size(coordsAssigned,1),3);
                    coordsToAssignAll = coords2assign;
                    
                   % assignDistance{subject}.(side{1}).(brainPart{1}) = zeros(size(coordsAssigned,1),1); 
                    
                    for j=1:size(coordsAssigned,1)
                        
                        assignDistance{subject}.(side{1}).(brainPart{1})(j,1) = pdist2(coordsAssigned(j,:), coords2assign(j,:));
                        
                        if assignDistance{subject}.(side{1}).(brainPart{1})(j,1)>distanceThreshold
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
                    
                    if distanceThreshold >35
                        
                        coordsAssignedALL{subject}.(side{1}).(brainPart{1}) = coordsAssigned;
                        
                    elseif distanceThreshold <35
                        
                        Int=nonzeros(data_parcel);
                        int = unique(Int);
                        
                        intensity_all=zeros(size(coordsAssignedALL{subject}.(side{1}).(brainPart{1}),1), 1);
                        
                        for i=1:size(coordsAssignedALL{subject}.(side{1}).(brainPart{1}),1)
                            
                            intensity = data_parcel(coordsAssignedALL{subject}.(side{1}).(brainPart{1})(i,1), coordsAssignedALL{subject}.(side{1}).(brainPart{1})(i,2), coordsAssignedALL{subject}.(side{1}).(brainPart{1})(i,3));
                            intensity_all(i) = intensity;
                            
                        end
                        indRemove = find(any(isnan(coordsAssigned),2));
                        intensity_all(indRemove, :) = [];
                        
                        MNIcoordinates = sampleInfoSubjects{subject}.(side{1}).(brainPart{1}).MMCoordinates;
                        
                        MNIcoordinates(indRemove, :) = [];
                        coordsAssigned(indRemove, :) = [];
                        
                        %------------------------------------------------------------------------------
                        % sort nonzero coordinates according to intensity value
                        %------------------------------------------------------------------------------
                        informationALL = [intensity_all coordsAssigned];
                        
                        
                        informationALLmni = [intensity_all MNIcoordinates];
                        information_all = informationALL;
                        informationmni = informationALLmni;
                        
                        informationALL(any(isnan(informationALL),2),:) = NaN;
                        information_all(any(isnan(informationALL),2),:) = [];
                        
                        informationALLmni(any(isnan(informationALL),2),:) = NaN;
                        informationmni(any(isnan(informationALL),2),:) = [];
                        
                        information_sorted = sortrows(information_all,1);
                        information_sortedmni = sortrows(informationmni,1);
                        
                        
                        sampleIntensity = informationALL(:,1);
                        expr = expressionSubjects{subject}.(side{1}).(brainPart{1});
                        expr(indRemove,:)=[];
                        IntANDExpression = [sampleIntensity expr];
                        % exclude nonassigned points
                        IntANDExpression(any(isnan(IntANDExpression),2),:) = NaN;
                        IntANDExpression(any(isnan(IntANDExpression),2),:) = [];
                        
                        fprintf(1,'Sorts expression data according to ROIs\n')
                        
                        % sort according to intensity (roi)
                        IntANDExpressionSorted = sortrows(IntANDExpression,1);
                        % assign outputs to structure
                        information.(side{1}).(brainPart{1}).raw = information_all;
                        information.(side{1}).(brainPart{1}).sorted = information_sorted;
                        coordinates.(side{1}).(brainPart{1}).assigned = coordsAssigned;
                        coordinates.(side{1}).(brainPart{1}).toassign = coords2assign;
                        coordinates.(side{1}).(brainPart{1}).toassignALL = coordsToAssignAll;
                        coordinates.(side{1}).(brainPart{1}).nonassigned = coordsNONassigned;
                        
                        data.(side{1}).(brainPart{1}).expression = IntANDExpressionSorted;
                        data.(side{1}).(brainPart{1}).informationMRI = information_sorted;
                        data.(side{1}).(brainPart{1}).informationMNI = information_sortedmni;
                        
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
            
            
            %------------------------------------------------------------------------------
            % Save output
            %------------------------------------------------------------------------------
            if distanceThreshold <35
                if strcmp(parcellation, 'HCP')
                    nSamples = size(data.left.Cortex.informationMRI,1)+size(data.right.Cortex.informationMRI,1);
                    Expression = cat(1,data.left.Cortex.expression, data.right.Cortex.expression);
                    CoordinatesMRI = cat(1,data.left.Cortex.informationMRI, data.right.Cortex.informationMRI);
                    CoordinatesMNI = cat(1,data.left.Cortex.informationMNI, data.right.Cortex.informationMNI);
                else
                    nSamples = size(data.left.Cortex.informationMRI,1)+size(data.left.Subcortex.informationMRI,1)+size(data.right.Cortex.informationMRI,1)+size(data.right.Subcortex.informationMRI,1);
                    Expression = cat(1,data.left.Cortex.expression,data.left.Subcortex.expression, data.right.Cortex.expression, data.right.Subcortex.expression);
                    CoordinatesMRI = cat(1,data.left.Cortex.informationMRI,data.left.Subcortex.informationMRI, data.right.Cortex.informationMRI, data.right.Subcortex.informationMRI);
                    CoordinatesMNI = cat(1,data.left.Cortex.informationMNI,data.left.Subcortex.informationMNI, data.right.Cortex.informationMNI, data.right.Subcortex.informationMNI);
                end
                SUBJECT = zeros(nSamples,1);
                SUBJECT(:,1) = subject;
                
                DataExpression{subject} = [SUBJECT, Expression];
                DataCoordinatesMRI{subject} = [SUBJECT, CoordinatesMRI];
                DataCoordinatesMNI{subject} = [SUBJECT, CoordinatesMNI];
                
                %save(sprintf('%s%s%dDistThresh%d_CoordsAssigned_S0%d.mat', startFileName, probeSelection, NumNodes, distanceThreshold, subject), ...
                %'data');
                
            elseif distanceThreshold > 35
                
                save(sprintf('CoordsAssignedAll%d.mat', NumNodes), 'coordsAssignedALL');
                
            end
            % cd ../../..
            cd ../../..
        end
        
        
        %% save data for all subjects
        cd ('data/genes/processedData')
        if distanceThreshold < 35
            save(sprintf('%s%s%dDistThresh%d_CoordsAssigned.mat', startFileName, probeSelections{t}, NumNodes, distanceThreshold), 'DataExpression', 'DataCoordinatesMRI', 'DataCoordinatesMNI', 'probeInformation', 'assignDistance');
        end
        cd ../../..
        
    end
end



