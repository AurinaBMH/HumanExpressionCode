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
parcellations = {'aparcaseg'}; %,'cust100', 'cust250', 'aparcaseg', HCP};
distanceThreshold = 35; % first run with 40 for one probeSelection (will make a file that fits all)
% then run with 2 for all probe selections.
cd ('data/genes/processedData')
load('MicroarrayDataProbesUpdatedRNAseq.mat')
[~, ~, parcelLabels] = xlsread('82parcelLabels.xlsx','Sheet1');
parcelLabels(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),parcelLabels)) = {''};
subjects = 1:6;
out = cell(6,1);
samplesLR = cell(6,1);
distLR = cell(6,1);
assignDistance = cell(6,1); 
cd ../../..

for subject = subjects
    cd ('data/genes/parcellations')
    subjectDir = sprintf('S0%d_H0351', subject);
    cd (subjectDir)
    
    parcName = 'default_NativeAnat';
    [~, data_parcel]=read('defaultparc_NativeAnat.nii');
    NumNodes = 82;
    LeftCortex = 1:34;
    LeftSubcortex = 35:41;
    RightCortex = 42:75;
    RightSubcortex = 76:82;
    
    coords2assign = sampleInfo{subject}.MRIvoxCoordinates;
    
    [Coordx,Coordy,Coordz] = ind2sub(size(data_parcel),find(ismember(data_parcel, 1:82)));
    coordsNonzeroParcel = cat(2,Coordx,Coordy,Coordz);
    
    coordsAssigned =  zeros(size(coords2assign));
    % find closest point
    k = dsearchn(coordsNonzeroParcel,coords2assign);
    
    for i = 1:size(k,1)
        coordsAssigned(i,:) = coordsNonzeroParcel(k(i),:);
    end
    
    Int=nonzeros(data_parcel);
    int = unique(Int);
    
    intensity_all=zeros(size(coordsAssigned, 1),1);
    
    for i=1:size(coordsAssigned,1)
        
        intensity = data_parcel(coordsAssigned(i,1),...
            coordsAssigned(i,2), ...
            coordsAssigned(i,3));
        intensity_all(i) = intensity;
        
    end
    
    for j=1:size(coordsAssigned,1)
        
        assignDistance{subject}(j,1) = pdist2(coordsAssigned(j,:), coords2assign(j,:));
        
    end
    
    indL = find(intensity_all<=41);
    indR = find(intensity_all>=42);
    
    % gest structure names for these selected
    samplesLR{subject,1} = sampleInfo{subject}.StructureNames(indL);
    distLR{subject,1} = assignDistance{subject}(indL);
    samplesLR{subject,2} = sampleInfo{subject}.StructureNames(indR);
    distLR{subject,2} = assignDistance{subject}(indR);
    
    % find samples with assignment distance >2 and compare their labels and
    % parcellation ID
    
    indCheck = find(assignDistance{subject}>2);
    % get their percel
    parcelLabel = intensity_all(indCheck);
    % structure name
    sampleLabel = sampleInfo{subject}.StructureNames(indCheck);
    A = [sampleLabel,num2cell(parcelLabel)];
    out{subject} = A; 
    
    for i=1:length(sampleLabel)
        out{subject}{i,3} = parcelLabels(A{i,2});
    end
    cd ../../../..
    
end



for p=1:2
    m=0;
    if p==1
        side = 'right';
    elseif p==2
        side = 'left';
    end
    for i=1:6
        % sheck how many right samples are assigned to left
        for j=1:size(samplesLR{i,p},1)
            if strfind(samplesLR{i,p}{j}, side)
                m=m+1;
            end
        end
        
    end
    total(p) = m;
end






