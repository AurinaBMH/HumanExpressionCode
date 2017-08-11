% for each subject make a .nii file with all zeros and assign indexes for
% the coordinates of mappped samples.
subjects = 1; 
brainPart = 'Leftcortex'; 
numSamples = 1; 
for subject = subjects
    cd ('data/genes/parcellations')
    subjectDir = sprintf('S0%d_H0351', subject);
    cd (subjectDir)
    %for parcellation = parcellations
    
    fprintf('Subject %u parcellation %s assignment distance threshold %u\n; ', subject, parcellation{1}, distanceThreshold )
    
    %------------------------------------------------------------------------------
    % Load parcellations
    %------------------------------------------------------------------------------
    if strcmp(parcellation, 'aparcaseg')
        parcName = 'default_NativeAnat';
        cd (parcName);
        [hdr, data_parcel]=read('defaultparc_NativeAnat.nii');
        NumNodes = 82;
        LeftCortex = 1:34;
        LeftSubcortex = 35:41;
        RightCortex = 42:75;
        RightSubcortex = 76:82;
    elseif strcmp(parcellation, 'cust100')
        parcName = 'custom100_NativeAnat';
        cd (parcName);
        [hdr, data_parcel]=read('customparc_NativeAnat.nii');
        NumNodes = 220;
        LeftCortex = 1:100;
        LeftSubcortex = 101:110;
        RightCortex = 111:210;
        RightSubcortex = 211:220;
    elseif strcmp(parcellation, 'cust250')
        parcName = 'custom250_NativeAnat';
        cd (parcName);
        [hdr, data_parcel]=read('customparc_NativeAnat.nii');
        NumNodes = 530;
        LeftCortex = 1:250;
        LeftSubcortex = 251:265;
        RightCortex = 266:515;
        RightSubcortex = 516:530;
        
    end
    % make an empty image
    image = zeros(size(data_parcel)); 
    coordinates = DataCoordinatesMRI{subject}(:,3:5);
    switch brainPart
        case 'Lcortex'
            nROIs = LeftCortex;
        case 'LcortexSubcortex'
            nROIs = 1:max(LeftSubcortex);
        case 'wholeBrain'
            nROIs = 1:NumNodes;
        case 'LRcortex'
            nROIs = [LeftCortex,RightCortex];
    end
    samplesIND = find(ismember(DataCoordinatesMRI{subject}(:,2),nROIs));
    index = numSamples:numSamples+max(samplesIND)-1; 
    numSamples = numSamples+length(index); 
    
    for samp=1:length(index)
        point = coordinates(samp,:); 
        image(point(1), point(2), point(3)) = index(samp); 
    end
    cd ../../..
    cd ('processedData')
    filename = sprintf('S%dsamples.nii', subject); 
    write(hdr,image,filename); 
    cd ../../..
   
end

