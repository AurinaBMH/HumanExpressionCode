%------------------------------------------------------------------------------
% Choose options
%------------------------------------------------------------------------------
parcellation = {'aparcaseg'};%, 'cust100', 'cust250'};
probeSelection = 'Mean';% (Variance', LessNoise', 'Mean')
distanceThreshold = 2; % first run 30, then with the final threshold 2
subjects = 1;

if strcmp(parcellation, 'aparcaseg')
    NumNodes = 82;
elseif strcmp(parcellation, 'cust100')
    NumNodes = 220;
elseif strcmp(parcellation, 'cust250')
    NumNodes = 530;
end

cd ('data/genes/processedData')
load(sprintf('MicroarrayDatad%s%dDistThresh%d_CoordsAssigned.mat', probeSelection, NumNodes, distanceThreshold));

%------------------------------------------------------------------------------
% Select variables according to side/brain part selections
%------------------------------------------------------------------------------
sides = {'left', 'right'};
brainParts = {'Cortex', 'Subcortex'};


for subject = subjects
    %------------------------------------------------------------------------------
    % Load microarray samples for distance calculations
    %------------------------------------------------------------------------------
    LcortexIND = find(DataCoordinates{subject}(:,2)<=34);
    samples = DataCoordinates{subject}(LcortexIND,3:5);
    
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
        [~, data_parcel]=read('defaultparc_NativeAnat.nii');
        NumNodes = 82;
        LeftCortex = 1:34;
        LeftSubcortex = 35:41;
        RightCortex = 42:75;
        RightSubcortex = 76:82;
    elseif strcmp(parcellation, 'cust100')
        parcName = 'custom100_NativeAnat';
        cd (parcName);
        [~, data_parcel]=read('customparc_NativeAnat.nii');
        NumNodes = 220;
        LeftCortex = 1:100;
        LeftSubcortex = 101:110;
        RightCortex = 111:210;
        RightSubcortex = 211:220;
    elseif strcmp(parcellation, 'cust250')
        parcName = 'custom250_NativeAnat';
        cd (parcName);
        [~, data_parcel]=read('customparc_NativeAnat.nii');
        NumNodes = 530;
        LeftCortex = 1:250;
        LeftSubcortex = 251:265;
        RightCortex = 266:515;
        RightSubcortex = 516:530;
        
    end
    
    %[~, data_parcel]=read('defaultparc_NativeAnat.nii');
    
    data = (0<data_parcel) & (data_parcel<=length(LeftCortex));
    ind = 1:(size(data,1)*size(data,2)*size(data,3));
    
    dataNew = reshape(ind,size(data));
    nodes = length(nonzeros(data(:)));
    arround = zeros(nodes,3,3,3);
    arroundDist = zeros(nodes,3,3,3);
    nodeList = zeros(nodes*3,3);
    
    
    n=1;
    a=1;
    i=1;
    c=sqrt(a^2+a^2);
    % assign distance values to the selected plain
    surrValues2 = [c,a,c; a,0,a; c,a,c];
    % assign distance values to the more distant planes
    surrValues13 = [c,c,c; c,a,c; c,c,c];
    
    % make a 3D cube of distances
    surrValues3D(1,:,:) = surrValues13;
    surrValues3D(2,:,:) = surrValues2;
    surrValues3D(3,:,:) = surrValues13;
    
    % loop for each voxel in the parcellation
    for x=1:size(data,1)
        for y=1:size(data,2)
            for z=1:size(data,3)
                node = data(x,y,z);
                % for every nonzero node in the parcellation
                if node~=0
                    % get edge weights based on distance between nodes
                    weights = surrValues3D.*data(x-1:x+1,y-1:y+1,z-1:z+1);
                    % select indexes for each neighbour for a list
                    neighbourIDs = dataNew(x-1:x+1,y-1:y+1,z-1:z+1);
                    listIDs = neighbourIDs(:);
                    % assign values for each nonzero node
                    nodeList(i:i+26,1) = dataNew(x,y,z);        % index of the first node
                    nodeList(i:i+26,2) = listIDs;               % index of the second node
                    nodeList(i:i+26,3) = weights(:);            % weight between them
                    i=i+27;
                end
            end
        end
        fprintf('X loop nr %d out of %d\n', x, size(data,1))
    end
    % remove edges with zero weight
    indices = find(nodeList(:,3)==0);
    nodeList(indices,:) = [];
    
    % remove duplicate edges
    
    A = nodeList(:,1:2);
    B = [A(:,2),A(:,1)];
    C = [A;B];
    W = [nodeList(:,3);nodeList(:,3)];
    
    [tf, loc]=ismember(C(:,[2 1]),C,'rows');
    C = C(loc>=(1:size(loc,1))',:);
    W = W(loc>=(1:size(loc,1))',:);
    nodeListUnique(:,1:2) = C;
    nodeListUnique(:,3) = W;

    
    % make a graph from a list of edges
    G = graph(nodeListUnique(:,1),nodeListUnique(:,2),nodeListUnique(:,3));
    
    % get index values for samples
    indSamples = zeros(size(samples,1),1);
    for j=1:size(samples,1)
        indSamples(j) = dataNew(samples(j,1), samples(j,2), samples(j,3));
    end
    
    % calculate distance between each pair of samples
    shortDist = zeros(size(samples,1));
    for sam1 = 1:size(samples,1)
        for sam2=sam1+1:size(samples,1)
            
            [~,shortDist(sam1,sam2)] = shortestpath(G,indSamples(sam1),indSamples(sam2));
            
        end
    end
    
    
    
    
    
    
    % for each pair of
    
end
cd ../../../


