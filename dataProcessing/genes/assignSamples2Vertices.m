% Just a set of scripts to get the code started
cd ('data/genes/processedData')
load('MicroarrayDatadPC82DistThresh2_CoordsAssigned.mat');
cd ../forFreesurfer
keepSamplesOrig = cell(6,1);
keepSamples = cell(6,1);

[verticesFSaverage,facesFSaverage] = read_surf('lhfsaverage.white');

k=1;
for sub = 1:6
    cd (sprintf('S0%d', sub))
    % Here first load up the freesurfer volume and lh pial surface
    data = MRIread('orig.mgz');
    [vertices,faces] = read_surf('lh.pial');
    
    vert2 = [vertices-1 ones(size(vertices,1),1)];
    h2 = inv(data.tkrvox2ras)*vert2.';
    
    data2 = MRIread('001.nii');
    h3 = inv(data.tkrvox2ras*inv(data.vox2ras)*(data2.niftihdr.sform))*vert2.';
    
    % Vertices in voxel space for samples and T1
    newVertices = h3(1:3,:).';
    [xx, yy, zz] = meshgrid(1:size(data2.vol,2),1:size(data2.vol,1),1:size(data2.vol,3));
    
    % load overlay for a subject
    dataOrig = MRIread(sprintf('S%dsamples.mgz', sub));
    
    coordinatesall = DataCoordinatesMRI{sub};
    ROI = coordinatesall(:,2); keep = find(ROI<=34);
    coordinates = coordinatesall(keep,3:5);
    
    coordinatesNEWvox = zeros(length(keep),3);
    coordinatesNEWvert = zeros(length(keep),3);
    overlay = zeros(size(vertices));
    %int=linspace(1,1249,1249);
    
    for i=1:length(keep)
        dataOrig.vol = zeros(size(dataOrig.vol));
    
        x = coordinates(i,1);
        y = coordinates(i,2);
        z = coordinates(i,3);
        
        [minval, ind] =min(sqrt((newVertices(:,1)-x).^2 + (newVertices(:,2)-y).^2 + (newVertices(:,3)-z).^2));
        vertexind = ind;
        
        % save coordinates
        coordinatesNEWvox(i,:) = [newVertices(vertexind,1),newVertices(vertexind,2),newVertices(vertexind,3)];
        coordinatesNEWvert(i,:) = [vertices(vertexind,1),vertices(vertexind,2),vertices(vertexind,3)];
        
        lab = i+k-1; 
        overlay(vertexind) = i;
        dataOrig.vol(vertexind) = i;
        
        MRIwrite(dataOrig,sprintf('S%dsample%d_singleVert.mgz', sub, i));
        % In freesurfecoordinatesNEWvertr space vertex co-ordinate is:
        % disp([vertices(vertexind,1),vertices(vertexind,2),vertices(vertexind,3)]);
    end
    
    [~,ia] = unique(coordinatesNEWvert, 'rows', 'stable');
    keepSamplesOrig{sub} = keep;
    k=k+length(ia);
    keepSamples{sub} = ia;
    sampleList = keepSamplesOrig{sub}; 
    
    fileID = fopen(sprintf('S%dsampleList.txt', sub),'w');
    nbytes = fprintf(fileID,'%1d\n',sampleList); 
    fclose(fileID);
    
    
    %MRIwrite(dataOrig,sprintf('S%dsamples_singleVert.mgz', sub));
    cd ..
    
end
save('keepSamples.mat', 'keepSamples');
