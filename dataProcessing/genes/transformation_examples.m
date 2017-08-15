% Just a set of scripts to get the code started
cd ('data/genes/processedData')
load('samples2keep.mat')
load('MicroarrayDatadPC82DistThresh2_CoordsAssigned.mat'); 
cd ../
for sub = 1:6
    cd (sprintf('forFreesurfer/S0%d', sub))
% Here first load up the freesurfer volume and lh pial surface
data = MRIread('orig.mgz');
[vertices,faces] = read_surf('lh.white');

vert2 = [vertices-1 ones(size(vertices,1),1)];
h2 = inv(data.tkrvox2ras)*vert2.';

ab = find(abs(h2(3,:)-120)<0.5);
%figure;
%imagesc(data.vol(:,:,120));
%hold on;
%plot(h2(1,ab),h2(2,ab),'*');

% Now to have the other files in the same space -- there is an extra
% transformation step -- have to check hdr information.

%%
data2 = MRIread('001.nii');
%h3 = inv(data2.vox2ras)*vert2.';
h3 = inv(data.tkrvox2ras*inv(data.vox2ras)*(data2.niftihdr.sform))*vert2.';

ab = find(abs(h3(3,:)-120)<0.5);
%figure;
%imagesc(data2.vol(:,:,120));
%hold on;
%plot(h3(1,ab),h3(2,ab),'*');

%% Vertices in voxel space for samples and T1
newVertices = h3(1:3,:).';
[xx, yy, zz] = meshgrid(1:size(data2.vol,2),1:size(data2.vol,1),1:size(data2.vol,3));
%figure;slice(xx,yy,zz,data2.vol,[120],[120],[120]);shading interp;
%hold on;
%patch('Vertices',newVertices,'faces',faces+1,'FaceColor','none');

%% 

% now a position on newVertices corresponds to a position in voxel space,
% so it is in voxel co-ordinates.
% load overlay for a subject
dataOrig = MRIread(sprintf('S%dsamples.mgz', sub));
dataOrig.vol = zeros(size(dataOrig.vol)); 

%keepSamples = sampleIND{sub}; 
Lcortex = 1:34; 
coordinatesall = DataCoordinatesMRI{sub}; 
ROI = coordinatesall(:,2); keep = find(ROI<=34); 
coordinates = coordinatesall(keep,3:5); 

coordinatesNEWvox = zeros(length(keep),3); 
coordinatesNEWvert = zeros(length(keep),3); 
overlay = zeros(size(vertices)); 
for i=1:length(keep)
x = coordinates(i,1); 
y = coordinates(i,2);
z = coordinates(i,3);

[minval, ind] =min(sqrt((newVertices(:,1)-x).^2 + (newVertices(:,2)-y).^2 + (newVertices(:,3)-z).^2));
vertexind = ind;

% In voxel space vertex coordinate is:
%disp([newVertices(vertexind,1),newVertices(vertexind,2),newVertices(vertexind,3)]);
 % save coordinates
coordinatesNEWvox(i,:) = [newVertices(vertexind,1),newVertices(vertexind,2),newVertices(vertexind,3)];
coordinatesNEWvert(i,:) = [vertices(vertexind,1),vertices(vertexind,2),vertices(vertexind,3)];

overlay(vertexind) = i; 
dataOrig.vol(vertexind) = 1; 
% In freesurfecoordinatesNEWvertr space vertex co-ordinate is:
%disp([vertices(vertexind,1),vertices(vertexind,2),vertices(vertexind,3)]);
end
MRIwrite(dataOrig,sprintf('S%dsamples_mod.mgz', sub)); 
end
% find index of each sample in vertex space based on the coordinates
%write_surf('S1samples_vertex.mgz', coordinatesNEWvert, faces); 


% Use this information to check if its working correctly.