[~, data_parcel]=read('defaultparc_NativeAnat.nii');

data = (0<data_parcel) & (data_parcel<=34); 
ind = 1:(size(data,1)*size(data,2)*size(data,3)); 
%data2 = ind; 
dataNew = reshape(ind,size(data)); 
nodes = length(nonzeros(data(:)));
arround = zeros(nodes,3,3,3);
arroundDist = zeros(nodes,3,3,3);
nodeList = zeros(nodes,3); 


n=1;
a=1;
i=1; 
c=sqrt(a^2+a^2);
surrValues1 = [c,a,c; a,0,a; c,a,c];
surrValues23 = [c,c,c; c,a,c; c,c,c];

surrValues3D(1,:,:) = surrValues23;
surrValues3D(2,:,:) = surrValues1; 
surrValues3D(3,:,:) = surrValues23;
for x=1:size(data,1)
    for y=1:size(data,2)
        for z=1:size(data,3)
            node = data(x,y,z);
            % for every nonzero node in the parcellation
            if node~=0
                % get edge weights based on distance between nodes    
                values = surrValues3D.*data(x-1:x+1,y-1:y+1,z-1:z+1);
                % select indexes for each neighbour for a list
                neighbours = dataNew(x-1:x+1,y-1:y+1,z-1:z+1); 
                N = neighbours(:); 
                % assign values for each nonzero node
                nodeList(i:i+26,1) = dataNew(x,y,z); % index of the first node
                nodeList(i:i+26,2) = N;              % index of the second node
                nodeList(i:i+26,3) = values(:);      % weight between them
                i=i+27; 
            end
        end
    end
end


