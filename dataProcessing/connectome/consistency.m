% check if consistency relates to volumes
% give grpup connectome and variation matrix
study = 'HCP'; % (HCP)
dens = 1; %0.5; 

switch study
    case 'HCP'
        
        load('HCPMMP1_acpc_SIFT2_FACT_volume.mat')
        load('HCPMMP1ANDfslatlas20_acpc_connectome_data.mat')
    case 'GenCog'
        
        
        load('GenCog_HCPMMP1ANDfslatlas20_default_11_volume.mat')
        load('HCPMMP1ANDfslatlas20_GenCOG_connectome_data.mat')
end

load('HCPMMP1ANDfslatlas20_MNILinear_COGflippedX.mat')

type = standard;

[groupAdj, consist] = giveMeGroupAdj(type, dens); 

if strcmp(study, 'HCP')
groupAdj = groupAdj([1:180,191:370],[1:180,191:370]);
consist = consist([1:180,191:370],[1:180,191:370]);
end



weightVariation = consist.*logical(groupAdj);
volumes = vol.*logical(groupAdj);

weightVariation(isnan(weightVariation))=0; 
data(:,1) = weightVariation(:); data(:,2)=volumes(:); 
data( ~any(data,2), : ) = []; 
figure; scatter(data(:,1), data(:,2));
xlabel('variation in weights'); ylabel('combined volume of 2 regions'); 
title(sprintf('%s', study)); 


% get distnces between rois
% give grpup connectome and variation matrix
%[groupLength, consistLength] = giveMeGroupAdj(SIFT2_length, 0.1);
[groupAdj, consist] = giveMeGroupAdj(type, dens); 
weightVariation = consist.*logical(groupAdj);

length = pdist2(coordinates, coordinates); %groupLength.*logical(groupAdj); 
dataLength(:,1) = weightVariation(:); dataLength(:,2)=length(:);
dataLength(any(isnan(dataLength), 2), :) = [];
dataLength = dataLength(all(dataLength,2),:);

dataLength( ~any(dataLength,2), : ) = []; 
figure; scatter(dataLength(:,1), dataLength(:,2));
xlabel('variation in weights'); ylabel('distance between regions'); 
title(sprintf('%s', study)); 