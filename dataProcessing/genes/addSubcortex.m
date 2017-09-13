parc = 'cust250';

switch parc
    case 'cust100'
        NumNodes=220;
        LeftCortex = 1:100;
        LeftSubcortex = 101:110;
        RightCortex = 111:210;
        RightSubcortex = 211:220;
        num1 = 100; 
        num2 = 200; 
        n=10;
    case 'cust250'
        
        NumNodes = 530;
        LeftCortex = 1:250;
        LeftSubcortex = 251:265;
        RightCortex = 266:515;
        RightSubcortex = 516:530;
        num1=250; 
        num2=500; 
        n=15; 
end
[hdr, data_parcel_subcortical] = read(sprintf('customparc%d_NativeAnatFixedOLD.nii', num1));

remove = [LeftCortex,RightCortex];
data_parcel_subcortical(ismember(data_parcel_subcortical, remove))=0;

[~,data_parcel_cortical] = read(sprintf('custom%d_acpc_uncorr.nii', num2));
% replace
data_parcel_cortical= single(data_parcel_cortical);
data_parcel_subcortical= single(data_parcel_subcortical);
% change values for cortical ROIs
RC = find(data_parcel_cortical>max(LeftCortex));
data_parcel_cortical(RC) = data_parcel_cortical(RC)+n;

% combina cortex+subcortex
data = data_parcel_cortical+data_parcel_subcortical;
mask = logical(data_parcel_cortical) + logical(data_parcel_subcortical);

% find overlap and assign those voxels to cortex (cortical parcellation is
% more accurate)
A = find(mask==2);
data(A)=data_parcel_cortical(A);
write(hdr, data, sprintf('customparc%d_NativeAnatFixedTEST.nii', num1))

