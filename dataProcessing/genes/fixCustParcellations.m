% script to fix sucrom parcellations
% it uses aparc+aseg parcellation as
subjects = 1:6;
parcellation = 250;
parcellationName = sprintf('customparc%d_NativeAnat.nii', parcellation);
cd ('data/genes/parcellations')
for subj=subjects
    folderName = sprintf('S0%d_H0351',subj);
    cd (folderName);
    [hdr, data_cust]=read(parcellationName);
    [~, data_aparc]=read('defaultparc_NativeAnat.nii');
    
    % mask the cyst parc with aparc parcellation to get rid of the borders
    data_custFixed = data_cust.*logical(data_aparc);
    
    % loop through each voxel
    for i=2:size(data_cust,1)-1
        for j=2:size(data_cust,2)-1
            for z=2:size(data_cust,3)-1
                % evaluate each voxel in mask (aparc+aseg) and cust
                % parcellation
                parcel = data_custFixed(i,j,z);
                mask = data_aparc(i,j,z);
                
                % if a voxel is 0 in cust and non-zero in aparc, then fix it.
                if parcel==0 && mask~=0
                    
                    % get the surroundings of thet voxel
                    s = data_custFixed(i-1:i+1, j-1:j+1, z-1:z+1);
                    % if the surroundings are not 0, do fix it
                    if ~isempty(nonzeros(s(:)))
                        % calculate the percentage of each non-zero intensity value
                        t = tabulate(nonzeros(s(:)));
                        % count then number of max values
                        rep = find(t(:,2) == max(t(:,2)));
                        % if 2 or more intensity values are found to be equally
                        % popular, increase the radius of searching and do that until one
                        % intensity value 'wins'
                        k=1;
                        while length(rep)>1
                            k=k+1;
                            s = data_custFixed(i-k:i+k, j-k:j+k, z-k:z+k);
                            t = tabulate(nonzeros(s(:)));
                            rep = find(t(:,2) == max(t(:,2)));
                        end
                        % assign thet most popular values to the empty voxel
                        [~,ind] = max(t(:,2));
                        val = t(ind,1);
                        data_custFixed(i,j,z) = val;
                    end
                    
                end
                
            end
        end
        
    end
    write(hdr,data_custFixed,sprintf('customparc%d_NativeAnatFixed.nii', parcellation));
    cd ..
end
cd ../../..
