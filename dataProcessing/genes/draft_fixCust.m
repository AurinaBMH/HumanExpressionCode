        [~, data_parcel]=read('HCPMMP1_acpc_uncorr.nii');
        [~, data_parcelmask]=read('defaultparc_NativeAnat.nii');
        
        for i=
            for j=
                for z=
                    
                    
                end
            end
        end
        
        A = find(data_parcel==0); 
        B = find(data_parcelmask~=0); 
        C = intersect(A,B); 