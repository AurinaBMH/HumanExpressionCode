% count numbers of samples in CX, BS, CB

cd ('data/genes/rawData');
for subj=1:6
    fprintf(1,'Loading data for %u subject\n', subj)
    folder = sprintf('normalized_microarray_donor0%d', subj);
    cd (folder);
    %%load information specific for each subject

    FileAnnot = 'SampleAnnot.xlsx';
                       % exclude probe IDs from expression matrix
    [~,~,SlabType] = xlsread(FileAnnot, 'D:D');

    SlabType(1) = [];                           % remove headline

    
    % To exclude expression data and coordinates for braintem (BS) and cerebellum (CB)
    % exclude columns in expression and rows in coordinates if slabtype is CB or BS

        fprintf('Excluding brainstem and cerebellum data\n')
        
        CX = strfind(SlabType, 'CX'); CXBSCBind(subj,1) = length(find(~cellfun(@isempty,CX)));
        BS = strfind(SlabType, 'BS'); CXBSCBind(subj,2) = length(find(~cellfun(@isempty,BS)));
        CB = strfind(SlabType, 'CB'); CXBSCBind(subj,3) = length(find(~cellfun(@isempty,CB)));
cd ..
end