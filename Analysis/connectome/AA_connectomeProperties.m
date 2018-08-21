
%------------------------------------------------------------------------------
% Compare FACT and iFOD2 connectomes at the same density with the samemethod
%------------------------------------------------------------------------------

parc = 'HCP';
tracts = {'FACT','iFOD2'};
weights = {'density', 'standard', 'FA'};
brainPart = 'wholeBrain';
groupConn = 'variance'; %  lengthCV, consistency
densThreshold = 0.1;
giveRC = true;

 
figure; 

i=1;
for weight = weights
    j=1;
    for tract = tracts
        [coexpData, A, matrices, coordinates, avWeight] = giveConnCoexp(parc,tract{1},weight{1},brainPart);
        [~, ~, ~, ~, avWeightLength] = giveConnCoexp(parc,tract{1},'length',brainPart);
        
        giveMeRichClub(matrices, coordinates, groupConn ,densThreshold, giveRC);
        [Gr] = giveMeRichClub(matrices, coordinates, groupConn ,densThreshold, false);
        title(sprintf('RC and distance curves, %s %s, %.2f density', tract{1}, weight{1}, densThreshold))
        binGr = logical(Gr);
        
        GroupMatrB{i,j} = binGr;
        GroupMatrW{i,j} = Gr;
        
        % calculate properties
        % degree distributions
        nodeDeg{i,j} = degrees_und(Gr);
        conMask = double(binGr);
        conMask(conMask==0) = NaN;
        % connection length distributions
        GroupMatrL{i,j} = maskuHalf(conMask.*avWeightLength);

        figure; histogram(GroupMatrL{i,j}, 50); title(sprintf('%s, %s', tract{1}, weight{1})); 
        % global efficiency
        Eglob(i,j) = efficiency_bin(binGr);
        
        j=j+1;
    end
    i=i+1;
end

figure; imagesc(Eglob); 





parc = 'HCP';

tract = 'iFOD2';
weights = {'density', 'standard', 'FA'};
brainPart = 'wholeBrain';
groupConn = 'variance'; %  lengthCV, consistency
densThreshold = 0.1;
giveRC = false;

i=1;
for weight = weights
    [coexpData, A, matrices, coordinates, avWeight] = giveConnCoexp(parc,tract,weight{1},brainPart);
    [Gr] = giveMeRichClub(matrices, coordinates, groupConn ,densThreshold, giveRC);
    %title(sprintf('RC and distance curves, %.2f density', densThreshold))
    nodeDeg{i} = degrees_und(Gr);
    GroupMatr{i} = logical(Gr);
    i=i+1;
end