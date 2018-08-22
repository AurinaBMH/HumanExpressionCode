
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

i=1;
for weight = weights
    j=1;
    for tract = tracts
        [coexpData, A, matrices, coordinates, avWeight] = giveConnCoexp(parc,tract{1},weight{1},brainPart);
         [~, Alength, matricesLength] = giveConnCoexp(parc,tract{1},'length',brainPart);
        
%         if strcmp(weight{1},'standard')
%             % for standard weight (where streamline count is used), divide wegith
%             % by the length of the streamline
% 
%             for s=1:length(matrices)
%                 matrices{s} = log(matrices{s}); %./matricesLength{s};
%                 matrices{s}(isnan(matrices{s})) = 0;
%                 matrices{s}(~isfinite(matrices{s}))=0; % all log values after log transformation are>2.3, because we first exclude links
%                 % with weight<10. So 0 value will not be assigned to any
%                 % exiting edges
%                 A(:,:,s) = matrices{s};
%             end
%         end

        [~, ~, ~, ~, avWeightLength] = giveConnCoexp(parc,tract{1},'length',brainPart);
        
        giveMeRichClub(matrices, coordinates, groupConn ,densThreshold, giveRC);
        title(sprintf('RC and distance curves, %s %s, %.2f density', tract{1}, weight{1}, densThreshold))
        [Gr] = giveMeRichClub(matrices, coordinates, groupConn ,densThreshold, false);
        
        binGr = logical(Gr);
        
        GroupMatrB{i,j} = binGr;
        GroupMatrW{i,j} = Gr;
        
        % calculate properties
        % degree distributions
        nodeDeg{i,j} = degrees_und(Gr);
        nodeStrength{i,j} = strengths_und(Gr);
        figure; histogram(nodeStrength{i,j}, 50); title(sprintf('%s, %s', tract{1}, weight{1}));
        xlabel('Connection strength')
        conMask = double(binGr);
        conMask(conMask==0) = NaN;
        % connection length distributions
        GroupMatrL{i,j} = maskuHalf(conMask.*avWeightLength);
        
        figure; histogram(GroupMatrL{i,j}, 50); title(sprintf('%s, %s', tract{1}, weight{1}));
        xlabel('Connection distance')
        % global efficiency
        Eglob(i,j) = efficiency_bin(binGr);
        
        j=j+1;
    end
    i=i+1;
end

figure; imagesc(Eglob); 

i=1;
for weight = weights

    A = GroupMatrB{i,1}-GroupMatrB{i,2}; 
    figure; imagesc(A); title (sprintf('FACT - iFOD (%s)', weight{1}))

    i=i+1;
end

FACTlog = log(GroupMatrW{2,1}); figure; histogram(FACTlog); xlim([0 12]); xlabel('log weight'); title('FACT')
iFODlog = log(GroupMatrW{2,2}); figure; histogram(iFODlog); xlim([0 12]); xlabel('log weight'); title('iFOD')


