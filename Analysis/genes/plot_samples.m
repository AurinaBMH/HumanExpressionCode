
clear all; close all; 
cd ('/Users/Aurina/GoogleDrive/Genetics_connectome/HumanExpression/data/genes/forFreesurfer')
[vertices,faces] = read_surf(['lhfsaverage.pial']);
figure;p = patch('Vertices',vertices,'faces',faces+1);
[vertex,label,ctab] = read_annotation('lh.aparc.annot');
overlay = zeros(size(vertex));
for i=1:36
    inds = find(label==ctab.table(i,5));
    overlay(inds) = i*1000;
    % for j=1:length(inds)
    % overlay(inds(j),:) = ctab.table(i,1:3);
    % end
    
end

for sub=1:6
fprintf('Processing subject %d\n', sub);
cd (sprintf('S0%d', sub))
samples = load(sprintf('S%dsampleList.txt', sub));
subjChosenVert = zeros(length(samples),1);
k=1;
for i=1:length(samples)
    sample = samples(i);
    
    dataOrig = MRIread(sprintf('S%dsample%dONfsaveragePial.mgz',sub, sample));
   
    [~,chosenVert] = max(dataOrig.vol);
    subjChosenVert(k) = chosenVert;
    k=k+1;
end
listChosenVert{sub} = subjChosenVert;
cd ..
end

cd (sprintf('S0%d', 1))

%[verticesFSaverage,facesFSaverage] = read_surf('lh.pial');
% make a combined list of vertices
INDS = vertcat(listChosenVert{1}, listChosenVert{2},listChosenVert{3}, listChosenVert{4}, listChosenVert{5}, listChosenVert{6}); 
overlay(INDS) = 50000; 


set(p,'FaceVertexCData',overlay,'FaceColor','interp','EdgeColor','none');
nice_cmap = [make_cmap('steelblue',50,30,0);flipud(make_cmap('orangered',50,30,0))];
set(gcf,'color','w');
colormap(nice_cmap)
axis image
axis off
camlight('left')
material dull
