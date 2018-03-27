
clear all; close all; 
cd ('/Users/Aurina/GoogleDrive/Genetics_connectome/HumanExpression/data/genes/forFreesurfer')
[vertices,faces] = read_surf(['lhfsaverage.pial']);

figure;p = patch('Vertices',vertices,'faces',faces+1);
[vertex,label,ctab] = read_annotation('lh.aparc.annot');
overlay = zeros(size(vertex));
inds = find(label==ctab.table(19,5));
overlay(inds) = 100; 

set(p,'FaceVertexCData',overlay,'FaceColor','interp','EdgeColor','none');
nice_cmap = [make_cmap('steelblue',50,30,0);flipud(make_cmap('orangered',50,30,0))];
set(gcf,'color','w');
colormap(nice_cmap)
axis image
axis off
camlight('left')
material dull

% count samples in the region
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

A1 = length(intersect(inds, listChosenVert{1}))
A2 = length(intersect(inds, listChosenVert{2}))
A3 = length(intersect(inds, listChosenVert{3}))
A4 = length(intersect(inds, listChosenVert{4}))
A5 = length(intersect(inds, listChosenVert{5}))
A6 = length(intersect(inds, listChosenVert{6}))
