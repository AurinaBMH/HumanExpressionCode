% choose single vertex per sample based on the max value
cd ('data/genes/forFreesurfer')
[verticesFSaverage,facesFSaverage] = read_surf('lhfsaverage.pial');
listChosenVert = cell(6,1); 

for sub=1:6
    fprintf('Processing subject %d\n', sub); 
    cd (sprintf('S0%d', sub))
    samples = load(sprintf('S%dsampleList.txt', sub));
    subjChosenVert = zeros(length(samples),1);
    k=1;
    for i=1:length(samples)
    sample = samples(i);
    dataOrig = MRIread(sprintf('S%dsample%dONfsaveragePial.mgz', sub, sample));
    [~,chosenVert] = max(dataOrig.vol);
    subjChosenVert(k) = chosenVert;
    k=k+1;
    end
    listChosenVert{sub} = subjChosenVert; 
    cd ..
end

% make a combined list of vertices
vertices = vertcat(listChosenVert{1}, listChosenVert{2},listChosenVert{3}, listChosenVert{4}, listChosenVert{5}, listChosenVert{6}); 
coordinates = verticesFSaverage(vertices,:); 

% calculate distances between each pair of vertices
