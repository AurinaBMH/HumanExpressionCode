% choose single vertex per sample

% load a list of samples
samples = load('S2sampleList.txt'); 
listChosenVert = zeros(length(samples),1); 
k=1; 
for i = 1:length(samples)
    sample = samples(i); 
    dataOrig = MRIread(sprintf('S2sample%dONfsaveragePial.mgz', sample));
    [~,chosenVert] = max(dataOrig.vol); 
    listChosenVert(k) = chosenVert; 
    k=k+1; 
end