% count how many orobes have 0 mismatches. 
probes = probes2annotateALLmergedreadAnnotation(2:end,:); 
nrMismatches = cell2mat(probes(:,4)); 
nr0 = length(find(nrMismatches==0)); 

% perc 0 mismatches. 

nrProbes = length(nrMismatches); 
prop0Mismatch = nr0/nrProbes; 