function groupAdj = giveMeGroupAdj(connectomes)

M = zeros(size(connectomes{1},1),size(connectomes{1},2),size(connectomes,2)); 
nSubs = size(connectomes,2); 
d = zeros(nSubs,1); 
for i=1:nSubs
    M(:,:,i) = connectomes{i}; 
    d(i) = density_und(connectomes{i}); 
end
dMean = mean(d); 
groupAdj = threshold_consistency(M, dMean); 

end