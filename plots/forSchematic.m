load('MicroarrayDataProbesUpdated.mat')

vectexp = Expressionall(:,1723)'; 
vectbin = noiseall(:,1723)'; 

nice_cmap = [make_cmap('steelblue',50,30,0);flipud(make_cmap('orangered',50,30,0))];

figure; imagesc(vect); 
colormap(nice_cmap); 
caxis([5 12]); 
figure; imagesc(vectbin); 
colormap(gray)


ind = find(DataTableProbe.EntrezID{1,1}==740);
E = Expressionall(ind,:); 
figure; box off; imagesc(E); colormap(nice_cmap); 