
% plot the proportion of assigned samples as a function of distanc
% threshold
nrSamples(:,1) = [393; 671; 714; 731; 739];
nrSamples(:,2) = [319; 538; 604; 619; 622];
nrSamples(:,3) = [156; 255; 281; 288; 295];
nrSamples(:,4) = [216; 350; 383; 397; 401];
nrSamples(:,5) = [179; 281; 318; 326; 329];
nrSamples(:,6) = [187; 310; 350; 358; 362];
nrSamples(:,7) = sum(nrSamples,2);

percSamples = nrSamples./nrSamples(5,:); 

figure;
imagesc(percSamples([1:4],:)); 
title('Percentage of assigned samples')
%colormap([;BF_getcmap('reds',9)]);
caxis([0 1]); 
xticks([1:7])
xticklabels({'S01','S02','S03','S04','S05','S06','all'})
xlabel('Subjects'); ylabel('Distance threshold')
yticks([1:4])
yticklabels({'0mm','2mm','5mm','10mm'})
colorbar; 
set(gcf,'color','w');
nice_cmap = [make_cmap('steelblue',50,30,0);flipud(make_cmap('orangered',50,30,0))];
colormap(nice_cmap)
caxis([0.5 1])