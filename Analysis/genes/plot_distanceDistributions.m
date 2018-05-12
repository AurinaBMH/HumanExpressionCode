cd ('data/genes/processedData');
load('DistancesONsurfaceXXX.mat'); 
surf = maskuHalf(distSamples);
surf(isnan(surf)) = []; 
load('distancesGM_MNIXXX.mat'); 
GM = maskuHalf(distSamples); 
GM(isnan(GM)) = []; 
load('distancesEuclideanXXX.mat')
eucl = maskuHalf(distSamples); 
eucl(isnan(eucl)) = []; 

figure; 
set(gcf,'color','w'); 
histogram(surf, 50, 'EdgeColor', [.45 .45 .45], 'FaceColor', [.75 .5 .51]); hold on;
histogram(GM, 50, 'EdgeColor', [.45 .45 .45], 'FaceColor', [.85 .54 .4]); hold on; 
histogram(eucl, 50, 'EdgeColor', [.45 .45 .45], 'FaceColor',[.43 .6 .47] ); hold on; 
legend({'On the surface', 'Within grey matter volume', 'Euclidean'})
box('off')
xlabel('Distance between samples (mm)'); 
ylabel('Number of sample pairs')
set(gca,'fontsize',15)

