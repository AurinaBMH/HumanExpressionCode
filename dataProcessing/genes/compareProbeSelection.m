% load initial data
clear all; 
cd ('data/genes/processedData')
load('MicroarrayDataProbesUpdated.mat')
% select genes that have multiple probes, so thay will be sub-selected for
% comparison
[~, ind] = unique(DataTableProbe.EntrezID{1});
% duplicate indices
duplicate_ind = setdiff(1:size(DataTableProbe.EntrezID{1}, 1), ind);
duplicate_value = DataTableProbe.EntrezID{1}(duplicate_ind);

% Load probes, selected using different methods
load('MicroarrayDataProbesUpdatedVariance.mat')
probes{1} = probeInformation; 
load('MicroarrayDataProbesUpdatedRandom.mat')
probes{2} = probeInformation; 
load('MicroarrayDataProbesUpdatedLessNoise.mat')
probes{3} = probeInformation; 
load('MicroarrayDataProbesUpdatedPC.mat')
probes{4} = probeInformation; 
load('MicroarrayDataProbesUpdatedRNAseq.mat')
probes{5} = probeInformation; 
load('MicroarrayDataProbesUpdatedRandom2.mat')
probes{6} = probeInformation; 

%select only genes that had more than one probe available
[~, indFilter] = intersect(probes{1}.EntrezID, duplicate_value); % separatelly for RNAseq and others as the numbef of genes is different
[~, indFilterRNA] = intersect(probes{5}.EntrezID, duplicate_value); 

% select only the ones that are left
probes{1}.ProbeID = probes{1}.ProbeID(indFilter); probes{1}.EntrezID = probes{1}.EntrezID(indFilter); 
probes{2}.ProbeID = probes{2}.ProbeID(indFilter); probes{2}.EntrezID = probes{2}.EntrezID(indFilter);  
probes{3}.ProbeID = probes{3}.ProbeID(indFilter); probes{3}.EntrezID = probes{3}.EntrezID(indFilter); 
probes{4}.ProbeID = probes{4}.ProbeID(indFilter); probes{4}.EntrezID = probes{4}.EntrezID(indFilter);   
probes{6}.ProbeID = probes{6}.ProbeID(indFilter); probes{6}.EntrezID = probes{6}.EntrezID(indFilter); 
probes{5}.ProbeID = probes{5}.ProbeID(indFilterRNA); probes{5}.EntrezID = probes{5}.EntrezID(indFilterRNA); 

% find genes that are both in RNAseq and in t=other methods
[both, indRNA, indOthers] = intersect(probes{5}.EntrezID, probes{1}.EntrezID); 

 perc = zeros(6,6);
 % calculate what percentage of probes match
 for j=1:6
     if j==5
         probe1 = probes{j}.ProbeID(indRNA);
     else
         probe1 = probes{j}.ProbeID(indOthers);
     end
     for i=1:6
         if i==5
             probe2 = probes{i}.ProbeID(indRNA);
         else
             probe2 = probes{i}.ProbeID(indOthers);
         end
         
         perc(i,j) = (length(intersect(probe1, probe2)))/length(probe2);
         
     end
 end
 % plot
nice_cmap = [make_cmap('steelblue',50,30,0);flipud(make_cmap('orangered',50,30,0))];

figure; imagesc(perc);set(gcf,'color','w'); 
colormap(nice_cmap)
caxis([0 1])
xticks([1 2 3 4 5 6])
xticklabels({'Variance', 'Random', 'Noise', 'PC','RNAseq', 'Random2'}); 
yticks([1 2 3 4 5 6])
yticklabels({'Variance', 'Random', 'Noise', 'PC','RNAseq', 'Random2'}); 
