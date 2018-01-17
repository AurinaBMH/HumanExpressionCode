
clear all; 
cd ('data/genes/processedData')
load('MicroarrayDataProbesUpdated.mat')
% select genes that have multiple probes, so thay will be sub-selected for
% comparison
[v, ind] = unique(DataTableProbe.EntrezID{1});
m=0; 

for p=1:length(ind)
    A = find(DataTableProbe.EntrezID{1}==v(p)); 
    %howMany = length(A); 
    if length(A)>1
        m=m+1;
    end
    
end
percentage = (m/length(unique(DataTableProbe.EntrezID{1})))*100; 
format compact
fprintf('%d genes have more than one probe\n', m); 
percentage 


[~, ind] = unique(DataTableProbe.EntrezID{1});
% duplicate indices
duplicate_ind = setdiff(1:size(DataTableProbe.EntrezID{1}, 1), ind);
duplicate_value = DataTableProbe.EntrezID{1}(duplicate_ind);

% Load probes, selected using different methods
load('MicroarrayDataProbesUpdatedMean.mat')
probes{1} = probeInformation; 
expression{1} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6}); 

load('MicroarrayDataProbesUpdatedVariance.mat')
probes{2} = probeInformation; 
expression{2} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6}); 

load('MicroarrayDataProbesUpdatedRandom.mat')
probes{6} = probeInformation; 
expression{6} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6}); 

load('MicroarrayDataProbesUpdatedLessNoise.mat')
probes{3} = probeInformation; 
expression{3} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6}); 

load('MicroarrayDataProbesUpdatedPC.mat')
probes{4} = probeInformation; 
expression{4} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6}); 

load('MicroarrayDataProbesUpdatedRNAseq5.000000e-01thr.mat')
probes{5} = probeInformation; 
expression{5} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6}); 

load('MicroarrayDataProbesUpdatedRandom2.mat')
probes{7} = probeInformation; 
expression{7} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6}); 

%select only genes that had more than one probe available
[genesMultiple, indFilter] = intersect(probes{1}.EntrezID, duplicate_value); % separatelly for RNAseq and others as the numbef of genes is different
[genesMultipleRNAseq, indFilterRNA] = intersect(probes{5}.EntrezID, duplicate_value); 

% filter genes that had multiple probes
for k=1:7
    if k==5
        expression{k} = expression{k}(:,indFilterRNA); 
    else 
        expression{k} = expression{k}(:, indFilter); 
    end
end

[v1, indALL] = intersect(genesMultiple, genesMultipleRNAseq); 
[v2, indRNA] = intersect(genesMultipleRNAseq,genesMultiple); 
% calculate correlation between each way of choosing a probe
avCorr = zeros(7,7); 

for i=1:7
%     if i==5
%         expr1 = expression{i}(:,indRNA);
%     else
%         expr1 = expression{i};
%     end
    
    for j=1:7
        if j==5 && i~=5
            expr1 = expression{i}(:,indALL);
            expr2 = expression{j}(:,indRNA);
            
        elseif j~=5 && i==5
            expr1 = expression{i}(:,indRNA);
            expr2 = expression{j}(:,indALL);
        elseif j~=5 && i~=5
            expr1 = expression{i};
            expr2 = expression{j};
        else
            expr1 = expression{i}(:,indRNA);
            expr2 = expression{j}(:,indRNA);
        end
        
correlation = zeros(size(expr1,2),1); 

        for g=1:size(expr1,2)
            correlation(g) = corr(expr1(:,g), expr2(:,g), 'type', 'Spearman');
        end
        
        
        avCorr(j,i) = mean(correlation);
    end
    
end

nice_cmap = [make_cmap('steelblue',50,30,0);flipud(make_cmap('orangered',50,30,0))];

figure; imagesc(avCorr);set(gcf,'color','w'); 
colormap(nice_cmap)
caxis([0.5 1])
xticks([1 2 3 4 5 6 7])
xticklabels({'Mean', 'Variance', 'Noise', 'PC','RNAseq', 'Random1', 'Random2'}); 
yticks([1 2 3 4 5 6 7])
yticklabels({'Mean', 'Variance', 'Noise', 'PC','RNAseq', 'Random1', 'Random2'}); 
