
clear all; 
cd ('data/genes/processedData')
load('MicroarrayDataProbesUpdated.mat')
% select genes that have multiple probes, so thay will be sub-selected for
% comparison
[v, ind] = unique(DataTableProbe.EntrezID{1});
% m=0; w=1; 
% for p=1:length(ind)
%     A = find(DataTableProbe.EntrezID{1}==v(p)); 
%     %howMany = length(A); 
%     if length(A)>1
%         m=m+1;
%         PL{p} = A;
%         for k=1:length(A)
%         IND2(w) = A(k);
%         w=w+1;
%         end
%         %w=w+length(A);
%     end
%     
% end

% duplicate_ind = IND2; 
% duplicate_value = unique(DataTableProbe.EntrezID{1}(duplicate_ind));

duplicate_ind = setdiff(1:size(DataTableProbe.EntrezID{1}, 1), ind);
duplicate_value = unique(DataTableProbe.EntrezID{1}(duplicate_ind));

percentage = (length(duplicate_value)/length(unique(DataTableProbe.EntrezID{1})))*100; 
format compact
%fprintf('%d genes have more than one probe\n', m); 
percentage 

% Load probes, selected using different methods
load('MicroarrayDataProbesUpdatedMean.mat')
probes{1} = probeInformation; 
expression{1} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6}); 

load('MicroarrayDataProbesUpdatedLessNoise.mat')
probes{2} = probeInformation; 
expression{2} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6}); 

load('MicroarrayDataProbesUpdatedPC.mat')
probes{3} = probeInformation; 
expression{3} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6}); 

load('MicroarrayDataProbesUpdatedDS.mat')
probes{4} = probeInformation; 
expression{4} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6}); 

%load('MicroarrayDataProbesUpdatedRNAseq3.000000e-01thr.mat')
load('MicroarrayDataProbesUpdatedRNAseq2.000000e-01RNAthr5.000000e-01noisethr.mat')
probes{5} = probeInformation; 
expression{5} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6}); 

load('MicroarrayDataProbesUpdatedRandom.mat')
probes{6} = probeInformation; 
expression{6} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6}); 

load('MicroarrayDataProbesUpdatedRandom2.mat')
probes{7} = probeInformation; 
expression{7} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6}); 

load('MicroarrayDataProbesUpdatedVariance.mat')
probes{8} = probeInformation; 
expression{8} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6}); 

load('MicroarrayDataProbesUpdatedCV.mat')
probes{9} = probeInformation; 
expression{9} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6}); 


%select only genes that had more than one probe available
[genesMultiple, indFilter] = intersect(probes{1}.EntrezID, duplicate_value); % separatelly for RNAseq and others as the numbef of genes is different
[genesMultipleRNAseq, indFilterRNA] = intersect(probes{5}.EntrezID, duplicate_value); 

% filter genes that had multiple probes
for k=1:9
    if k==5
        expression{k} = expression{k}(:,indFilterRNA); 
    else 
        expression{k} = expression{k}(:, indFilter); 
    end
end

[v1, indALL] = intersect(genesMultiple, genesMultipleRNAseq); 
[v2, indRNA] = intersect(genesMultipleRNAseq,genesMultiple); 
% calculate correlation between each way of choosing a probe
avCorr = zeros(9,9); 

for i=1:9
%     if i==5
%         expr1 = expression{i}(:,indRNA);
%     else
%         expr1 = expression{i};
%     end
    
    for j=1:9
        if j==5 && i~=5
            expr1 = expression{i}(:,indALL);
            expr2 = expression{j}(:,indRNA);
            
        elseif j~=5 && i==5
            expr1 = expression{i}(:,indRNA);
            expr2 = expression{j}(:,indALL);
        elseif j~=5 && i~=5
            expr1 = expression{i}(:,indALL);
            expr2 = expression{j}(:,indALL);
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

% reorder according to similarity
R = BF_pdist(avCorr); 
[ord,R,keepers] = BF_ClusterReorder(avCorr,R); 
avCorrPlot = avCorr(ord, ord); 

nice_cmap = [make_cmap('steelblue',50,30,0);flipud(make_cmap('orangered',50,30,0))];

figure; imagesc(avCorrPlot);set(gcf,'color','w'); 
colormap(nice_cmap)
caxis([0.5 1])
tickNames = {'Mean', 'Noise', 'PC', 'DS', 'RNAseq', 'Random1', 'Random2','Variance', 'CV'}; 
tickNames = tickNames(ord); 

xticks([1 2 3 4 5 6 7 8 9])
xticklabels(tickNames); 
yticks([1 2 3 4 5 6 7 8 9])
yticklabels(tickNames); 
