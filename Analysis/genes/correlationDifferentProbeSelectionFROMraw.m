% make a figure for probe correlation from scratch

%1. regenerate data for different probe selection methods
options.signalThreshold = -1; % no QC filtering
options.useCUSTprobes = true;

% files will be saved with standart name noQC
Regeneratedata = false;
if Regeneratedata
    options.probeSelections = {'Mean'};
    S2_probes(options);
    
    options.probeSelections = {'maxIntensity'};
    S2_probes(options);
    
    options.probeSelections = {'LessNoise'};
    S2_probes(options);
    
    options.probeSelections = {'PC'};
    S2_probes(options);
    
    options.probeSelections = {'DS'};
    S2_probes(options);
    
    options.probeSelections = {'Random'};
    S2_probes(options)
    % rename into Random2
    
    options.probeSelections = {'Random'};
    S2_probes(options)
    
    options.probeSelections = {'Variance'};
    S2_probes(options)
    
    options.probeSelections = {'CV'};
    S2_probes(options)
end
% load the initial data
% created with there options
% options.ExcludeCBandBS =  true;
% options.useCUSTprobes = true;
% options.updateProbes = 'reannotator';

cd ('data/genes/processedData')
load('MicroarrayDataWITHcustProbesUpdatedXXX.mat')
[v, ind] = unique(DataTableProbe.EntrezID{1});

duplicate_ind = setdiff(1:size(DataTableProbe.EntrezID{1}, 1), ind);
duplicate_value = unique(DataTableProbe.EntrezID{1}(duplicate_ind));

percentage = (length(duplicate_value)/length(unique(DataTableProbe.EntrezID{1})))*100;
format compact
percentage


% Load probes, selected using different methods that were just generated
load('MicroarrayDataWITHcustProbesUpdatedXXXMeannoQC.mat')
probes{1} = probeInformation;
expression{1} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

load('MicroarrayDataWITHcustProbesUpdatedXXXLessNoisenoQC.mat')
probes{2} = probeInformation;
expression{2} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

load('MicroarrayDataWITHcustProbesUpdatedXXXPCnoQC.mat')
probes{3} = probeInformation;
expression{3} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

load('MicroarrayDataWITHcustProbesUpdatedXXXDSnoQC.mat')
probes{4} = probeInformation;
expression{4} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

load('MicroarrayDataWITHcustProbesUpdatedXXXRandomnoQC.mat')
probes{5} = probeInformation;
expression{5} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

load('MicroarrayDataWITHcustProbesUpdatedXXXRandom2noQC.mat')
probes{6} = probeInformation;
expression{6} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

load('MicroarrayDataWITHcustProbesUpdatedXXXVariancenoQC.mat')
probes{7} = probeInformation;
expression{7} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

load('MicroarrayDataWITHcustProbesUpdatedXXXCVnoQC.mat')
probes{8} = probeInformation;
expression{8} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

load('MicroarrayDataWITHcustProbesUpdatedXXXmaxIntensitynoQC.mat')
probes{9} = probeInformation;
expression{9} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

%select only genes that had more than one probe available
[genesMultiple, indFilter] = intersect(probes{1}.EntrezID, duplicate_value); % separatelly for RNAseq and others as the numbef of genes is different

% filter genes that had multiple probes
for k=1:length(probes)
    
    expression{k} = expression{k}(:, indFilter);
    
end

avCorr = zeros(9,9);
stdCorr = zeros(9,9); 

for i=1:9
    
    
    for j=i+1:9
        
        expr1 = expression{i};
        expr2 = expression{j};
        
        correlation = zeros(size(expr1,2),1);
        
        for g=1:size(expr1,2)
            correlation(g) = corr(expr1(:,g), expr2(:,g), 'type', 'Spearman');
        end
        
        avCorr(j,i) = mean(correlation);
        stdCorr(i,j) = std(correlation); 
    end
    
end
avCorrfull = avCorr+avCorr'; 
avCorrfull(logical(eye(size(avCorrfull)))) = 1; 

stdCorrfull = stdCorr+stdCorr'; 

% reorder according to similarity
R = BF_pdist(avCorrfull);
[ord,R,keepers] = BF_ClusterReorder(avCorrfull,R);
avCorrPlot = avCorrfull(ord, ord);
stdCorrPlot = stdCorrfull(ord, ord); 

nice_cmap = [flipud(make_cmap('orangered',50,30,0))];

figure; imagesc(avCorrPlot);set(gcf,'color','w');
colormap(flipud(pink))
caxis([0.5 1])
tickNames = {'mean', 'signal proportion', 'PC', 'DS', 'Rand_1', 'Rand_2','Variance', 'CV', 'max intensity'};
tickNamesORD = tickNames(ord);

xticks([1 2 3 4 5 6 7 8 9])
xticklabels(tickNamesORD);
xtickangle(45)
yticks([1 2 3 4 5 6 7 8 9])
yticklabels(tickNamesORD);
ytickangle(45)
set(gca,'FontSize', 14)

cd ../../..
% make a plot with RNA seq
options.probeSelections = {'RNAseq'};
options.RNAseqThreshold = 0.2;
options.signalThreshold = -1; % no QC filtering
options.useCUSTprobes = true;
if Regeneratedata
S2_probes(options)
end

% load data
cd ('data/genes/processedData')
load('MicroarrayDataWITHcustProbesUpdatedXXXRNAseqnoQC.mat')
probes{10} = probeInformation;
expression{10} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});
% select genes in other versions and in this one  that overlap
entrezIDs = probes{1}.EntrezID(indFilter);

[genesMultiple, indFilter] = intersect(entrezIDs, duplicate_value); % separatelly for RNAseq and others as the numbef of genes is different
[genesMultipleRNAseq, indFilterRNA] = intersect(probes{10}.EntrezID, duplicate_value);

% filter genes that had multiple probes
for k=1:10
    if k==10
        expression{k} = expression{k}(:,indFilterRNA);
    else
        expression{k} = expression{k}; %; (:,indFilter);
    end
end

[v1, indALL] = intersect(genesMultiple, genesMultipleRNAseq);
[v2, indRNA] = intersect(genesMultipleRNAseq,genesMultiple);
% calculate correlation between each way of choosing a probe

RNAcorr = zeros(9,1);
stDEV = zeros(9,1); 
for i=10
    for j=1:9
        
        expr1 = expression{i}(:,indRNA);
        expr2 = expression{j}(:,indALL);
        
        
        correlation = zeros(size(expr1,2),1);
        
        for g=1:size(expr1,2)
            correlation(g) = corr(expr1(:,g), expr2(:,g), 'type', 'Spearman');
        end
        
        
        RNAcorr(j) = mean(correlation);
        stDEV(j) = std(correlation); 
    end
end

[valS, indS] = sort(RNAcorr, 'descend');
namesS = tickNames(indS);

nice_cmap = [make_cmap('steelblue',50,30,0);flipud(make_cmap('orangered',50,30,0))];
colors2use = nice_cmap([5 15 25 35 45 55 65 75 85 95],:); 


figure; 
set(gcf,'color','w');
hold on
for i = 1:length(valS)
    h=bar(i,valS(i));
    %errorbar(i,valS(i),stDEV(i),'.')
    if strcmp(namesS(i), 'CV')
        set(h,'FaceColor',[1 0.8 0.6],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
    elseif strcmp(namesS(i), 'Variance')
        set(h,'FaceColor',[1 0.8 0.6],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
    elseif strcmp(namesS(i), 'Rand_1')
        set(h,'FaceColor',[.96 .63 .55],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
    elseif strcmp(namesS(i), 'Rand_2')
        set(h,'FaceColor',[.96 .63 .55],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
    elseif strcmp(namesS(i), 'PC')
        set(h,'FaceColor',[1 0.8 0.6],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
    elseif strcmp(namesS(i), 'mean')
        set(h,'FaceColor',[.95 .6 .6],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
    elseif strcmp(namesS(i), 'signal proportion')
        set(h,'FaceColor',[.72 .43 .47],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
    elseif strcmp(namesS(i), 'DS')
        set(h,'FaceColor',[.72 .43 .47],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);                 
    elseif strcmp(namesS(i), 'max intensity')
        set(h,'FaceColor',[.72 .43 .47],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);                 
    end
end
hold off
ylim([0.5 1]); ylabel('Spearman correlation'); xlabel('Probe selection methods');
%title(sprintf(' (%d)', length(indRNA)))
xticks([1 2 3 4 5 6 7 8 9])
%legend(sprintf('%d genes', length(indRNA)))
xticklabels(namesS);
xtickangle(45)
set(gca,'FontSize', 14)


%% now do the same thing only with QC filtering done first
% make a figure for probe correlation from scratch
% make sure files are saved without "noQC" label at the end
clear all;

Regeneratedata = false;
if Regeneratedata
    %1. regenerate data for different probe selection methods
    
    % files will be saved with standart name noQC
    options.probeSelections = {'Mean'};
    S2_probes(options);
    
    options.probeSelections = {'LessNoise'};
    S2_probes(options);
    
    options.probeSelections = {'PC'};
    S2_probes(options);
    
    options.probeSelections = {'DS'};
    S2_probes(options);
    
    options.probeSelections = {'Random'};
    S2_probes(options)
    % rename into Random2
    
    options.probeSelections = {'Random'};
    S2_probes(options)
    
    options.probeSelections = {'Variance'};
    S2_probes(options)
    
    options.probeSelections = {'CV'};
    S2_probes(options)
end

% load the initial data
% created with there options
% options.ExcludeCBandBS =  true;
% options.useCUSTprobes = true;
% options.updateProbes = 'reannotator';

%cd ('data/genes/processedData')
load('MicroarrayDataWITHcustProbesUpdatedXXX.mat')
options.signalThreshold = 0.5; % no QC filtering
options.useCUSTprobes = true;

signalLevel = sum(noiseall,2)./size(noiseall,2);
indKeepProbes = find(signalLevel>=options.signalThreshold);

% sduplicated values calculated only on those probes that pass QC
% threshold;
[v, ind] = unique(DataTableProbe.EntrezID{1}(indKeepProbes));
duplicate_ind = setdiff(1:size(DataTableProbe.EntrezID{1}(indKeepProbes), 1), ind);
entrezSelected = DataTableProbe.EntrezID{1}(indKeepProbes);
duplicate_value = unique(entrezSelected(duplicate_ind));

percentage = (length(duplicate_value)/length(unique(entrezSelected)))*100;
format compact
percentage

% Load probes, selected using different methods that were just generated
load('MicroarrayDataWITHcustProbesUpdatedXXXMean.mat')
probes{1} = probeInformation;
expression{1} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

load('MicroarrayDataWITHcustProbesUpdatedXXXLessNoise.mat')
probes{2} = probeInformation;
expression{2} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

load('MicroarrayDataWITHcustProbesUpdatedXXXPC.mat')
probes{3} = probeInformation;
expression{3} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

load('MicroarrayDataWITHcustProbesUpdatedXXXDS.mat')
probes{4} = probeInformation;
expression{4} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

load('MicroarrayDataWITHcustProbesUpdatedXXXRandom.mat')
probes{5} = probeInformation;
expression{5} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

load('MicroarrayDataWITHcustProbesUpdatedXXXRandom2.mat')
probes{6} = probeInformation;
expression{6} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

load('MicroarrayDataWITHcustProbesUpdatedXXXVariance.mat')
probes{7} = probeInformation;
expression{7} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

load('MicroarrayDataWITHcustProbesUpdatedXXXCV.mat')
probes{8} = probeInformation;
expression{8} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

%select only genes that had more than one probe available
[~, indFilter] = intersect(probes{1}.EntrezID, duplicate_value);
% separatelly for RNAseq and others as the numbef of genes is different

% filter genes that had multiple probes
for k=1:length(probes)
    
    expression{k} = expression{k}(:,indFilter);
end

avCorr = zeros(8,8);

for i=1:8
    
    for j=1:8
        
        expr1 = expression{i};
        expr2 = expression{j};
        
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

nice_cmap = [flipud(make_cmap('orangered',50,30,0))];

figure; imagesc(avCorrPlot);set(gcf,'color','w');
colormap(nice_cmap)
caxis([0.5 1])
tickNames = {'Mean', 'Noise', 'PC', 'DS', 'Random1', 'Random2','Variance', 'CV'};
tickNamesORD = tickNames(ord);

xticks([1 2 3 4 5 6 7 8])
xticklabels(tickNamesORD);
yticks([1 2 3 4 5 6 7 8])
yticklabels(tickNamesORD);

cd ../../..
% make a plot with RNA seq
options.probeSelections = {'RNAseq'};
options.RNAseqThreshold = 0.2;
options.signalThreshold = 0.5; % filter QC
options.useCUSTprobes = true;

if Regeneratedata
S2_probes(options)
end

% load data
cd ('data/genes/processedData')
load('MicroarrayDataWITHcustProbesUpdatedXXXRNAseq.mat')
probes{9} = probeInformation;
expression{9} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});
% select genes in other versions and in this one  that overlap
entrezIDs = probes{1}.EntrezID(indFilter);

[genesMultiple, indFilter] = intersect(entrezIDs, duplicate_value); % separatelly for RNAseq and others as the numbef of genes is different
[genesMultipleRNAseq, indFilterRNA] = intersect(probes{9}.EntrezID, duplicate_value);

% filter genes that had multiple probes
for k=1:9
    if k==9
        expression{k} = expression{k}(:,indFilterRNA);
    else
        expression{k} = expression{k}(:,indFilter);
    end
end

[v1, indALL] = intersect(genesMultiple, genesMultipleRNAseq);
[v2, indRNA] = intersect(genesMultipleRNAseq,genesMultiple);
% calculate correlation between each way of choosing a probe

RNAcorr = zeros(8,1);
for i=9
    for j=1:8
        
        expr1 = expression{i}(:,indRNA);
        expr2 = expression{j}(:,indALL);
        
        correlation = zeros(size(expr1,2),1);
        
        for g=1:size(expr1,2)
            correlation(g) = corr(expr1(:,g), expr2(:,g), 'type', 'Spearman');
        end
        
        RNAcorr(j) = mean(correlation);
    end
end

[valS, indS] = sort(RNAcorr, 'descend');
namesS = tickNames(indS);

figure; 
set(gcf,'color','w');
hold on
for i = 1:length(valS)
    h=bar(i,valS(i));
    %errorbar(i,valS(i),stDEV(i),'.')
    if strcmp(namesS(i), 'CV')
        set(h,'FaceColor',[1 .87 .68],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
    elseif strcmp(namesS(i), 'Variance')
        set(h,'FaceColor',[1 0.8 0.6],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
    elseif strcmp(namesS(i), 'Random1')
        set(h,'FaceColor',[.96 .63 .55],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
    elseif strcmp(namesS(i), 'Random2')
        set(h,'FaceColor',[.96 .63 .55],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
    elseif strcmp(namesS(i), 'PC')
        set(h,'FaceColor',[.95 .6 .6],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
    elseif strcmp(namesS(i), 'Mean')
        set(h,'FaceColor',[.88 .56 .59],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
    elseif strcmp(namesS(i), 'Noise')
        set(h,'FaceColor',[.72 .43 .47],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
    elseif strcmp(namesS(i), 'DS')
        set(h,'FaceColor',[.75 .31 .38],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);                 
    end
end
hold off
ylim([0.5 1]); ylabel('Spearman correlation'); xlabel('Probe selection methods');
title(sprintf('Correlation to probes selected based on highest correlation to RNA-seq (%d)', length(indRNA)))
xticks([1 2 3 4 5 6 7 8])
xticklabels(namesS);
set(gca,'fontsize',15)

