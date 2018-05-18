% make a figure for probe correlation from scratch

%1. regenerate data for different probe selection methods

options.useCUSTprobes = true;
afterQC = true; 
Regeneratedata = false;

if afterQC
    options.signalThreshold = 0.5; % no QC filtering
    fileName = '';
else
    options.signalThreshold = -1; % no QC filtering
    fileName = 'noQC'; 
end

% files will be saved with standart name noQC

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
    
    options.probeSelections = {'Random1'};
    S2_probes(options)
    
    options.probeSelections = {'Random2'};
    S2_probes(options)
    
    options.probeSelections = {'Variance'};
    S2_probes(options)
    
    options.probeSelections = {'CV'};
    S2_probes(options)
    
    options.probeSelections = {'maxCorrelation_intensity'};
    S2_probes(options)
    
    options.probeSelections = {'maxCorrelation_variance'};
    S2_probes(options)
end
% load the initial data
% created with there options
% options.ExcludeCBandBS =  true;
% options.useCUSTprobes = true;
% options.updateProbes = 'reannotator';

cd ('data/genes/processedData')
load('MicroarrayDataWITHcustProbesUpdatedXXX.mat')
if afterQC

options.signalThreshold = 0.5; %  QC filtering
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

else
    
[v, ind] = unique(DataTableProbe.EntrezID{1});
duplicate_ind = setdiff(1:size(DataTableProbe.EntrezID{1}, 1), ind);
duplicate_value = unique(DataTableProbe.EntrezID{1}(duplicate_ind));
percentage = (length(duplicate_value)/length(unique(DataTableProbe.EntrezID{1})))*100;

end

format compact
percentage

% Load probes, selected using different methods that were just generated
load(sprintf('MicroarrayDataWITHcustProbesUpdatedXXXMean%s.mat', fileName))
probes{1} = probeInformation;
expression{1} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

load(sprintf('MicroarrayDataWITHcustProbesUpdatedXXXLessNoise%s.mat', fileName))
probes{2} = probeInformation;
expression{2} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

load(sprintf('MicroarrayDataWITHcustProbesUpdatedXXXPC%s.mat', fileName))
probes{3} = probeInformation;
expression{3} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

load(sprintf('MicroarrayDataWITHcustProbesUpdatedXXXDS%s.mat', fileName))
probes{4} = probeInformation;
expression{4} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

load(sprintf('MicroarrayDataWITHcustProbesUpdatedXXXRandom1%s.mat', fileName))
probes{5} = probeInformation;
expression{5} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

load(sprintf('MicroarrayDataWITHcustProbesUpdatedXXXRandom2%s.mat', fileName))
probes{6} = probeInformation;
expression{6} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

load(sprintf('MicroarrayDataWITHcustProbesUpdatedXXXVariance%s.mat', fileName))
probes{7} = probeInformation;
expression{7} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

load(sprintf('MicroarrayDataWITHcustProbesUpdatedXXXCV%s.mat', fileName))
probes{8} = probeInformation;
expression{8} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

load(sprintf('MicroarrayDataWITHcustProbesUpdatedXXXmaxIntensity%s.mat', fileName))
probes{9} = probeInformation;
expression{9} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

load(sprintf('MicroarrayDataWITHcustProbesUpdatedXXXmaxCorrelation_intensity%s.mat', fileName))
probes{10} = probeInformation;
expression{10} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});

load(sprintf('MicroarrayDataWITHcustProbesUpdatedXXXmaxCorrelation_variance%s.mat', fileName))
probes{11} = probeInformation;
expression{11} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});
%select only genes that had more than one probe available
[genesMultiple, indFilter] = intersect(probes{1}.EntrezID, duplicate_value); % separatelly for RNAseq and others as the numbef of genes is different

[~, indFilter] = intersect(probes{1}.EntrezID, duplicate_value);
% filter genes that had multiple probes
for k=1:length(probes)
    
    expression{k} = expression{k}(:, indFilter);
    
end

avCorr = zeros(11,11);
stdCorr = zeros(11,11);

for i=1:11
    
    
    for j=i+1:11
        
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

h = figure; imagesc(avCorrPlot);
set(gcf,'color','w');
colormap(nice_cmap)
caxis([0.5 1])
%alpha(0.5)
tickNames = {'mean', 'signal proportion', 'PC', 'DS', 'Rand_1', 'Rand_2','Variance', 'CV', 'max intensity', 'correlation intensity', 'correlation variance'};
tickNamesORD = tickNames(ord);

xticks([1 2 3 4 5 6 7 8 9 10 11])
xticklabels(tickNamesORD);
xtickangle(45)
yticks([1 2 3 4 5 6 7 8 9 10 11])
yticklabels(tickNamesORD);
ytickangle(45)
set(gca,'FontSize', 14)


% make a plot with RNA seq
options.probeSelections = {'RNAseq'};
options.RNAseqThreshold = 0.2;
options.useCUSTprobes = true;
if Regeneratedata
    S2_probes(options)
end

% load data
load(sprintf('MicroarrayDataWITHcustProbesUpdatedXXXRNAseq%s.mat', fileName))
probes{12} = probeInformation;
expression{12} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, expressionAll{4}, expressionAll{5},expressionAll{6});
% select genes in other versions and in this one  that overlap
entrezIDs = probes{1}.EntrezID(indFilter);

[genesMultiple, indFilter] = intersect(entrezIDs, duplicate_value); % separatelly for RNAseq and others as the numbef of genes is different
[genesMultipleRNAseq, indFilterRNA] = intersect(probes{12}.EntrezID, duplicate_value);

% filter genes that had multiple probes
for k=1:12
    if k==12
        expression{k} = expression{k}(:,indFilterRNA);
    else
        expression{k} = expression{k}; %; (:,indFilter);
    end
end

[v1, indALL] = intersect(genesMultiple, genesMultipleRNAseq);
[v2, indRNA] = intersect(genesMultipleRNAseq,genesMultiple);
% calculate correlation between each way of choosing a probe

RNAcorr = zeros(11,1);
stDEV = zeros(11,1);
% calculate for RNAseq
i=12;
for j=1:11
    
    expr1 = expression{i}(:,indRNA);
    expr2 = expression{j}(:,indALL);
    
    
    correlation = zeros(size(expr1,2),1);
    
    for g=1:size(expr1,2)
        correlation(g) = corr(expr1(:,g), expr2(:,g), 'type', 'Spearman');
    end
    
    
    RNAcorr(j) = mean(correlation);
    stDEV(j) = std(correlation);
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
    elseif strcmp(namesS(i), 'correlation intensity')
        set(h,'FaceColor',[.72 .43 .47],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
    elseif strcmp(namesS(i), 'correlation variance')
        
        set(h,'FaceColor',[1 0.8 0.6],'EdgeColor',[0.45 0.45 0.45],'LineWidth',1.5);
    end
end
hold off
ylim([0.5 1]); ylabel('Spearman correlation'); %xlabel('Probe selection methods');
%title(sprintf(' (%d)', length(indRNA)))
xticks([1 2 3 4 5 6 7 8 9 10 11])
%legend(sprintf('%d genes', length(indRNA)))
xticklabels(namesS);
xtickangle(45)
set(gca,'FontSize', 14)