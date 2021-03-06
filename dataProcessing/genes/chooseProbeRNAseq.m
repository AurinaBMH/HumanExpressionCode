% function to select probe based on RNA seq values
function indProbe = chooseProbeRNAseq(DataTable, DataTableProbe, threshold)

numGenes = length(unique(DataTableProbe.EntrezID{1}));
correlations = cell(numGenes,2);
for subject=1:2
    
    folderName = sprintf('rnaseq_donor0%d', subject);
    cd (folderName)
    
    RNAseqinfo = importRNAseqInfo('Genes.csv');
    RNAseqTPM = importRNAseq('RNAseqTPM.csv');
    
    % replace gene names with IDs in the expression file
    RNAseqTPM = cell2mat(RNAseqTPM(:,2:end));
    RNAseqGene = RNAseqinfo.NCBIgeneID;
    
    % import info about RNAseq tissue samples
    structureIDrna = importRNAseqStructInf('SampleAnnot.csv');
    % import information about microarray samples
    structureIDmic = DataTable.SampleID{subject};
    
    % import expression values for microarray
    microarray = DataTable.Expression{subject};
    microarrayGene = DataTableProbe.EntrezID{1};
    uniqueGenes = unique(microarrayGene);
    
    % find overlapping structures between microarray and RNAseq
    overlapStructures = intersect(structureIDmic, structureIDrna);
    % average if there are multiple samples in the same structure
    expmic = zeros(size(microarray,1), length(overlapStructures));
    exprna = zeros(size(RNAseqTPM,1), length(overlapStructures));
    for struct=1:length(overlapStructures)
        % find ind for microarray and average expression for those elements
        indmic = find(structureIDmic==overlapStructures(struct));
        expmic(:,struct) = mean(microarray(:,indmic),2);
        % find ind for RNAseq
        indrna = find(structureIDrna==overlapStructures(struct));
        exprna(:,struct) = mean(RNAseqTPM(:,indrna),2);
    end
    
    % correlate each microarray probe with RNAseq of the corresponding gene, if
    % gene is not available put NaNs.
    
    for g=1:numGenes
        indRNA = find(RNAseqGene==uniqueGenes(g));
        indMIC = find(microarrayGene==uniqueGenes(g));
        % for each probe found in microarray correlate it with RNAseq
        if isempty(indRNA)
            correlations{g,subject} = NaN;
        else
            C = zeros(length(indMIC),1);
            for p=1:length(indMIC)
                C(p) = corr(exprna(indRNA,:)',expmic(indMIC(p),:)', 'type', 'Spearman');
            end

            correlations{g, subject} = C;

        end
    end
    cd ..
end

% take an average of results between two subjects:
% if multiple values are available, choose the one which is on average
% higher
% keep the probe if on average between two subjects correlation is higher
% than 0.1 (ir some other threshold)
avgCorr = cell(numGenes,1);
indProbe = zeros(numGenes,1);
for gene = 1:numGenes
    cors = [correlations{gene,1}, correlations{gene,2}];
    avc = mean(cors,2);
    
    avgCorr{gene} = avc;
    [chosenVal, chosenInd] = max(avc);
    if chosenVal > threshold
        indProbe(gene) = chosenInd;
    else
        indProbe(gene) = NaN;
    end
    
end
end

