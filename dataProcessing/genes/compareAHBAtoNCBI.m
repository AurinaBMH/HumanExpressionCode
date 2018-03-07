
% script to check how many gene IDs/symbols/names do not match in the original AHBA dataset compared with NSCI data (downloaded 05/03/2018)
clear all; 
cd ('data/genes/rawData');

fprintf(1,'Loading Probes.xlsx file\n')
FileProbes = 'Probes.xlsx';
ProbeTable = readtable(FileProbes);
ProbeID = ProbeTable.probe_id;
EntrezID = ProbeTable.entrez_id;
ProbeName =  ProbeTable.probe_name;
GeneID = ProbeTable.gene_id;
GeneSymbol = ProbeTable.gene_symbol;
GeneName = ProbeTable.gene_name;

fprintf(1,'Loading ncbi data\n')
Homosapiens = importNCBIgenefile('Homo_sapiens_gene_info.txt');
%load('HomoSampiens_geneInfo20171111.mat'); 

[matches, GeneSymbol, GeneName, nrUpdated, checkEntrezID, MissingProbes] = checkGene(Homosapiens, GeneSymbol, GeneName, EntrezID); 


% [uGene, ind] = unique(GeneSymbol); 
% 
% gOK = zeros(length(uGene),1); 
% 
% for i=1:length(uGene)
%     gene = uGene(i); 
%     geneInd = find(strcmp(gene, GeneSymbol));
%     if length(geneInd) >1
%         gID = unique(EntrezID(geneInd));
%         gName = unique(GeneName(geneInd));
%         
%         if length(gID) == length(gName)
%             gOK(i) = 1;
%         else
%             gOK(i) = 2;
%  
%         end
%         
%         
%     end
% end
% 
