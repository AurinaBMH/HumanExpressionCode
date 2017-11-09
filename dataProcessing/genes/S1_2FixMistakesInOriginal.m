clear all;
useCUSTprobes = true;
probeSelection = 'PC';% probeSelection = {'Mean', 'Variance', 'LessNoise', 'Random', 'PC'};
signalThreshold = -1; % percentage of samples that a selected probe has expression levels that are higher than background
%------------------------------------------------------------------------------
% Load the data
%------------------------------------------------------------------------------
cd ('data/genes/processedData');

if useCUSTprobes
    fprintf(1,'Loading the data with CUST probes and assigning variables\n')
    
    startFileName = 'MicroarrayDataWITHcust';
else
    fprintf(1,'Loading the data without CUST probes and assigning variables\n')
    startFileName = 'MicroarrayData';
end
fprintf(1,sprintf('Probe selection based on %s is chosen\n', probeSelection))
load(sprintf('%s.mat', startFileName));

cd ..
cd ('rawData');

[~, ind] = unique(DataTableProbe.GeneSymbol{1}, 'stable');
% find duplicate indices
duplicate_ind = setdiff(1:size(DataTableProbe.GeneSymbol{1},1), ind);

k=1; 
for gene = 1:length(duplicate_ind)
    
    % if gene names are not matching, exclude both of them
    symb = DataTableProbe.GeneSymbol{1}(duplicate_ind(gene));
    test_ind = find(strcmp(DataTableProbe.GeneSymbol{1}, symb{1}));
    % check if for all duplicated instances gene name is the same
    
    entrezID = DataTableProbe.EntrezID{1}(test_ind);
    doMatch = zeros(length(test_ind));
    for j=1:length(test_ind)
        for i=1:length(test_ind)
            
            doMatch(j,i) = entrezID(j)~=entrezID(i);
            
        end
    end
    % if they're not matching, record indexes to exclude later
    if sum(doMatch(:))~=0
        toCheckOriginal{k,:} = test_ind;
        k=k+1; 
    end

end

% remove duplicated rows in the cell
[r1, r2] = ndgrid(1:size(toCheckOriginal, 1));
duplicates = any(triu(arrayfun(@(r1, r2) isequal(toCheckOriginal(r1, :), toCheckOriginal(r2, :)), r1, r2), 1)); 
toCheckOriginal(duplicates, :) = []; 

% for each row in the cell check names and symbols of the data manually and
% confirm with the database
i=toCheckOriginal{20}
symb = DataTableProbe.GeneSymbol{1}(i)
name = DataTableProbe.GeneName{1}(i)
probe = DataTableProbe.ProbeName{1}(i)
gene = DataTableProbe.GeneID{1}(i)
entrez = DataTableProbe.EntrezID{1}(i)

% check each instance manually and rename or delete a probe
DataTableProbe.GeneSymbol{1}{2174} = 'GYPA';
DataTableProbe.GeneSymbol{1}{2175} = 'GYPA';

DataTableProbe.GeneSymbol{1}{34143} = 'INS-IGF2'; 

DataTableProbe.GeneSymbol{1}{4373} = 'NPHS1'; 
DataTableProbe.GeneSymbol{1}{4374} = 'NPHS1'; 

DataTableProbe.GeneSymbol{1}{26942} = 'PRSS3P2'; 

DataTableProbe.GeneSymbol{1}{39302} = 'SEBOX';
DataTableProbe.GeneSymbol{1}{42587} = 'SEBOX';
DataTableProbe.GeneSymbol{1}{47956} = 'SEBOX';

DataTableProbe.GeneSymbol{1}{7375} = 'ZNF221';
DataTableProbe.GeneSymbol{1}{7376} = 'ZNF221';

DataTableProbe.GeneSymbol{1}{7454} = 'ZNF230';
DataTableProbe.GeneSymbol{1}{7455} = 'ZNF230';

DataTableProbe.GeneSymbol{1}{7798} = 'SIGLEC5'; 
DataTableProbe.GeneSymbol{1}{7799} = 'SIGLEC5'; 
DataTableProbe.GeneSymbol{1}{47823} = 'SIGLEC5'; 

DataTableProbe.GeneSymbol{1}{7888} = 'UBE2M';
DataTableProbe.GeneSymbol{1}{7889} = 'UBE2M';
 
DataTableProbe.GeneSymbol{1}{17389} = 'AHRR'; 
DataTableProbe.GeneSymbol{1}{17390} = 'AHRR'; 

DataTableProbe.GeneSymbol{1}{24332} = []; % uncharacterized protein
DataTableProbe.GeneSymbol{1}{24333} = []; % uncharacterized protein

DataTableProbe.GeneSymbol{1}{36610} = 'UGT2A2'; 

DataTableProbe.GeneSymbol{1}{11621} = []; % maybe HYPK huntingtin interacting protein K, but not sure
DataTableProbe.GeneSymbol{1}{11622} = []; % maybe HYPK huntingtin interacting protein K, but not sure
DataTableProbe.GeneSymbol{1}{11623} = []; % maybe HYPK huntingtin interacting protein K, but not sure

DataTableProbe.GeneSymbol{1}{37878} = 'ZNF283';
DataTableProbe.GeneSymbol{1}{41133} = 'ZNF283';
DataTableProbe.GeneSymbol{1}{41930} = 'ZNF283';
DataTableProbe.GeneSymbol{1}{47952} = 'ZNF283';

DataTableProbe.GeneSymbol{1}{12185} = 'OR4F4'; 
DataTableProbe.GeneSymbol{1}{12186} = 'OR4F4'; 

DataTableProbe.GeneSymbol{1}{32243} = 'PCDHGC5';
DataTableProbe.GeneSymbol{1}{32244} = 'PCDHGC5';

DataTableProbe.GeneSymbol{1}{12606} = 'RABGEF1'; 
DataTableProbe.GeneSymbol{1}{12607} = 'RABGEF1';    
DataTableProbe.GeneSymbol{1}{40862} = 'RABGEF1'; 

DataTableProbe.GeneSymbol{1}{25502} = 'ANKHD1-EIF4EBP3';
DataTableProbe.GeneSymbol{1}{25503} = 'ANKHD1-EIF4EBP3';

DataTableProbe.GeneSymbol{1}{14542} = 'SARS2';
DataTableProbe.GeneSymbol{1}{14543} = 'SARS2'; 

DataTableProbe.GeneSymbol{1}{14947} = 'MARC2';
DataTableProbe.GeneSymbol{1}{14948} = 'MARC2';

DataTableProbe.GeneSymbol{1}{26205} = 'MARCH2'; 
DataTableProbe.GeneSymbol{1}{26206} = 'MARCH2'; 

% find all that are 
% '01-Mar-2002' and other dates and rename them
i=toCheckOriginal{20}
symb = DataTableProbe.GeneSymbol{1}(i)
name = DataTableProbe.GeneName{1}(i)
probe = DataTableProbe.ProbeName{1}(i)
gene = DataTableProbe.GeneID{1}(i)
entrez = DataTableProbe.EntrezID{1}(i)