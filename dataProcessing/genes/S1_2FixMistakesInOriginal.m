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

% check each instance manually and rename or delete a probe
% 1
DataTableProbe.GeneSymbol{1}{2174} = 'GYPA';
DataTableProbe.GeneSymbol{1}{2175} = 'GYPA';
% 2
DataTableProbe.GeneSymbol{1}{34143} = 'INS-IGF2';
% 3
DataTableProbe.GeneSymbol{1}{4373} = 'NPHS1';
DataTableProbe.GeneSymbol{1}{4374} = 'NPHS1';
% 4
DataTableProbe.GeneSymbol{1}{26942} = 'PRSS3P2';
% 5
DataTableProbe.GeneSymbol{1}{39302} = 'SEBOX';
DataTableProbe.GeneSymbol{1}{42587} = 'SEBOX';
DataTableProbe.GeneSymbol{1}{47956} = 'SEBOX';
% 6
DataTableProbe.GeneSymbol{1}{7375} = 'ZNF221';
DataTableProbe.GeneSymbol{1}{7376} = 'ZNF221';
% 7
DataTableProbe.GeneSymbol{1}{7454} = 'ZNF230';
DataTableProbe.GeneSymbol{1}{7455} = 'ZNF230';
% 8
DataTableProbe.GeneSymbol{1}{7798} = 'SIGLEC5';
DataTableProbe.GeneSymbol{1}{7799} = 'SIGLEC5';
DataTableProbe.GeneSymbol{1}{47823} = 'SIGLEC5';
% 9
DataTableProbe.GeneSymbol{1}{7888} = 'UBE2M';
DataTableProbe.GeneSymbol{1}{7889} = 'UBE2M';
% 10
DataTableProbe.GeneSymbol{1}{17389} = 'AHRR';
DataTableProbe.GeneSymbol{1}{17390} = 'AHRR';
% 11
DataTableProbe.GeneSymbol{1}{24332} = []; % uncharacterized protein
DataTableProbe.GeneSymbol{1}{24333} = []; % uncharacterized protein
% 12
DataTableProbe.GeneSymbol{1}{36610} = 'UGT2A2';
% 13
DataTableProbe.GeneSymbol{1}{11621} = []; % maybe HYPK huntingtin interacting protein K, but not sure
DataTableProbe.GeneSymbol{1}{11622} = []; % maybe HYPK huntingtin interacting protein K, but not sure
DataTableProbe.GeneSymbol{1}{11623} = []; % maybe HYPK huntingtin interacting protein K, but not sure
% 14
DataTableProbe.GeneSymbol{1}{37878} = 'ZNF283';
DataTableProbe.GeneSymbol{1}{41133} = 'ZNF283';
DataTableProbe.GeneSymbol{1}{41930} = 'ZNF283';
DataTableProbe.GeneSymbol{1}{47952} = 'ZNF283';
% 15
DataTableProbe.GeneSymbol{1}{12185} = 'OR4F4';
DataTableProbe.GeneSymbol{1}{12186} = 'OR4F4';

% 16
DataTableProbe.GeneSymbol{1}{32243} = 'PCDHGC5';
DataTableProbe.GeneSymbol{1}{32244} = 'PCDHGC5';
% 17
DataTableProbe.GeneSymbol{1}{12606} = 'RABGEF1';
DataTableProbe.GeneSymbol{1}{12607} = 'RABGEF1';
DataTableProbe.GeneSymbol{1}{40862} = 'RABGEF1';
% 18
DataTableProbe.GeneSymbol{1}{25502} = 'ANKHD1-EIF4EBP3';
DataTableProbe.GeneSymbol{1}{25503} = 'ANKHD1-EIF4EBP3';
% 19
DataTableProbe.GeneSymbol{1}{14542} = 'SARS2';
DataTableProbe.GeneSymbol{1}{14543} = 'SARS2';
% 20
DataTableProbe.GeneSymbol{1}{14947} = 'MARC2';
DataTableProbe.GeneSymbol{1}{14948} = 'MARC2';
% 21
DataTableProbe.GeneSymbol{1}{26205} = 'MARCH2';
DataTableProbe.GeneSymbol{1}{26206} = 'MARCH2';

% 21 is correct 'non-protein coding RNA 185' ir replaced with
% testis-specific transcript, Y-linked 14 (non-protein coding) in https://www.ncbi.nlm.nih.gov/gene/83869
% 22
DataTableProbe.GeneSymbol{1}{38137} = [];

% 23 is correct 'uncharacterized LOC100130093' ir replaced with
% synaptosomal-associated protein, 47kDa in https://www.ncbi.nlm.nih.gov/gene/?term=100130093
% 24
DataTableProbe.GeneSymbol{1}{21625} = 'MYL6B';
DataTableProbe.GeneSymbol{1}{21626} = 'MYL6B';
% 25
DataTableProbe.GeneSymbol{1}{22929} = 'ZNF584';
DataTableProbe.GeneSymbol{1}{22930} = 'ZNF584';
% 26
DataTableProbe.GeneSymbol{1}{24505} = 'DCDC1';
DataTableProbe.GeneSymbol{1}{24506} = 'DCDC1';
% 27
DataTableProbe.GeneSymbol{1}{24670} = 'TAS2R46';
DataTableProbe.GeneSymbol{1}{24671} = 'TAS2R46';
DataTableProbe.GeneSymbol{1}{37046} = 'TAS2R46';

% 28 very ambiguous annotations in the website, discontinued or changed to
% something else.
DataTableProbe.GeneSymbol{1}{25885} = [];
DataTableProbe.GeneSymbol{1}{25886} = [];
DataTableProbe.GeneSymbol{1}{25927} = [];
DataTableProbe.GeneSymbol{1}{25928} = [];
DataTableProbe.GeneSymbol{1}{29847} = [];
DataTableProbe.GeneSymbol{1}{29848} = [];
DataTableProbe.GeneSymbol{1}{35196} = [];
% 29
DataTableProbe.GeneSymbol{1}{26499} = 'PRR5-ARHGAP8';
DataTableProbe.GeneSymbol{1}{26500} = 'PRR5-ARHGAP8';
% 30
DataTableProbe.GeneSymbol{1}{27893} = 'OR9G1';
DataTableProbe.GeneSymbol{1}{27894} = 'OR9G1';
% 31
% MRC1L1 replaced with MRC1 in the database
DataTableProbe.GeneSymbol{1}{3480} = 'MRC1';
DataTableProbe.GeneSymbol{1}{28579} = 'MRC1';
DataTableProbe.GeneSymbol{1}{46284} = 'MRC1';
% 32
DataTableProbe.GeneSymbol{1}{30453} = 'GABARAPL1';
DataTableProbe.GeneSymbol{1}{30454} = 'GABARAPL1';
% 33
DataTableProbe.GeneSymbol{1}{10522} = 'NUDT4';
% 34
DataTableProbe.GeneSymbol{1}{30764} = 'TMEM189-UBE2V1';
DataTableProbe.GeneSymbol{1}{30765} = 'TMEM189-UBE2V1';
DataTableProbe.GeneSymbol{1}{30766} = 'TMEM189-UBE2V1';
DataTableProbe.GeneSymbol{1}{30767} = 'TMEM189-UBE2V1';
DataTableProbe.GeneSymbol{1}{30768} = 'TMEM189-UBE2V1';
DataTableProbe.GeneSymbol{1}{30769} = 'TMEM189-UBE2V1';
DataTableProbe.GeneSymbol{1}{30770} = 'TMEM189-UBE2V1';
% 35
DataTableProbe.GeneName{32506} = 'ELMO domain containing 1';
DataTableProbe.GeneName{32507} = 'ELMO domain containing 1';
DataTableProbe.GeneName{38909} = 'ELMO domain containing 1';
DataTableProbe.GeneName{40316} = 'ELMO domain containing 1';
DataTableProbe.GeneName{40447} = 'ELMO domain containing 1';
% 36
DataTableProbe.GeneSymbol{32505} = 'NBPF10';
DataTableProbe.GeneSymbol{35693} = 'NBPF11';
DataTableProbe.GeneSymbol{36886} = 'NBPF11';
% 37
DataTableProbe.GeneSymbol{38097} = 'CT47A3';
DataTableProbe.GeneSymbol{38120} = 'CT47A9';
% 38
DataTableProbe.GeneSymbol{36691} = 'NPIPA5';
DataTableProbe.GeneName{36691} = 'nuclear pore complex interacting protein family member A5';
% 39
DataTableProbe.GeneSymbol{34939} = 'HHLA1';
DataTableProbe.GeneSymbol{40045} = 'HHLA1';
% 40
DataTableProbe.GeneSymbol{42431} = 'JMJD7';
% 41
DataTableProbe.GeneSymbol{39835} = 'AGAP9';
DataTableProbe.GeneSymbol{43502} = 'AGAP9';
DataTableProbe.GeneSymbol{44951} = 'AGAP9';
% 42
DataTableProbe.GeneSymbol{44093} = 'APOC2';
DataTableProbe.GeneSymbol{44410} = 'APOC2';
% 43
DataTableProbe.GeneSymbol{44220} = 'ADH1A';
DataTableProbe.GeneSymbol{44221} = 'ADH1A';
% 44
DataTableProbe.GeneSymbol{45661} = 'OR8U1';
DataTableProbe.GeneSymbol{45662} = 'OR8U1';
DataTableProbe.GeneSymbol{45663} = 'OR8U1';
% 45
DataTableProbe.EntrezID(47836) = 401024;
DataTableProbe.GeneID(47836) = 125404;
% 46
DataTableProbe.GeneSymbol{48065} = 'ZNF93';
DataTableProbe.GeneSymbol{48066} = 'ZNF93';
DataTableProbe.GeneSymbol{48067} = 'ZNF93';
DataTableProbe.GeneSymbol{48068} = 'ZNF93';

%i=toCheckOriginal{46}
%symb = DataTableProbe.GeneSymbol{1}(i)
%name = DataTableProbe.GeneName{1}(i)
%probe = DataTableProbe.ProbeName{1}(i)
%gene = DataTableProbe.GeneID{1}(i)
%entrez = DataTableProbe.EntrezID{1}(i)

% find all that are
% '01-Mar-2002' and other dates and rename them
% to do that import ProbesReordered.xlsx file (first 69 items for entrez IDs - they correspond to gene symbols that are incorrectly read)
% take the unique values of entrez IDs and manually confirm their geneSymbols
% 27 unique values
uniqueEntrezID = unique(entrez_id);

DataTableProbe.GeneSymbol{170} = 'SEPT7';
DataTableProbe.GeneSymbol{171} = 'SEPT7';
DataTableProbe.GeneSymbol{172} = 'SEPT7';
DataTableProbe.GeneSymbol{35705} = 'SEPT7';

DataTableProbe.GeneSymbol{804} = 'SEPT1';
DataTableProbe.GeneSymbol{805} = 'SEPT1';

DataTableProbe.GeneSymbol{29849} = 'SEPT2';
DataTableProbe.GeneSymbol{29850} = 'SEPT2';
DataTableProbe.GeneSymbol{29851} = 'SEPT2';

DataTableProbe.GeneSymbol{4635} = 'SEPT5';
DataTableProbe.GeneSymbol{4636} = 'SEPT5';

DataTableProbe.GeneSymbol{4637} = 'SEPT4';
DataTableProbe.GeneSymbol{4638} = 'SEPT4';
DataTableProbe.GeneSymbol{4639} = 'SEPT4';
DataTableProbe.GeneSymbol{4640} = 'SEPT4';

DataTableProbe.GeneSymbol{9516} = 'MARCH6';
DataTableProbe.GeneSymbol{9517} = 'MARCH6';

DataTableProbe.GeneSymbol{39263} = 'SEPT9';
DataTableProbe.GeneSymbol{44929} = 'SEPT9';

DataTableProbe.GeneSymbol{11458} = 'SEPT6';
DataTableProbe.GeneSymbol{11459} = 'SEPT6';
DataTableProbe.GeneSymbol{11460} = 'SEPT6';

DataTableProbe.GeneSymbol{7383} = 'SEPT8';
DataTableProbe.GeneSymbol{38984} = 'SEPT8';
DataTableProbe.GeneSymbol{39630} = 'SEPT8';
DataTableProbe.GeneSymbol{42931} = 'SEPT8';

DataTableProbe.GeneSymbol{13311} = 'DEC1';
DataTableProbe.GeneSymbol{13312} = 'DEC1';

DataTableProbe.GeneSymbol{26205} = 'MARCH2';
DataTableProbe.GeneSymbol{26206} = 'MARCH2';

DataTableProbe.GeneSymbol{14111} = 'MARCH5';
DataTableProbe.GeneSymbol{14112} = 'MARCH5';

DataTableProbe.GeneSymbol{14947} = 'MARC2';
DataTableProbe.GeneSymbol{14948} = 'MARC2';

DataTableProbe.GeneSymbol{14977} = 'MARCH1';
DataTableProbe.GeneSymbol{14978} = 'MARCH1';
DataTableProbe.GeneSymbol{48100} = 'MARCH1';
DataTableProbe.GeneSymbol{48101} = 'MARCH1';
DataTableProbe.GeneSymbol{48102} = 'MARCH1';

DataTableProbe.GeneSymbol{15722} = 'SEPT11';
DataTableProbe.GeneSymbol{15723} = 'SEPT11';

DataTableProbe.GeneSymbol{31049} = 'SEPT3';
DataTableProbe.GeneSymbol{31050} = 'SEPT3';

DataTableProbe.GeneSymbol{25752} = 'MARCH4';
DataTableProbe.GeneSymbol{25753} = 'MARCH4';

DataTableProbe.GeneSymbol{19136} = 'MARCH7';
DataTableProbe.GeneSymbol{19137} = 'MARCH7';

DataTableProbe.GeneSymbol{20599} = 'MARCH9';
DataTableProbe.GeneSymbol{20600} = 'MARCH9';

DataTableProbe.GeneSymbol{20429} = 'MARCH3';
DataTableProbe.GeneSymbol{20430} = 'MARCH3';
DataTableProbe.GeneSymbol{20431} = 'MARCH3';

DataTableProbe.GeneSymbol{38720} = 'SEPT12';
DataTableProbe.GeneSymbol{42381} = 'SEPT12';

DataTableProbe.GeneSymbol{28284} = 'SEPT10';
DataTableProbe.GeneSymbol{28285} = 'SEPT10';
DataTableProbe.GeneSymbol{28286} = 'SEPT10';
DataTableProbe.GeneSymbol{28287} = 'SEPT10';

DataTableProbe.GeneSymbol{38051} = 'MARCH10';
DataTableProbe.GeneSymbol{42485} = 'MARCH10';

DataTableProbe.GeneSymbol{27259} = 'MARCH8';
DataTableProbe.GeneSymbol{27260} = 'MARCH8';

DataTableProbe.GeneSymbol{24601} = 'SEPT14';

DataTableProbe.GeneSymbol{37753} = 'MARCH11';
DataTableProbe.GeneSymbol{43033} = 'MARCH11';

DataTableProbe.GeneSymbol{34367} = 'SEPT7P2';
DataTableProbe.GeneSymbol{40501} = 'SEPT7P2';
DataTableProbe.GeneSymbol{41640} = 'SEPT7P2';
DataTableProbe.GeneSymbol{44995} = 'SEPT7P2';

%ix = find(DataTableProbe.EntrezID{1}==uniqueEntrezID(27))
%symb = DataTableProbe.GeneSymbol{1}(ix)
%name = DataTableProbe.GeneName{1}(ix)
%probe = DataTableProbe.ProbeName{1}(ix)
%gene = DataTableProbe.GeneID{1}(ix)
%entrez = DataTableProbe.EntrezID{1}(ix)
