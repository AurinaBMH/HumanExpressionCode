%% Author: Aurina
%% Date modified: 2016-03-22
%% Date modified: 2017-07-14
%% Date modified: 2017-08-29 - noise level data added
%% Date modified: 2017-11-10 - manual gene symbol naming fixed
%% This script:
%   1. Loads all microarray data from excell files for each subject
%   2. Excludes custom probes;
%   3. Excludes probes with missing entrezIDs
%   4. Saves expression data, coordinates, sample structure names for all samples
%   5. Saves data for separate subjects to DataTable
%   6. Saves data for all subjects combined as variables 'MicorarrayData.mat' file
%%
%------------------------------------------------------------------------------
% Choose options
%------------------------------------------------------------------------------
% RemoveCUSTProbes = true; will exclude CUST probes
% ExcludeCBandBS = true; will exclude samples from brainstem and cerebellum
useCUSTprobes = true;
ExcludeCBandBS = true;

if useCUSTprobes
    fprintf(1,'CUST probes will be included\n')
    startFileName = 'MicroarrayDataWITHcust';
else
    fprintf(1,'CUST probes will be excluded\n')
    startFileName = 'MicroarrayData';
end

cd ('data/genes/rawData');
%% load probe information (same for all subjects)
fprintf(1,'Loading Probes.xlsx file\n')
FileProbes = 'Probes.xlsx';
% in Probes.xlsx file probes from 1070378 down are deleted as there are no information about gene samples.
ProbeTable = readtable(FileProbes);
ProbeID = ProbeTable.probe_id;
EntrezID = ProbeTable.entrez_id;
ProbeName =  ProbeTable.probe_name;
GeneID = ProbeTable.gene_id;
GeneSymbol = ProbeTable.gene_symbol;
GeneName = ProbeTable.gene_name;


% assign NaN values for all probes with missing entrezIDs
% this is the final list of probes to be used in max var calculations
ProbeID(isnan(EntrezID)) = NaN;
fprintf(1,'%d probes with missing entrez IDs\n', sum(isnan(EntrezID)));
% creat a Data cell to store the output
headerdata = {'Expression' , 'MMcoordinates', 'StructureName', 'MRIvoxCoordinates', 'Noise'};
headerprobe = { 'ProbeID', 'EntrezID','ProbeName', 'GeneID', 'GeneSymbol', 'GeneName'};
Data = cell(6,5);
DataProbe = cell(1,6);
%%
%------------------------------------------------------------------------------
%Go to each subject's directory and take the data
%------------------------------------------------------------------------------

for subj=1:6
    fprintf(1,'Loading data for %u subject\n', subj)
    folder = sprintf('normalized_microarray_donor0%d', subj);
    cd (folder);
    %%load information specific for each subject
    FileMicroarray = 'MicroarrayExpression.csv';
    FileAnnot = 'SampleAnnot.xlsx';
    FileNoise = 'PACall.csv';
    Expression = csvread(FileMicroarray);
    noise = csvread(FileNoise);
    % rows from 57860 are removed corresponding to probes from 1070378 down as there are no information about gene samples.
    % Expression = removerows(Expression,57860:size(Expression,1));
    Expression(:,1) = [];                         % exclude probe IDs from expression matrix
    [~,~,SlabType] = xlsread(FileAnnot, 'D:D');
    [~,~, StructureName] = xlsread(FileAnnot, 'F:F');
    SlabType(1) = [];                           % remove headline
    StructureName(1) = [];                      % remove headline
    MMcoordinates = xlsread(FileAnnot, 'K:M');
    MRIvoxCoordinates = xlsread(FileAnnot, 'H:J');
    [~,probeList] = intersect(noise(:,1),ProbeID, 'stable');
    noise = noise(probeList,2:end);
    % keep only non custom probes with existing entrezIDs
    Expression(isnan(ProbeID),:) = [];

    % To exclude expression data and coordinates for braintem (BS) and cerebellum (CB)
    % exclude columns in expression and rows in coordinates if slabtype is CB or BS
    if ExcludeCBandBS
        fprintf('Excluding brainstem and cerebellum data\n')

        BS = strfind(SlabType, 'BS'); BSind = find(~cellfun(@isempty,BS));
        CB = strfind(SlabType, 'CB'); CBind = find(~cellfun(@isempty,CB));
        BSandCBind = [BSind;CBind];

        fprintf(1,'%d cerebellum and brainstem samples to remove\n', length(BSandCBind))

        Expression(:,BSandCBind) = NaN;
        MMcoordinates(BSandCBind,:) = NaN;
        MRIvoxCoordinates(BSandCBind,:) = NaN;
        noise(:,BSandCBind) = NaN;
        StructureName(BSandCBind) = {NaN};
    end


    % for nan columns
    % keep only existing expression values
    Expression = Expression(:,all(~isnan(Expression)));
    % keep only existing coordinates
    MMcoordinates = MMcoordinates(all(~isnan(MMcoordinates),2),:); % for nan rows
    MRIvoxCoordinates = MRIvoxCoordinates(all(~isnan(MRIvoxCoordinates),2),:); % for nan rows
    noise = noise(:,all(~isnan(noise)));
    % keep only existing structure names
    StructureName(cellfun(@(StructureName) any(isnan(StructureName)),StructureName)) = [];

    % assign output to Data cell;
    Data{subj,1} = Expression;
    Data{subj,2} = MMcoordinates;
    Data{subj,3} = StructureName;
    Data{subj,4} = MRIvoxCoordinates;
    Data{subj,5} = noise;
    cd ..

end
%%
%------------------------------------------------------------------------------
% Keep only existing ProbeNames, EntrezIDs and ProbeIDs and other gene related information.
%------------------------------------------------------------------------------
fprintf(1,'Removing irrelevant probes\n')
fprintf(1,'Removing %d irrelevant probes\n', sum(isnan(ProbeID)))

ProbeName(isnan(ProbeID)) = [];
EntrezID(isnan(ProbeID)) = [];
GeneID(isnan(ProbeID)) = [];
GeneSymbol(isnan(ProbeID)) = [];
GeneName(isnan(ProbeID)) = [];
ProbeID(isnan(ProbeID)) = [];

%------------------------------------------------------------------------------
% Check for mistakes in probe naming and fix them plus update expression values for those probes that are excluded
%------------------------------------------------------------------------------
% 
% % find all that are
% % '01-Mar-2002' and other dates and rename them
% % to do that import ProbesReordered.xlsx file (first 69 items for entrez IDs - they correspond to gene symbols that are incorrectly read)
% % take the unique values of entrez IDs and manually confirm their geneSymbols
% % 27 unique values (variable called entrez_id)
% % uniqueEntrezID = unique(entrez_id);
% 
% GeneSymbol{170} = 'SEPT7';
% GeneSymbol{171} = 'SEPT7';
% GeneSymbol{172} = 'SEPT7';
% GeneSymbol{35705} = 'SEPT7';
% 
% GeneSymbol{804} = 'SEPT1';
% GeneSymbol{805} = 'SEPT1';
% 
% GeneSymbol{29849} = 'SEPT2';
% GeneSymbol{29850} = 'SEPT2';
% GeneSymbol{29851} = 'SEPT2';
% 
% GeneSymbol{4635} = 'SEPT5';
% GeneSymbol{4636} = 'SEPT5';
% 
% GeneSymbol{4637} = 'SEPT4';
% GeneSymbol{4638} = 'SEPT4';
% GeneSymbol{4639} = 'SEPT4';
% GeneSymbol{4640} = 'SEPT4';
% 
% GeneSymbol{9516} = 'MARCH6';
% GeneSymbol{9517} = 'MARCH6';
% 
% GeneSymbol{39263} = 'SEPT9';
% GeneSymbol{44929} = 'SEPT9';
% 
% GeneSymbol{11458} = 'SEPT6';
% GeneSymbol{11459} = 'SEPT6';
% GeneSymbol{11460} = 'SEPT6';
% 
% GeneSymbol{37383} = 'SEPT8';
% GeneSymbol{38984} = 'SEPT8';
% GeneSymbol{39630} = 'SEPT8';
% GeneSymbol{42931} = 'SEPT8';
% 
% GeneSymbol{13311} = 'DEC1';
% GeneSymbol{13312} = 'DEC1';
% 
% GeneSymbol{26205} = 'MARCH2';
% GeneSymbol{26206} = 'MARCH2';
% 
% GeneSymbol{14111} = 'MARCH5';
% GeneSymbol{14112} = 'MARCH5';
% 
% GeneSymbol{14947} = 'MARC2';
% GeneSymbol{14948} = 'MARC2';
% 
% GeneSymbol{14977} = 'MARCH1';
% GeneSymbol{14978} = 'MARCH1';
% GeneSymbol{48100} = 'MARCH1';
% GeneSymbol{48101} = 'MARCH1';
% GeneSymbol{48102} = 'MARCH1';
% 
% GeneSymbol{15722} = 'SEPT11';
% GeneSymbol{15723} = 'SEPT11';
% 
% GeneSymbol{31049} = 'SEPT3';
% GeneSymbol{31050} = 'SEPT3';
% 
% GeneSymbol{25752} = 'MARCH4';
% GeneSymbol{25753} = 'MARCH4';
% 
% GeneSymbol{19136} = 'MARCH7';
% GeneSymbol{19137} = 'MARCH7';
% 
% GeneSymbol{20599} = 'MARCH9';
% GeneSymbol{20600} = 'MARCH9';
% 
% GeneSymbol{20429} = 'MARCH3';
% GeneSymbol{20430} = 'MARCH3';
% GeneSymbol{20431} = 'MARCH3';
% 
% GeneSymbol{38720} = 'SEPT12';
% GeneSymbol{42381} = 'SEPT12';
% 
% GeneSymbol{28284} = 'SEPT10';
% GeneSymbol{28285} = 'SEPT10';
% GeneSymbol{28286} = 'SEPT10';
% GeneSymbol{28287} = 'SEPT10';
% 
% GeneSymbol{38051} = 'MARCH10';
% GeneSymbol{42485} = 'MARCH10';
% 
% GeneSymbol{27259} = 'MARCH8';
% GeneSymbol{27260} = 'MARCH8';
% 
% GeneSymbol{24601} = 'SEPT14';
% 
% GeneSymbol{37753} = 'MARCH11';
% GeneSymbol{43033} = 'MARCH11';
% 
% GeneSymbol{34367} = 'SEPT7P2';
% GeneSymbol{40501} = 'SEPT7P2';
% GeneSymbol{41640} = 'SEPT7P2';
% GeneSymbol{44995} = 'SEPT7P2';
% 
% % example of code how each entry was checked
% %ix = find(EntrezID{1}==uniqueEntrezID(27))
% %symb = GeneSymbol{1}(ix)
% %name = GeneName{1}(ix)
% %probe = ProbeName{1}(ix)
% %gene = GeneID{1}(ix)
% %entrez = EntrezID{1}(ix)
% 
% %------------------------------------------------------------------------------
% % Check for mistakes in probe naming 
% %------------------------------------------------------------------------------
% %changetoNumber = find(cell2mat(cellfun(@(x)any(isempty(x)),GeneSymbol,'UniformOutput',false))); 
% 
% toCheckSymbols = selectDuplicates(GeneSymbol, EntrezID);
% 
% % for each row in the cell check names and symbols of the data manually and
% % confirm with the database
% 
% % check each instance manually and rename or delete a probe
% % 1
% GeneSymbol{2174} = 'GYPA';
% GeneSymbol{2175} = 'GYPA';
% % 2
% GeneSymbol{34143} = 'INS-IGF2';
% % 3
% GeneSymbol{4373} = 'NPHS1';
% GeneSymbol{4374} = 'NPHS1';
% % 4
% GeneSymbol{26942} = 'PRSS3P2';
% % 5
% GeneSymbol{39302} = 'SEBOX';
% GeneSymbol{42587} = 'SEBOX';
% GeneSymbol{47956} = 'SEBOX';
% % 6
% GeneSymbol{7375} = 'ZNF221';
% GeneSymbol{7376} = 'ZNF221';
% % 7%%%%%%%
% GeneSymbol{7454} = 'ZNF230';
% GeneSymbol{7455} = 'ZNF230';
% % 8
% GeneSymbol{7798} = 'SIGLEC5';
% GeneSymbol{7799} = 'SIGLEC5';
% GeneSymbol{47823} = 'SIGLEC5';
% % 9
% GeneSymbol{7888} = 'UBE2M';
% GeneSymbol{7889} = 'UBE2M';
% % 10
% GeneSymbol{17389} = 'AHRR';
% GeneSymbol{17390} = 'AHRR';
% % 11
% GeneSymbol{24332} = NaN; % uncharacterized protein
% GeneSymbol{24333} = NaN;% uncharacterized protein
% % 12
% GeneSymbol{36610} = 'UGT2A2';
% % 13
% GeneSymbol{11621} = NaN;% maybe HYPK huntingtin interacting protein K, but not sure
% GeneSymbol{11622} = NaN;% maybe HYPK huntingtin interacting protein K, but not sure
% GeneSymbol{11623} = NaN;% maybe HYPK huntingtin interacting protein K, but not sure
% % 14
% GeneSymbol{37878} = 'ZNF283';
% GeneSymbol{41133} = 'ZNF283';
% GeneSymbol{41930} = 'ZNF283';
% GeneSymbol{47952} = 'ZNF283';
% % 15
% GeneSymbol{12185} = 'OR4F4';
% GeneSymbol{12186} = 'OR4F4';
% 
% % 16
% GeneSymbol{32243} = 'PCDHGC5';
% GeneSymbol{32244} = 'PCDHGC5';
% % 17
% GeneSymbol{12606} = 'RABGEF1';
% GeneSymbol{12607} = 'RABGEF1';
% GeneSymbol{40862} = 'RABGEF1';
% % 18
% GeneSymbol{25502} = 'ANKHD1-EIF4EBP3';
% GeneSymbol{25503} = 'ANKHD1-EIF4EBP3';
% % 19
% GeneSymbol{14542} = 'SARS2';
% GeneSymbol{14543} = 'SARS2';
% 
% % 20 %bothed edited before
% % GeneSymbol{14947} = 'MARC2';
% % GeneSymbol{14948} = 'MARC2';
% % % 21
% % GeneSymbol{26205} = 'MARCH2';
% % GeneSymbol{26206} = 'MARCH2';
% 
% % 21 is correct 'non-protein coding RNA 185' ir replaced with
% % testis-specific transcript, Y-linked 14 (non-protein coding) in https://www.ncbi.nlm.nih.gov/gene/83869
% % change entrez id to the updated one
% EntrezID(15057) = 83869; 
% EntrezID(15058) = 83869; 
% % 22
% GeneSymbol{38137} = NaN;
% 
% % 23 is correct 'uncharacterized LOC100130093' ir replaced with
% % synaptosomal-associated protein, 47kDa in https://www.ncbi.nlm.nih.gov/gene/?term=100130093
% % change entrez id to the updated one
% EntrezID(35442) = 116841; 
% 
% % 24
% GeneSymbol{21625} = 'MYL6B';
% GeneSymbol{21626} = 'MYL6B';
% % 25
% GeneSymbol{22929} = 'ZNF584';
% GeneSymbol{22930} = 'ZNF584';
% % 26
% GeneSymbol{24505} = 'DCDC1';
% GeneSymbol{24506} = 'DCDC1';
% % 27
% GeneSymbol{24670} = 'TAS2R46';
% GeneSymbol{24671} = 'TAS2R46';
% GeneSymbol{37046} = 'TAS2R46';
% 
% % 28 very ambiguous annotations in the website, discontinued or changed to
% % something else.
% GeneSymbol{25885} = NaN; 
% GeneSymbol{25886} = NaN; 
% GeneSymbol{25927} = NaN; 
% GeneSymbol{25928} = NaN; 
% GeneSymbol{29847} = NaN; 
% GeneSymbol{29848} = NaN; 
% GeneSymbol{35196} = NaN; 
% % 29
% GeneSymbol{26499} = 'PRR5-ARHGAP8';
% GeneSymbol{26500} = 'PRR5-ARHGAP8';
% % 30
% GeneSymbol{27893} = 'OR9G1';
% GeneSymbol{27894} = 'OR9G1';
% % 31
% % MRC1L1 replaced with MRC1 in the database
% GeneSymbol{3480} = 'MRC1';
% GeneSymbol{28579} = 'MRC1';
% EntrezID(28579) = 4360;
% GeneSymbol{46284} = 'MRC1';
% % 32
% GeneSymbol{30453} = 'GABARAPL1';
% GeneSymbol{30454} = 'GABARAPL1';
% % 33
% GeneSymbol{10522} = 'NUDT4';
% % 34
% GeneSymbol{30764} = 'TMEM189-UBE2V1';
% GeneSymbol{30765} = 'TMEM189-UBE2V1';
% GeneSymbol{30766} = 'TMEM189-UBE2V1';
% GeneSymbol{30767} = 'TMEM189-UBE2V1';
% GeneSymbol{30768} = 'TMEM189-UBE2V1';
% GeneSymbol{30769} = 'TMEM189-UBE2V1';
% GeneSymbol{30770} = 'TMEM189-UBE2V1';
% % 35
% GeneName{32506} = 'ELMO domain containing 1';
% GeneName{32507} = 'ELMO domain containing 1';
% GeneName{38909} = 'ELMO domain containing 1';
% GeneName{40316} = 'ELMO domain containing 1';
% GeneName{40447} = 'ELMO domain containing 1';
% EntrezID(32506) = 55531; %643923; 
% EntrezID(32507) = 55531; %643923; 
% GeneID(32506) = 34820; 
% GeneID(32507) = 34820; 
% % 36
% GeneSymbol{32505} = 'NBPF10';
% GeneSymbol{35693} = 'NBPF11';
% GeneSymbol{36886} = 'NBPF11';
% % 37
% GeneSymbol{38097} = 'CT47A3';
% GeneSymbol{38120} = 'CT47A9';
% % 38
% GeneSymbol{36691} = 'NPIPA5';
% GeneName{36691} = 'nuclear pore complex interacting protein family member A5';
% % 39
% GeneSymbol{34939} = 'HHLA1';
% GeneSymbol{40045} = 'HHLA1';
% % 40
% GeneSymbol{42431} = 'JMJD7';
% % 41
% GeneSymbol{39835} = 'AGAP9';
% GeneSymbol{43502} = 'AGAP9';
% GeneSymbol{44951} = 'AGAP9';
% % 42
% GeneSymbol{44093} = 'APOC2';
% GeneSymbol{44410} = 'APOC2';
% % 43
% GeneSymbol{44220} = 'ADH1A';
% GeneSymbol{44221} = 'ADH1A';
% % 44
% GeneSymbol{45661} = 'OR8U1';
% GeneSymbol{45662} = 'OR8U1';
% GeneSymbol{45663} = 'OR8U1';
% % 45
% EntrezID(47836) = 401024;
% GeneID(47836) = 125404;
% % 46
% GeneSymbol{48065} = 'ZNF93';
% GeneSymbol{48066} = 'ZNF93';
% GeneSymbol{48067} = 'ZNF93';
% GeneSymbol{48068} = 'ZNF93';
% % % example of code how each entry was checked
% % i=toCheckSymbols{6}
% % symb = GeneSymbol(i)
% % name = GeneName(i)
% % probe = ProbeName(i)
% % gene = GeneID(i)
% % entrez = EntrezID(i)
% %------------------------------------------------------------------------------
% % when these changes are done, do the second round of rhecks for
% % duplicates, because geneSymbols are changed, there might be more
% % duplicated items now 
% %------------------------------------------------------------------------------
% 
% changetoNumber = find(cell2mat(cellfun(@(x)any(isnan(x)),GeneSymbol,'UniformOutput',false))); 
% toCheckSymbolsAdditional = selectDuplicates(GeneSymbol, EntrezID, changetoNumber);
% 
% GeneSymbol{42970} = 'ZNF222';
% GeneSymbol{42996} = 'ZNF222';
%  
% GeneSymbol{15901} = 'PCDHGA2';
% GeneSymbol{15902} = 'PCDHGA2';
% 
% %------------------------------------------------------------------------------
% % Check for mistakes in GeneIDs
% %------------------------------------------------------------------------------
% 
% toCheckGeneIDs = selectDuplicates(EntrezID, GeneID); 
% 
% % add corrections
% 
% GeneID(15057) = 58103;
% GeneID(15058) = 58103;
% 
% GeneID(35442) = 78153; 
% 
% GeneID(28579) = 4334; 
% 
% %------------------------------------------------------------------------------
% % for those where gene symbol is [], remove them
% %------------------------------------------------------------------------------
% 
% emptyCells = find(cell2mat(cellfun(@(x)any(isnan(x)),GeneSymbol,'UniformOutput',false))); 
% fprintf(1,'Removing %d mislabeled probes \n', sum(emptyCells))
% fprintf(1,'Renaming misslabeled genes \n')
% ProbeID(emptyCells) = [];
% EntrezID(emptyCells) = [];
% ProbeName(emptyCells) = [];
% GeneID(emptyCells) = [];
% GeneSymbol(emptyCells) = [];
% GeneName(emptyCells) = [];
% 
% 
% % do the checks again to make sure all fixed
% toCheckSymbols2 = selectDuplicates(GeneSymbol, EntrezID);
% toCheckGeneIDs2 = selectDuplicates(EntrezID, GeneID); 

%------------------------------------------------------------------------------
% Check for mistakes in probe naming and fix them plus update expression values for those probes that are excluded
%------------------------------------------------------------------------------
% load ncbi data
fprintf(1,'Loading ncbi data\n')
load('HomoSampiens_geneInfo20171111.mat'); 

[matches, GeneSymbol, GeneName, nrUpdated, checkEntrezID, MissingProbes] = checkGene(Homosapiens, GeneSymbol, GeneName, EntrezID); 

% save a list of entrezIDs that were not found 
fileID = fopen('entrezIDs.txt','w');
fprintf(fileID, '%d\n',checkEntrezID);
fclose(fileID);

% check if entrezIDs that were not found were discontinued
%% run this from the terminal
% awk 'NR==FNR{pats[$0]; next} $3 in pats' entrezIDs.txt gene_history.txt > listOfentrezIDs.txt
%% load the list 

listOfentrezIDs = readtable('listOfentrezIDs.txt'); 
% name columns according to history file
listOfentrezIDs.Properties.VariableNames = {'tax_id' 'GeneID' 'Discontinued_GeneID' 'Discontinued_Symbol' 'Discontinue_Date'}; 
% remove those that were discontinued and replace those that are changed
% with new ones. 
% find entrezIDs that were discintinued

entrez2change = listOfentrezIDs.Discontinued_GeneID(~strcmp(listOfentrezIDs.GeneID, '-')); 
changeInd = find((ismember(EntrezID,entrez2change)));
fprintf('EntrezIDs for %d probes should be changed\n', length(changeInd));

% first change entrezIDs to updated ones, then keep only the valid ones.
for j=1:length(changeInd)
    entrez = EntrezID(changeInd(j)); 
    index = find(listOfentrezIDs.Discontinued_GeneID==entrez); 
    EntrezID(changeInd(j)) = str2double(listOfentrezIDs.GeneID{index}); 
end

% select entrezIDs to keep (meaning - remove the ones that were
% discontinued)
entrez2rem = listOfentrezIDs.Discontinued_GeneID(strcmp(listOfentrezIDs.GeneID, '-'));
keepInd = find((~ismember(EntrezID,entrez2rem)));
fprintf('Removing %d entrezIDs - %d corresponding probes\n', length(entrez2rem), (length(EntrezID) - length(keepInd)));
% keep only existing ones
ProbeName = ProbeName(keepInd); 
EntrezID = EntrezID(keepInd); 
GeneID = GeneID(keepInd);
GeneSymbol = GeneSymbol(keepInd);
GeneName = GeneName(keepInd);
ProbeID = ProbeID(keepInd);

% run the function again to make changes for new entries. 
[matches2, GeneSymbol, GeneName, nrUpdated, checkEntrezID, MissingProbes] = checkGene(Homosapiens, GeneSymbol, GeneName, EntrezID);
%------------------------------------------------------------------------------
% % then check if entrezID in allen data matches probeIDs compared to agilent
% annotations for agilent probes (CUST probes can't be compared). 
%------------------------------------------------------------------------------
% load annotation file downloaded from
% https://earray.chem.agilent.com/earray/ on 14/11/2017
load('annotations20150612.mat')

% find all probes that are in both files
[probes, iallen, iannot] = intersect(ProbeName, annotations20150612.ProbeID);
% getEntrez allen and annot
entrezAllen = EntrezID(iallen); 
entrezAnnot = annotations20150612.EntrezGeneID(iannot); 

% check if they match
doesMatch = entrezAllen==entrezAnnot; 
% save non-matching in the table
EntrezIDnomatchAllen = entrezAllen(~doesMatch); 
ProbeIDnomatchAllen = probes(~doesMatch);
EntrezIDnomatchAnnot = entrezAnnot(~doesMatch);

% remove those entrez that don't exist in annotation file 
nonmatchingEntrez = table; 
nonmatchingEntrez.probe = ProbeIDnomatchAllen(~isnan(EntrezIDnomatchAnnot)); 
nonmatchingEntrez.entrezAllen = EntrezIDnomatchAllen(~isnan(EntrezIDnomatchAnnot)); 
nonmatchingEntrez.entrezAnnot = EntrezIDnomatchAnnot(~isnan(EntrezIDnomatchAnnot)); 

%------------------------------------------------------------------------------
% Make a table from all the data
%------------------------------------------------------------------------------
DataTable = dataset({Data, headerdata{:}});
% delete probes that were excluded during the quality control from expression and noise matrices
DataTable.Expression{1,1} = DataTable.Expression{1,1}(keepInd,:); DataTable.Noise{1,1} = DataTable.Noise{1,1}(keepInd,:);
DataTable.Expression{2,1} = DataTable.Expression{2,1}(keepInd,:); DataTable.Noise{2,1} = DataTable.Noise{2,1}(keepInd,:);
DataTable.Expression{3,1} = DataTable.Expression{3,1}(keepInd,:); DataTable.Noise{3,1} = DataTable.Noise{3,1}(keepInd,:);
DataTable.Expression{4,1} = DataTable.Expression{4,1}(keepInd,:); DataTable.Noise{4,1} = DataTable.Noise{4,1}(keepInd,:);
DataTable.Expression{5,1} = DataTable.Expression{5,1}(keepInd,:); DataTable.Noise{5,1} = DataTable.Noise{5,1}(keepInd,:);
DataTable.Expression{6,1} = DataTable.Expression{6,1}(keepInd,:); DataTable.Noise{6,1} = DataTable.Noise{6,1}(keepInd,:);
%------------------------------------------------------------------------------
% Remove CUST probes:
%------------------------------------------------------------------------------
if ~useCUSTprobes
    % Remove all CUST probes (assign NaN values for all custom probes)
    fprintf(1,'Removing CUST probes\n')
    cust = strfind(ProbeName, 'CUST');
    remInd = find(~cellfun(@isempty,cust));
    fprintf(1,'%d CUST probes removed\n', length(remInd))
    ProbeName(remInd) = {NaN};
    ProbeID(remInd) = NaN;
    %fprintf(1,'Removing CUST probes\n')
    %fprintf(1,'Removing %d CUST probes\n', sum(isnan(ProbeID)))

ProbeName(isnan(ProbeID)) = [];
EntrezID(isnan(ProbeID)) = [];
GeneID(isnan(ProbeID)) = [];
GeneSymbol(isnan(ProbeID)) = [];
GeneName(isnan(ProbeID)) = [];


DataTable.Expression{1,1}(isnan(ProbeID),:) = []; DataTable.Noise{1,1}(isnan(ProbeID),:) = [];
DataTable.Expression{2,1}(isnan(ProbeID),:) = []; DataTable.Noise{2,1}(isnan(ProbeID),:) = [];
DataTable.Expression{3,1}(isnan(ProbeID),:) = []; DataTable.Noise{3,1}(isnan(ProbeID),:) = [];
DataTable.Expression{4,1}(isnan(ProbeID),:) = []; DataTable.Noise{4,1}(isnan(ProbeID),:) = [];
DataTable.Expression{5,1}(isnan(ProbeID),:) = []; DataTable.Noise{5,1}(isnan(ProbeID),:) = [];
DataTable.Expression{6,1}(isnan(ProbeID),:) = []; DataTable.Noise{6,1}(isnan(ProbeID),:) = [];

ProbeID(isnan(ProbeID)) = [];
end
%------------------------------------------------------------------------------
% Assign ProbeIDs, EntrezIDs and ProbeNames to Data cell.
%------------------------------------------------------------------------------

DataProbe{1,1} = ProbeID;
DataProbe{1,2} = EntrezID;
DataProbe{1,3} = ProbeName;
DataProbe{1,4} = GeneID;
DataProbe{1,5} = GeneSymbol;
DataProbe{1,6} = GeneName;
DataTableProbe = dataset({DataProbe, headerprobe{:}});
%------------------------------------------------------------------------------
% Combine data for all subjects
%------------------------------------------------------------------------------
fprintf(1,'Combining data for all subjects\n')
Expressionall = horzcat(DataTable{1,1}, DataTable{2,1}, DataTable{3,1}, DataTable{4,1}, DataTable{5,1}, DataTable{6,1});
Coordinatesall = vertcat(DataTable{1,2}, DataTable{2,2}, DataTable{3,2}, DataTable{4,2}, DataTable{5,2}, DataTable{6,2});
StructureNamesall = vertcat(DataTable{1,3}, DataTable{2,3}, DataTable{3,3}, DataTable{4,3}, DataTable{5,3}, DataTable{6,3});
MRIvoxCoordinatesAll = vertcat(DataTable{1,4}, DataTable{2,4}, DataTable{3,4}, DataTable{4,4}, DataTable{5,4}, DataTable{6,4});
noiseall = horzcat(DataTable{1,5}, DataTable{2,5}, DataTable{3,5}, DataTable{4,5}, DataTable{5,5}, DataTable{6,5});

%------------------------------------------------------------------------------
% Save relevant variables to a MicroarrayData.mat file
%------------------------------------------------------------------------------
cd ..
cd ('processedData');

fprintf(1,'Saving data to the file\n')
save(sprintf('%sTEST.mat', startFileName), 'DataTable','DataTableProbe', 'Expressionall', 'Coordinatesall', 'StructureNamesall', 'MRIvoxCoordinatesAll', 'noiseall');
cd ../../..
