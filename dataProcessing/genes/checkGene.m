function [matches, GeneSymbol, GeneName, nrUpdated, checkEntrezID, MissingProbes] = checkGene(Homosapiens, GeneSymbol, GeneName, EntrezID)
matches = zeros(length(EntrezID),2);
k=1;
l=1;
for gene = 1:length(EntrezID)
    
    gene_ind = find(Homosapiens.GeneID==EntrezID(gene));
    if isempty(gene_ind)
        %fprintf('entrezID %d not found\n', EntrezID(gene))
        matches(gene,1) = NaN;
        matches(gene,2) = NaN;
        % record that entrez ID
        checkEntrezID(k) = EntrezID(gene);
        k=k+1;
    else
        % check if other information matches
        symbolAllen = GeneSymbol{gene};
        symbolAllen(regexp(symbolAllen,',')) = [];
        symbolAllen(regexp(symbolAllen,'-')) = [];
        symbolAllen(regexp(symbolAllen,' ')) = [];
        
        symbolNCBI = Homosapiens.Symbol{gene_ind};
        symbolNCBI(regexp(symbolNCBI,',')) = [];
        symbolNCBI(regexp(symbolNCBI,'-')) = [];
        symbolNCBI(regexp(symbolNCBI,' ')) = [];
        % delete commas, lines and spaces from the gene name
        
        % 1. gene symbol
        matches(gene, 1) = strcmp(symbolNCBI, symbolAllen);
        % 1.2. if a match is not found among symbols, check for all synonims
        
        if matches(gene,1)==0
            symbolNCBIother = [Homosapiens.Synonyms{gene_ind}, Homosapiens.Symbol_from_nomenclature_authority{gene_ind}];
            symbolNCBIother(regexp(symbolNCBIother,',')) = [];
            symbolNCBIother(regexp(symbolNCBIother,'-')) = [];
            symbolNCBIother(regexp(symbolNCBIother,' ')) = [];
            matches(gene, 1) = logical(length(strfind(symbolNCBIother,symbolAllen)));
        end
        
        
        % 2. gene name (check both with comma and witohout)
        nameAllen = GeneName{gene};
        % delete commas, lines and spaces from the gene name
        nameAllen(regexp(nameAllen,',')) = [];
        nameAllen(regexp(nameAllen,'-')) = [];
        nameAllen(regexp(nameAllen,' ')) = [];
        
        
        nameNCBI = Homosapiens.Full_name_from_nomenclature_authority{gene_ind};
        % delete commas, lines and spaces from the gene name
        nameNCBI(regexp(nameNCBI,',')) = [];
        nameNCBI(regexp(nameNCBI,'-')) = [];
        nameNCBI(regexp(nameNCBI,' ')) = [];
        
        nameNCBIother = [Homosapiens.Other_designations{gene_ind}, Homosapiens.description{gene_ind}];
        % delete commas, lines and spaces from the gene name
        nameNCBIother(regexp(nameNCBIother,',')) = [];
        nameNCBIother(regexp(nameNCBIother,'-')) = [];
        nameNCBIother(regexp(nameNCBIother,' ')) = [];
        
        % compare to the database entry
        %1. first try original name
        matches(gene, 2) = strcmp(nameNCBI, nameAllen);
        % if no match, compare to other descriptions
        if matches(gene,2)==0
            
            matches(gene, 2) = logical(length(strfind(nameNCBIother,nameAllen)));
            
        end
        % make a way to record nr matches updated.
        
        % if one of 2 matches (gene symbol or gene ID) update the other one
        % according to the main name or gene symbol.
        allData = sum(logical(matches(gene,:)));
        if allData == 1
            nrUpdated(l,1) = EntrezID(gene);
            nrUpdated(l,2) = gene; 
            l=l+1;
            if matches(gene,1)==1 && matches(gene,2)==0
                
                GeneName{gene} = Homosapiens.Full_name_from_nomenclature_authority{gene_ind};
                % update match information as it matches now
                matches(gene,2) = 1;
                
            elseif matches(gene,1)==0 && matches(gene,2)==1
                
                GeneSymbol{gene} = Homosapiens.Symbol{gene_ind};
                % update match information as it matches now
                matches(gene,1) = 1;
                
            end
        end
        
    end
    
end
doesExistCheck = exist('checkEntrezID');
if doesExistCheck
    checkEntrezID = unique(checkEntrezID);
    fprintf('%d entrezIDs not found\n', length(checkEntrezID));
else
    checkEntrezID = [];
end
% calculate how many don't match after corrections
%fprintf('%d probes initially do not match gene symbol\n', (length(EntrezID) - nansum(matches(:,1))));
%fprintf('%d probes initially do not match gene names\n', (length(EntrezID) - nansum(matches(:,2))));
doesExist = exist('nrUpdated');
if doesExist
    fprintf('%d probes with updated information according symbol or gene name\n', length(nrUpdated));
else
    fprintf('NO probes with updated information \n');
    nrUpdated = [];
end

onlyExisting = matches;
E = EntrezID;
onlyExisting(any(isnan(matches), 2), :) = [];
E(any(isnan(matches), 2), :) = [];

ind = sum(onlyExisting,2)==0;
MissingProbes = E(ind);

fprintf('%d probes where entrezID exists, but STILL do not match gene symbol\n', length(find(matches(:,1)==0)));
fprintf('%d probes where entrezID exists, but STILL do not match gene symbol\n', length(find(matches(:,2)==0)))

end
