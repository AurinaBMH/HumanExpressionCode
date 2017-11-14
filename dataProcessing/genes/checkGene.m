function [matches, GeneSymbol, GeneName, nrUpdated, checkEntrezID] = checkGene(Homosapiens, GeneSymbol, GeneName, EntrezID)
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
            matches(gene, 1) = length(strfind(symbolNCBIother,symbolAllen));
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
            
            matches(gene, 2) = length(strfind(nameNCBIother,nameAllen));
            
        end
        % make a way to record nr matches updated.
        
        % if one of 2 matches (gene symbol or gene ID) update the other one
        % according to the main name or gene symbol.
        allData = sum(logical(matches(gene,:)));
        if allData == 1
            nrUpdated(l) = EntrezID(gene);
            l=l+1;
            if matches(gene,1)==1 && matches(gene,2)==0
                
                GeneName{gene} = Homosapiens.Full_name_from_nomenclature_authority{gene_ind};
                % update match information as it matches now
                matches(gene,2) = gene;
                
            elseif matches(gene,1)==0 && matches(gene,2)==1
                
                GeneSymbol{gene} = Homosapiens.Symbol{gene_ind};
                % update match information as it matches now
                matches(gene,1) = 1;
                
            end
        end
        
    end
    
end

checkEntrezID = unique(checkEntrezID);
fprintf('%d entrezIDs not found\n', length(checkEntrezID));
% calculate how many don't match after corrections
%fprintf('%d probes initially do not match gene symbol\n', (length(EntrezID) - nansum(matches(:,1))));
%fprintf('%d probes initially do not match gene names\n', (length(EntrezID) - nansum(matches(:,2))));

fprintf('%d probes with updated information according symbol on gene name\n', length(nrUpdated));

fprintf('%d probes STILL do not match gene symbol\n', (length(EntrezID) - nansum(matches(:,1))));
fprintf('%d probes STILL do not match gene names\n', (length(EntrezID) - nansum(matches(:,2))));

end
