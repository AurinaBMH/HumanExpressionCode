function toCheckOriginal = selectDuplicates(IDone, IDtwo, changetoNumber)
% 

if nargin>2
    for i=1:length(changetoNumber)
        IDone{changetoNumber(i)} = '999999';
    end
end
[~, ind] = unique(IDone, 'stable');
% find duplicate indices
duplicate_ind = setdiff(1:size(IDone,1), ind);

k=1;
toCheckOriginal = cell(1,1);
for gene = 1:length(duplicate_ind)
    
    symb = IDone(duplicate_ind(gene));
    if isa(IDone, 'double')
        test_ind = find(IDone==symb);
    elseif isa(IDone, 'cell')
        test_ind = find(strcmp(IDone, symb{1}));
    end

    % check if for all duplicated instances gene name is the same
    geneID = IDtwo(test_ind);
    doMatch = zeros(length(test_ind));
    for j=1:length(test_ind)
        for i=1:length(test_ind)
            
            doMatch(j,i) = geneID(j)~=geneID(i);
            
        end
    end
    % if they're not matching, record indexes to exclude later
    if sum(doMatch(:))~=0
        toCheckOriginal{k,:} = test_ind;
        k=k+1;
    end
end

% remove duplicated rows in the cell

if ~isempty(toCheckOriginal)
    [r1, r2] = ndgrid(1:size(toCheckOriginal, 1));
    duplicates = any(triu(arrayfun(@(r1, r2) isequal(toCheckOriginal(r1, :), toCheckOriginal(r2, :)), r1, r2), 1));
    toCheckOriginal(duplicates, :) = [];
end

end

