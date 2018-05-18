
load('MicroarrayDataWITHcustProbesUpdatedXXXmaxCorrelationnoQCPearson.mat')
expression{1} = vertcat(expressionAll{1}, expressionAll{2}, expressionAll{3}, ...
    expressionAll{4}, expressionAll{5},expressionAll{6});

% make matrix from cell
varargout = csvimport('MyDataWGCNA.csv'); 
expression = varargout(2:end, 2:end);
expression2 = cell2mat(expression);
entrezID = varargout(:,1); 

i=1; 
entrezIDWGCNA = zeros(20232,1); 

for g=2:20233
    st = entrezID{g}; 
    entrezIDWGCNA(i) = str2double(regexp(st,'\d*','match'));
    i=i+1; 
end

entrezIDWGCNA = cell2mat(entrezID);
[entrezIDWGCNAsor, indWGCNA] = sort(entrezIDWGCNA, 'ascend'); 
% reorder expression values

expressionWGCNA = expression2(indWGCNA,:); 
expressionWGCNA = expressionWGCNA'; 

r = zeros(size(expression{1},2),1); 
for p1=1:size(expression{1},2)
    r(p1) = corr(expression{1}(:,p1), expressionWGCNA(:,p1));     
end


