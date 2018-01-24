load('MicroarrayDataProbesUpdatedRNAseq360DistThresh2_CoordsAssigned.mat')

expression = vertcat(DataExpression{1},DataExpression{2},DataExpression{3}, ...
    DataExpression{4},DataExpression{5},DataExpression{6}); 

% the the manuscript reported on non-normalised data
% probably better calculate on normalised, because we don't want to average raq values ecross subjects
% (effect stronger) for normalised
expData = BF_NormalizeMatrix(expression(:,3:end),'scaledRobustSigmoid');

rois = unique(expression(:,2)); 
l=1; k=1; 

for r=1:length(rois)
    
    ind = find(expression(:,2)==rois(r));
    values = expData(ind,:);
    if length(ind)>1
        m = mean(values);
    else
        m = values;
    end
    if rois(r)<=180
        Vl(k) = var(m);
        k=k+1;
    elseif rois(r)>180
        Vr(l) = var(m);
        l=l+1;
    end
    
end

[p,h,stats] = ranksum(Vl,Vr); 