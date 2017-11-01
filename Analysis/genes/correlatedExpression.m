
% test the fact that expression is correlated between brains as stated in
% allens 2012 paper. select regions mutual for both brains and compare
% expression values (all regions, one gene at the time for 10000 (or all) genes) 

load('MicroarrayDataWITHcustMean82DistThresh2_CoordsAssigned.mat')
% mean for region
for i=1:length(unique(DataExpression{1}(:,2)))
a = find(DataExpression{1}(:,2)==i); 
brain1(i,:) = mean(DataExpression{1}(a,:),1); 
end
for i=1:length(unique(DataExpression{2}(:,2)))
b = find(DataExpression{2}(:,2)==i);
brain2(i,:) = mean(DataExpression{2}(b,:),1);
end
% select regions that have data in both brains
[both, ind1, ind2] = intersect(brain1(:,2), brain2(:,2)); 
final1 = brain1(ind1,:);
final2 = brain2(ind2,:);
% select only expression values
brainA = final1(:,3:end); 
brainB = final2(:,3:end);

% plot
figure;
for i=1:size(brainA,2)
    scatter(mean(brainA(:,i),1), mean(brainB(:,i),1)); 
    hold on; 
end