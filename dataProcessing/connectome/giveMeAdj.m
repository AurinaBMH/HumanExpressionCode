
allMatrrices = zeros(380,380,123); 
for sub = 1:100
    allMatrices(:,:,sub) = logical(standard{sub}); 
end
consist = mean(allMatrices,3); 

thresholds = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]; 
for j=1:length(thresholds)

keep = consist>thresholds(j); 

Adj = double(keep); 
Adj = Adj(1:180,1:180); 
con = averageCoexpression(Adj==1); 
uncon = averageCoexpression(Adj==0); 
dataCell{1} = con; 
dataCell{2} = uncon; 

JitteredParallelScatter(dataCell); 
end

figure; imagesc(Adj); 

subjects = [10 20 30 40 50 60 70 80 90 100]; 
for i=1:length(subjects)
    Adj = logical(standard{subjects(i)}(1:180,1:180)); 
    [S,P] = RichClubHuman(Adj,averageCoexpression);
end