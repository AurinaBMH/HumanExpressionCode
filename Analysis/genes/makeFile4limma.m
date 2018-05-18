% make a matrix for limma normalisation
clear all; 
cd ('data/genes/processedData'); 
%load('MicroarrayDataProbesUpdatedRNAseq360DistThresh2_CoordsAssigned.mat')
load('MicroarrayDataWITHcustProbesUpdatedXXXRNAseq82DistThresh2.mat') % this is better if attempting to 
cort = 1:34; 

for s=1:6
data = DataExpression{s}; 
    select = ismember(data(:,2), cort);
    cortexData = data(select==1,3:end);
    D{s} = cortexData; 
    subjNr{s} = data(select==1,1);
   
end
 
expressionALL = [D{1}; D{2}; D{3}; D{4}; D{5}; D{6}]; 
expressionALL = expressionALL';
batch = [subjNr{1}; subjNr{2}; subjNr{3}; subjNr{4}; subjNr{5}; subjNr{6}]; 

save('limmaExpression.mat', 'expressionALL'); 
save('limmaBatch.mat', 'batch'); 