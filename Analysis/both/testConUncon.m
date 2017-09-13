% load group connectome
load('HCP_GenCog_Group.mat')

% load distance uncorrected coexpression
load('averageCoexpressionHCPuncorrected.mat')

Adj = logical(HCP(1:180,1:180)); 
con = averageCoexpression(Adj==1); 
uncon = averageCoexpression(Adj==0); 
dataCell{1} = con; 
dataCell{2} = uncon; 

JitteredParallelScatter(dataCell); 

Adj = logical(HCP(1:180,1:180)); 
con = averageCoexpression(Adj==1); 
uncon = averageCoexpression(Adj==0); 
dataCell{1} = con; 
dataCell{2} = uncon; 

JitteredParallelScatter(dataCell); 

% load distance corrected data
load('geneROI_HCP.mat')
