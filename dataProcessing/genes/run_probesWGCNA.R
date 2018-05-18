# select probes WGCNA
install.packages("R.matlab")
library(R.matlab)
library(WGCNA)
  

data = readMat('probes2WGCNA.mat')
datET = as.data.frame(data$datET); 

rowGroup = data$rowGroup; 
rowID = unlist(data$rowID); 

row.names(datET)<-rowID

outdatED = collapseRows(datET, rowGroup, rowID,
             method="MaxMean", connectivityBasedCollapsing=FALSE,
             methodFunction=NULL, connectivityPower=1,
             selectFewestMissing=TRUE, thresholdCombine=NA)

