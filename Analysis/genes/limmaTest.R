 
  install.packages("R.matlab")
  library(R.matlab)
  dataOrig <- readMat('limmaExpression.mat')
  dataExp <- as.data.frame(dataOrig) #
  dataExp <- data.matrix(dataExp, rownames.force = NA)
  

  batchNr <- readMat('limmaBatch.mat')
  b<- matrix(unlist(batchNr), ncol = 1, byrow = TRUE)
  
  source("https://bioconductor.org/biocLite.R")
  biocLite()

  
  library(limma) 
  
  normT <- removeBatchEffect(dataExp, b) 
  write.table(normT, file="normalisedExpressionALL.txt", row.names=FALSE, col.names=FALSE)