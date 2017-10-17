#Mean F-measure calculation. Used in Levine et al, where the data comes from. For calculation procedure, see Critical assessment of automated flow cytometry data analysis techniques, Aghaeepour, N et al nature methods 2013 available at http://www.nature.com/nmeth/journal/v10/n3/full/nmeth.2365.html?foxtrotcallback=true#methods

fMeasureCalculation <- function(clusterVector, idVector, writeFiles=TRUE){
	
	basicMatrix <- table(clusterVector, idVector)
	
	idsTable <- table(idVector)
clusterTable <- table(clusterVector)

recallMatrix <- basicMatrix

for(i in 1:length(idsTable)){
	recallMatrix[,i] <- basicMatrix[,i]/idsTable[i]
}

probabilityMatrix <- basicMatrix

for(i in 1:length(clusterTable)){
	probabilityMatrix[i,] <- basicMatrix[i,]/clusterTable[i]
}

FMatrix <- (2*recallMatrix*probabilityMatrix)/(recallMatrix+probabilityMatrix)

idsFractionVector <- colSums(basicMatrix)/sum(basicMatrix)

FMatrixColMax <- apply(FMatrix, 2, max, na.rm = TRUE)
FMatrixColMaxtimesidsFractionVector <- idsFractionVector*FMatrixColMax

#Here the F-measure is finally calculated. 
fMeasure <- sum(FMatrixColMaxtimesidsFractionVector)
if(writeFiles==TRUE){
	write.csv(recallMatrix, "RecallMatrix.csv")
write.csv(probabilityMatrix, "ProbabilityMatrix.csv")
write.csv(FMatrix, "FMatrix.csv")
}


return(fMeasure)
	
}

