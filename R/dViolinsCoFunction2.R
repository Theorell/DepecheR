# Function needed to create the right kind of object for the dViolinsCoFunctionOneClust function
dViolinsCoFunction2 <- function(var, oneClustOneMu, varName, clust, cols, clustNum) {
  allClustOneVar <- data.frame(clust, var)

  colnames(allClustOneVar) <- c("Cluster", "var")

  allClustOneVarOneMuList <- list(allClustOneVar, oneClustOneMu, varName, clustNum, cols)
  return(allClustOneVarOneMuList)
}
