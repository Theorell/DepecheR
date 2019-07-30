#' @importFrom FNN knnx.index
microClust <- function(dataCenter, dataNeigh, dataReturn, method="median", k=11, trim=0){
  
  
  closest10Pos <- knnx.index(dataNeigh, dataCenter, k=k, algorithm="cover_tree")
  if(method=="median"){
    if(class(dataReturn)=="matrix"){
      closest10Result <- lapply(seq_len(nrow(closest10Pos)), 
                                function(x) 
                                    apply(dataReturn[closest10Pos[x,],], 2, 
                                          median))
      resultDf <- do.call("rbind", closest10Result)
    } else {
      closest10Result <- sapply(seq_len(nrow(closest10Pos)), 
                                function(x) 
                                    median(dataReturn[closest10Pos[x,]]))
    }
    
  } else {
    if(class(dataReturn)=="matrix"){
    closest10Result <- lapply(seq_len(nrow(closest10Pos)), 
                              function(x) 
                                  apply(dataReturn[closest10Pos[x,],], 2, 
                                        mean, trim=trim))
    resultDf <- do.call("rbind", closest10Result)
    } else {
      closest10Result <- sapply(seq_len(nrow(closest10Pos)), 
                                function(x) mean(dataReturn[closest10Pos[x,]], 
                                                 trim=trim))
      
    }  
  }
  
  
}