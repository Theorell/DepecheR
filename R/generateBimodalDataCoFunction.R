generateBimodalDataCoFunction <- function(probabilities, stdevs, observations){

  components <- sample(1:2,prob=probabilities,size=observations,replace=TRUE)

  samples <- rnorm(n=observations,mean=c(0,10)[components],sd=stdevs[components])
  return(samples)

}
