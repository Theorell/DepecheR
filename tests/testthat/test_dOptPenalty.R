######################
#dOptPenalty
context('dOptPenatly')
x <- DepecheR:::generateBimodalData()[,-1]
x_scaled <- DepecheR:::dScale(x)
dOptPenaltyResult <- DepecheR:::dOptPenalty(x_scaled, k=30, maxIter=20, bootstrapObservations=500, penalties=c(2,4), makeGraph=TRUE, disableWarnings=TRUE, minARI=0.95)