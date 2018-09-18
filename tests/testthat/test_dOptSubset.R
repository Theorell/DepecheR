######################
#dOptSubset
context('dOptSubset')
x <- DepecheR:::generateBimodalData()[,-1]
x_scaled <- DepecheR:::dScale(x)
depecheResult <- DepecheR:::dOptSubset(inDataFrameScaled=x_scaled, firstClusterNumber=1, sampleSizes=500, k=30, maxIter=20, minARI=0.95, minARIImprovement=0.01, penalties=c(2,4), selectionSampleSize=500)
