#' @importFrom gplots heatmap.2
#' @importFrom stats as.dendrogram order.dendrogram hclust dist
# This function is part of the depeche function. See depeche and dAllocate for 
# most flags. 
# inDataFrameUsed: the scaled indata, where a zoom-in has been made based on
# the samplingSubset.
# reducedClusterCenters: the selected cluster centers, where only variables
# that have not been completely sparsed out are kept. 
depecheHeatMap <- function(inDataFrameUsed, reducedClusterCenters, 
                           plotDir, logCenterSd, createOutput){
    # Here, a heatmap over the cluster
    # centers is saved. Only true if the
    # number of clusters exceeds one.
    if (ncol(reducedClusterCenters) > 500) {
        message("As the number of variables in the result exceeds 500, ",
                "it is not meaningful to produce a cluster center heatmap, ",
                "so it is omitted")
    } else if (nrow(reducedClusterCenters) > 1 &&
               ncol(reducedClusterCenters) > 1) {
        graphicClusterCenters <- reducedClusterCenters
        # Here we scale each center value to the
        # range between the lowest and highest
        # permille of the observations in the
        # inDataScaled for that variable
        for (i in seq_len(ncol(graphicClusterCenters))) {
            scaledFocus <-
                inDataFrameUsed[, colnames(inDataFrameUsed) ==
                                    colnames(graphicClusterCenters)[i]]
            graphicClusterCenters[, i] <- dScale(reducedClusterCenters[, i],
                                                 scaledFocus,
                                                 robustVarScale = FALSE,
                                                 center = FALSE,
                                                 truncate = TRUE)
        }
        
        # Here, the order of the columns are changed, so that the
        # most similar ones are the closest to each other
        colOrder <-
            order.dendrogram(
                as.dendrogram(hclust(dist(t(graphicClusterCenters)))))
        
        graphicClusterCenters[reducedClusterCenters == 0] <- NA
        graphicClusterCenters <- graphicClusterCenters[, colOrder]
        
        colorLadder <- dColorVector(seq_len(11),
                                    colorScale = c(
                                        "#0D0887FF", "#6A00A8FF", "#900DA4FF",
                                        "#B12A90FF", "#CC4678FF", "#E16462FF", 
                                        "#F1844BFF", "#FCA636FF", "#FCCE25FF"
                                    )
        )
        
        if (createOutput) {
            # First, the name is defined, depending on two criteria: if the data
            # was log transformed, and if a directory should be created or not.
            logOrNoLogName <- ifelse(logCenterSd[[1]] == FALSE, "Cluster",
                                     "Log2_transformed_cluster"
            )
            clusterCenterName <- file.path(plotDir, paste0(
                logOrNoLogName,
                "_centers.pdf"
            ))
            pdf(clusterCenterName)
            heatmap.2(as.matrix(graphicClusterCenters),
                      Rowv = FALSE,
                      Colv = FALSE, dendrogram = "none", scale = "none",
                      col = colorLadder, breaks = seq(0, 1, length.out = 12),
                      trace = "none", keysize = 1.5, density.info = "none",
                      key.xlab =
                          "0=no expression, 1=high expression\ngrey=penalized",
                      na.color = "#A2A2A2"
            )
            dev.off()
        }
    }
}
