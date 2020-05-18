#' @importFrom beanplot beanplot
#' @importFrom reshape2 melt
# This function is called by dViolins. It provides the core functionality,
# and is called for each cluster separately.
dViolinsCoFunction <- function(clusterNum, plotElement, clusterVector,
                               inDataFrameTrunc, clusterCol, plotDir,
                               createOutput) {
    # First, the plot elements of interest are selected
    inDataFrameRelevant <- inDataFrameTrunc[, plotElement]

    # In the case that this turns the dataframe into a vector, this is corrected
    # here for systematic reasons.
    if (is.vector(inDataFrameRelevant)) {
        inDataFrameRelevant <- as.data.frame(inDataFrameRelevant)
    }

    # Here, a splitting factor, turning all non-relevant clusters into 0 is
    # created
    clusterVector[clusterVector != clusterNum] <- 0

    # Now, the dataframe is split based on this
    inDataFrameSplit <- split(inDataFrameRelevant, f = clusterVector)

    # Now, if any of the two dataframes contain more than 5000 events, a random
    # subsample is used downstream
    for (i in seq_len(2)) {
        if (nrow(inDataFrameSplit[[i]]) > 5000) {
            inDataFrameSplit[[i]] <-
                inDataFrameSplit[[i]][sample(
                    seq_len(nrow(inDataFrameSplit[[i]])), 5000
                ), ]
        }
    }

    # And now, again, to avoid errors due to the exception with only one
    # variable, a new small transformaiton is made in that case
    if (ncol(inDataFrameRelevant) == 1) {
        for (i in seq_len(2)) {
            inDataFrameSplit[[i]] <- as.data.frame(inDataFrameSplit[[i]])
            colnames(inDataFrameSplit[[i]]) <- plotElement
        }
    }

    inDataFrameMelt <- lapply(inDataFrameSplit, melt)

    for (i in seq(1, 2)) {
        inDataFrameMelt[[i]]$variable <- paste0(
            inDataFrameMelt[[i]]$variable,
            " ", names(inDataFrameMelt[i])
        )
    }

    # And now, the data frames are merged again
    inDataFrameUnSplit <- do.call("rbind", inDataFrameMelt)

    # And so it is all plotted
    if (createOutput) {
        pdf(file.path(plotDir, paste0("Cluster_", clusterNum, "_violins.pdf")))

        par(las = 3, mar = c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
        beanplot(value ~ variable,
            data = inDataFrameUnSplit, what = c(0, 1, 0, 0),
            main = paste0(
                "Feature distribution for cluster ",
                clusterNum, " vs all others"
            ),
            ylab = "Expression level, arbitrary units", side = "both",
            border = NA, col = list("grey", clusterCol), maxwidth = 0.8
        )
        legend("topright",
            inset = c(-0.33, 0), legend = c(
                paste0("All - ", clusterNum),
                clusterNum
            ),
            fill = c("grey", clusterCol), title = "Cluster"
        )

        dev.off()
    }
}
