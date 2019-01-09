# This function is the core of this whole thing, creating the plot
dViolinsCoFunction3 <- function(allClustOneVarOneMuList, plotAll, createOutput) {
  if (plotAll == FALSE) {
    if (allClustOneVarOneMuList[[2]] != 0) {
      plotname <- paste("Cluster", "_", allClustOneVarOneMuList[[4]], "_", allClustOneVarOneMuList[[3]], ".pdf", sep = "")
      dp <- ggplot(allClustOneVarOneMuList[[1]], aes(x = Cluster, y = var, fill = Cluster)) +
        geom_violin(trim = FALSE) +
        scale_color_manual(values = allClustOneVarOneMuList[[5]]) +
        scale_fill_manual(values = allClustOneVarOneMuList[[5]]) +
        labs(title = paste("Plot of", allClustOneVarOneMuList[[3]], "for cluster", allClustOneVarOneMuList[[4]]), x = "Cluster", y = "Intensity")
      dp + theme_classic()
      if (createOutput == TRUE) {
        ggsave(filename = plotname, dpi = 300)
      }
    }
  } else {
    plotname <- paste("Cluster", "_", allClustOneVarOneMuList[[4]], "_", allClustOneVarOneMuList[[3]], ".pdf", sep = "")
    dp <- ggplot(allClustOneVarOneMuList[[1]], aes(x = Cluster, y = var, fill = Cluster)) +
      geom_violin(trim = FALSE) +
      scale_color_manual(values = allClustOneVarOneMuList[[5]]) +
      scale_fill_manual(values = allClustOneVarOneMuList[[5]]) +
      labs(title = paste("Plot of", allClustOneVarOneMuList[[3]], "for cluster", allClustOneVarOneMuList[[4]]), x = "Cluster", y = "Intensity")
    dp + theme_classic()
    if (createOutput == TRUE) {
      ggsave(filename = plotname, dpi = 300)
    }
  }
}
