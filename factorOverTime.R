library(ChIPseeker)
library(ggplot2)

# Utils ====
generatePeakAnnoList <- function(objectA, objectB, factorName) {
  # Filter annotations based on the presence of factorName in the names
  filtered_A <-
    objectA$annotations[grep(factorName, names(objectA$annotations))]
  filtered_B <-
    objectB$annotations[grep(factorName, names(objectB$annotations))]
  
  # Combine the filtered annotations
  peakAnnoList <- unlist(list(filtered_A, filtered_B))
  
  return(peakAnnoList)
}

# Main script ====

peakAnno_3d <- readRDS("results/3_5_days_regular/peakAnno.rds")
peakAnno_4h <- readRDS("results/4_5_hours_regular/peakAnno.rds")

factor <- "KLF2"
peakAnnoList <-
  generatePeakAnnoList(peakAnno_4h, peakAnno_3d, factor)
plot <- plotAnnoBar(peakAnnoList)
ggsave(
  filename = paste0(factor, "_plot.png"),
  plot = plot,
  device = "png",
  width = 8,
  height = 6
)
