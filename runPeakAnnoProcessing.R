library(magrittr)
source("R6Classes/UtilityFunctions.R")
source("R6Classes/PeakAnno.R")

processExperiment <- function(experimentName) {
  basePath <- paste0("results/", experimentName)
  
  peakAnno <- PeakAnno$new(file.path(basePath, "02_peak_calling/called_peaks"))
  
  # Add names to the bedFiles
  names(peakAnno$bedFiles) <- 
    tools::file_path_sans_ext(basename(peakAnno$bedFiles)) %>%
    sub(".seacr.peaks.stringent$", "", .)
  
  # Process peakAnno object
  processPeakAnnotations(
    peakAnno = peakAnno,
    coveragePlotDir = file.path(basePath, "02_peak_calling/coverage_plots"),
    annotationsDir = file.path(basePath, "03_peak_annotations/all_annotations"),
    featureDistPlotDir = file.path(basePath, "03_peak_annotations/feature_distribution_plots"),
    featureSummaryFile = file.path(basePath, "03_peak_annotations/feature_summary.csv"),
    subsets = c("Promoter", "distal_intergenic"),
    subsetsDir = file.path(basePath, "03_peak_annotations/subsetted_peaks"),
    peakAnnoRDSFile = file.path(basePath, "peakAnno.rds")
  )
}

# Usage:
experimentNames <- c("3_5_days_regular", "4_5_hours_regular")
lapply(experimentNames, processExperiment)
