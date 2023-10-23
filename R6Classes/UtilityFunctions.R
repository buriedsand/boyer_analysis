ensureDirExists <- function(directory) {
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
}

processPeakAnnotations <- function(peakAnno,
                                   coveragePlotDir,
                                   annotationsDir,
                                   featureDistPlotDir,
                                   featureSummaryFile,
                                   subsets,
                                   subsetsDir,
                                   peakAnnoRDSFile) {
  # Ensure directories exist
  ensureDirExists(coveragePlotDir)
  ensureDirExists(annotationsDir)
  ensureDirExists(featureDistPlotDir)
  ensureDirExists(dirname(featureSummaryFile))
  ensureDirExists(subsetsDir)
  ensureDirExists(dirname(peakAnnoRDSFile))
  
  # Create coverage plots
  peakAnno$createCoveragePlots(directory = coveragePlotDir)
  
  # Annotate peaks
  peakAnno$annotateFiles()
  
  # Save annotations as dataframes
  lapply(names(peakAnno$annotations), function(name) {
    filepath <- file.path(annotationsDir, paste0(name, ".csv"))
    write.csv(as.data.frame(peakAnno$annotations[[name]]),
              filepath,
              row.names = FALSE)
  })
  
  # Create feature distribution plots
  peakAnno$createFeatureDistPlots(directory = featureDistPlotDir)
  
  # Create and save feature summary
  feature_summary <- peakAnno$summarizeAll()
  write.csv(feature_summary, featureSummaryFile, row.names = FALSE)
  
  # Extract and save specified subsets
  for (subset_name in subsets) {
    current_subsetDir <- file.path(subsetsDir, subset_name)
    ensureDirExists(current_subsetDir)
    peakAnno$subsetAnnotations(subset_name, directory = current_subsetDir)
  }
  
  # Save peakAnno object as RDS
  saveRDS(peakAnno, peakAnnoRDSFile)
}

# Example usage:
# processPeakAnnotations(
#   peakAnno = peakAnno,
#   coveragePlotDir = "path/to/coveragePlots",
#   annotationsDir = "path/to/annotations",
#   featureDistPlotDir = "path/to/featureDistPlots",
#   featureSummaryFile = "path/to/featureSummary.csv",
#   subsets = c("Promoter", "distal_intergenic"),
#   subsetsDir = "path/to/subsets",
#   peakAnnoRDSFile = "path/to/peakAnno.rds"
# )
