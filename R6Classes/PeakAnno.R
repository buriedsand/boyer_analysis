library(R6)
library(ChIPseeker)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(pbapply)

PeakAnno <- R6Class("PeakAnno",
    public = list(
        dirPath = NULL,
        bedFiles = NULL,
        annotations = NULL,
        coveragePlots = NULL,
        featureDistPlots = NULL,
        summaryFeatureDist = NULL,

        initialize = function(dirPath = "") {
            self$dirPath <- dirPath
            self$bedFiles <- list.files(dirPath, pattern="*.bed$", full.names=TRUE)
        },

        createCoveragePlots = function(directory = NULL) {
            self$coveragePlots <- pblapply(self$bedFiles, function(file) {
                peak <- readPeakFile(file)
                plot <- covplot(peak, weightCol="V5")
                if (!is.null(directory)) {
                    filepath <- file.path(directory, paste0(basename(file), "_covplot.png"))
                    ggsave(filepath, plot=plot)
                }
                return(plot)
            })
        },

        annotateFiles = function() {
            txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
            self$annotations <- pblapply(self$bedFiles, function(file) {
                annotatePeak(file, tssRegion=c(-5000, 5000), TxDb=txdb, annoDb="org.Hs.eg.db")
            }, cl = NULL)
        },

        subsetAnnotations = function(feature, directory = NULL) {
            subset_data <- lapply(self$annotations, function(ann) {
                feature_idx <- ann@detailGenomicAnnotation[[feature]]
                df <- as.data.frame(ann)[feature_idx, c("seqnames", "start", "end", "V4", "V5", "V6")]
                df$start <- df$start - 1
                return(df)
            })
            if (!is.null(directory)) {
                lapply(seq_along(subset_data), function(i) {
                    current_name <- names(subset_data)[[i]]
                    current_data <- subset_data[[i]]

                    filepath <- file.path(directory, paste0(current_name, ".bed"))
                    write.table(
                        current_data,
                        filepath,
                        quote=FALSE,
                        sep="\t",
                        row.names=FALSE,
                        col.names=FALSE,
                    )
                })
            }
            return(subset_data)
        },

        createFeatureDistPlots = function(directory = NULL) {
            self$featureDistPlots <- lapply(seq_along(self$annotations), function(i) {
                current_name <- names(self$annotations)[[i]]
                current_ann <- self$annotations[[i]]

                plot <- plotAnnoBar(current_ann)
                if (!is.null(directory)) {
                    filepath <- file.path(directory, paste0(current_name, "_featureDist.png"))
                    ggsave(filepath, plot=plot)
                }
                return(plot)
            })
        },

        summarizeAll = function() {
            df_list <- list()
            
            for (i in seq_along(self$annotations)) {
                current_name <- names(self$annotations)[i]
                current_ann <- self$annotations[i]
                
                current_df <- current_ann[[current_name]]@annoStat
                current_df <- data.frame(
                    name=current_name,
                    Feature=current_df$Feature,
                    Frequency=current_df$Frequency
                )
                
                df_list[[i]] <- current_df
            }
            
            # Combine all dataframes vertically
            do.call(rbind, df_list)
        }
    )
)

# Usage example:
# processor <- PeakAnno$new("/path/to/bedfiles/")
# processor$annotateFiles()
# processor$createCoveragePlots(directory = "/path/to/coverage_plots/")
# promoters <- processor$subsetAnnotations("Promoter")
# enhancers <- processor$subsetAnnotations("Enhancer")
# processor$createFeatureDistPlots(directory = "/path/to/featureDist_plots/")
