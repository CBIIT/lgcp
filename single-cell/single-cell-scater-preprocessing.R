library(DropletUtils)
library(scater)
library(mbkmeans)
library(scran)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BiocSingular)
library(annotables)
library(tidyverse)

sce <- read10xCounts("/Volumes/group05/CCBB/CS024892_Kelly_Beshiri/02_PrimaryAnalysisOutput/00_FullCellrangerOutputs/SCAF730_190328_G7/outs/filtered_feature_bc_matrix/")

# ADD FEATURE FOR REF ORG HERE!
location_tidy <- rowData(sce) %>% 
  as.data.frame() %>% 
  left_join(grch38 %>% 
              select(ensgene, chr), by = c("ID" = "ensgene")) %>% 
  distinct()
is.mito <- which(location_tidy$chr == "MT")

sce <- calculateQCMetrics(sce, feature_controls=list(Mito=is.mito), BPPARAM = MulticoreParam(7))

# remove cells with log-library sizes that are more than 3 MADs below the median
qc.lib <- isOutlier(log(sce$total_counts), nmads=3, type="lower")
# remove cells where the log-transformed number of expressed genes is 3 MADs below the median
qc.nexprs <- isOutlier(log(sce$total_features_by_counts), nmads=3,
                        type="lower")
# remove cells where the number of mito genes is 3 MADs above the median
qc.mito <- isOutlier(sce$pct_counts_Mito, nmads=3, type="higher")

discard <- qc.lib | qc.nexprs | qc.mito

# Summarize the number of cells removed for each reason.
DataFrame(LibSize=sum(qc.lib),
          NExprs=sum(qc.nexprs),
          MitoProp=sum(qc.mito),
          Total=sum(discard))

# Retain only high-quality cells in the SingleCellExperiment.
sce <- sce[,!discard]

## normalization
sce <- normalize(sce)

## feature selection
fit <- trendVar(sce, use.spikes = FALSE)
dec <- decomposeVar(sce, fit)
hvg <- rownames(dec[dec$bio > 0, ]) # ~4k genes

## dimensionality reduction
set.seed(1234)
sce <- runPCA(sce, 
              feature_set = hvg)

sce <- runUMAP(sce)
