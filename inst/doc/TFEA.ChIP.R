## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----eval=TRUE,echo=TRUE,message=FALSE----------------------------------------
library(TFEA.ChIP)
library(dplyr)

data( "hypoxia_DESeq", "hypoxia", package="TFEA.ChIP" ) # Load example datasets
hypoxia_df <- preprocessInputData( hypoxia_DESeq )

# Display the first few rows
head(hypoxia_df)

## ----eval=TRUE,echo=TRUE,message=FALSE----------------------------------------
# Preview the raw hypoxia dataset
head(hypoxia)
hypoxia <- preprocessInputData( hypoxia )

# Display the first few rows
head(hypoxia)


## ----eval=TRUE,echo=TRUE,message=FALSE----------------------------------------
library(ExperimentHub)
eh <- ExperimentHub()
ChIPDB <- eh[['EH9854']]  # rE2G300d

## ----eval=TRUE----------------------------------------------------------------
query(eh, "ChIPDBData")

## ----eval=TRUE,echo=TRUE------------------------------------------------------
# Extract vector with names of upregulated genes
Genes.Upreg <- Select_genes( hypoxia_df, min_LFC = 1, max_pval = 0.05 )

# Extract vector with names of non-responsive genes
Genes.Control <- Select_genes( hypoxia_df,
 min_pval = 0.5, max_pval = 1,
 min_LFC = -0.25, max_LFC = 0.25 )

## ----eval=TRUE,echo=TRUE,message=FALSE----------------------------------------
# Conversion of hgnc to ENTREZ IDs
GeneID2entrez( gene.IDs = c("EGLN3","NFYA","ALS2","MYC","ARNT" ) )

# To translate from mouse IDs:
# GeneID2entrez( gene.IDs = c( "Hmmr", "Tlx3", "Cpeb4" ), mode = "m2h" ) # To get the equivalent human gene IDs

## ----eval=TRUE,echo=TRUE------------------------------------------------------
CM_list_UP <- contingency_matrix( Genes.Upreg, Genes.Control ) # Generates list of contingency tables, one per dataset
pval_mat_UP <- getCMstats( CM_list_UP ) # Generates list of p-values and OR from association test
head( pval_mat_UP )

## ----eval=TRUE,echo=TRUE------------------------------------------------------
chip_index <- get_chip_index( TFfilter = c( "HIF1A","EPAS1","ARNT" ) ) # Restrict the analysis to datasets assaying these factors

CM_list_UPe <- contingency_matrix( Genes.Upreg, chip_index = chip_index ) # Generates list of contingency tables
pval_mat_UPe <- getCMstats( CM_list_UPe, chip_index ) # Generates list of p-values and ORs
head( pval_mat_UPe )

## ----eval=TRUE,echo=TRUE------------------------------------------------------
data("MetaData", package = "TFEA.ChIP")
head(MetaData)

## ----eval=FALSE,echo=TRUE-----------------------------------------------------
# pval_mat_UPe <- analysis_from_table(hypoxia_df,
#  method = 'ora', # Overrepresentation Analysis
#  interest_min_LFC = 1, # min LFC
#  expressed = TRUE, # only take expressed TFs
#  control_min_pval = 0.5, control_max_pval = 1, # control pval limits
#  control_min_LFC = -0.25, control_max_LFC = 0.25) # control LFC limits
# 

## ----eval=FALSE, echo=TRUE----------------------------------------------------
# TF_ranking <- rankTFs( pval_mat_UP, rankMethod = "gsea", makePlot = TRUE )
# head( TF_ranking[[ "TF_ranking" ]] )
# TF_ranking[[ "TFranking_plot" ]]

## ----eval=FALSE, echo=TRUE----------------------------------------------------
# library(meta)
# TF_ranking2 <- metaanalysis_fx( pval_mat_UP )
# 
# # Get HIF results
# filter(TF_ranking2$summary,
#   TF %in% c("ARNT", "HIF1A", "EPAS1"))

## ----eval=FALSE, echo=TRUE----------------------------------------------------
# # Result for HIF1A
# forest(TF_ranking2[['results']]$HIF1A)

## ----eval=FALSE, echo=TRUE----------------------------------------------------
# plot_CM( pval_mat_UP ) # plot p-values against ORs

## ----eval=FALSE, echo=TRUE----------------------------------------------------
# HIFs <- c( "EPAS1","HIF1A","ARNT" )
# plot_CM( pval_mat_UP, specialTF = HIFs ) # Plot p-values against ORs highlighting indicated TFs

## ----eval=TRUE,echo=TRUE------------------------------------------------------
chip_index <- get_chip_index( TFfilter = c( "HIF1A","EPAS1","ARNT" ) ) # Restrict the analysis to datasets assaying these factors

## ----eval=TRUE,echo=TRUE,results='hide'---------------------------------------
# run GSEA analysis
GSEA.result <- GSEA_run( hypoxia_df$Genes, hypoxia_df$log2FoldChange, chip_index, get.RES = TRUE) 

## ----eval=FALSE,echo=TRUE-----------------------------------------------------
# GSEA.result <- analysis_from_table(hypoxia_df,
#                                    TFfilter = c('ARNT', 'HIF1A', 'EPAS1'),
#                                    expressed = TRUE,  # Takes only expressed TFs
#                                    method = 'gsea') # Gene Set Enrichment Analysis

## ----eval=FALSE, echo=TRUE----------------------------------------------------
# TF.hightlight <- c( "EPAS1","ARNT","HIF1A" )
# 
# plot_ES( GSEA.result$result, LFC = GSEA.result$processed_table$log2FoldChange, specialTF = TF.hightlight)

## ----eval=FALSE, echo=TRUE----------------------------------------------------
# plot_RES(
#  GSEA_result = GSEA.result$result, LFC = hypoxia_df$log2FoldChange,
#  Accession = c("GSE89836.ARNT.HUVEC-C", "GSE89836.EPAS1.HUVEC-C" ) )

## ----eval=FALSE,echo=TRUE-----------------------------------------------------
# folder <- "~/peak.files.folder"
# File.list<-dir( folder )
# format <- "macs"
# 
# gr.list <- lapply(
#  seq_along( File.list ),
#  function( File.list, myMetaData, format ){
# 
#  tmp<-read.table( File.list[i], ..., stringsAsFactors = FALSE )
# 
#  file.metadata <- myMetaData[ myMetaData$Name == File.list[i], ]
# 
#  ChIP.dataset.gr<-txt2GR(tmp, format, file.metadata)
# 
#  return(ChIP.dataset.gr)
#  },
#  File.list = File.list,
#  myMetadata = myMetadata,
#  format = format
# )

## ----eval=TRUE,echo=TRUE------------------------------------------------------
# As an example of the output
data( "ARNT.peaks.bed","ARNT.metadata", package = "TFEA.ChIP" ) # Loading example datasets for this function
ARNT.gr <- txt2GR( ARNT.peaks.bed, "macs1.4", ARNT.metadata )
head( ARNT.gr, n=2 )

## ----eval=TRUE,echo=TRUE------------------------------------------------------
data( "DnaseHS_db", "gr.list", package="TFEA.ChIP" ) # Loading example datasets for this function
TF.gene.binding.db <- makeChIPGeneDB( DnaseHS_db, gr.list ) 
str( TF.gene.binding.db )

## ----eval=TRUE,echo=TRUE------------------------------------------------------
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
data( "gr.list", package="TFEA.ChIP") # Loading example datasets for this function
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
Genes <- genes( txdb )
TF.gene.binding.db <- makeChIPGeneDB( Genes, gr.list, distanceMargin = 0 )
str( TF.gene.binding.db )

## ----eval=FALSE,echo=TRUE-----------------------------------------------------
# set_user_data( binary_matrix = myTFBSmatrix, metadata = myMetaData )

## ----eval=TRUE,echo=TRUE------------------------------------------------------
sessionInfo()

