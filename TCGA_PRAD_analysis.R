# analyzing expression data from TCGA 

# install TCGA bioconductor tools
source("https://bioconductor.org/biocLite.R")
biocLite("TCGAbiolinks")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)

# view data available for prostate cancer (TCGA-PRAD)
TCGAbiolinks:::getProjectSummary("TCGA-PRAD")

## clinical data
# download and parse clinical data directly from TCGA
clinical <- GDCquery_clinic(project = "TCGA-PRAD", type = "clinical")

# inspecting variables of interest
str(clinical) # 500 total records
table(clinical$race) # 147 white, 7 black, 2 asian, 344 not reported
table(clinical$vital_status) # 490 alive, 10 dead
table(clinical$morphology) 
clinical$days_to_death
clinical$bcr_patient_barcode # patient 
clinical$

## FPKM expression data (includes clinical metadata)
## aligned against Hg38

# data preparation and download (if fpkm.RData not in GDCdata/)
# identify desired data
#query_fpkm <- GDCquery(project = "TCGA-PRAD", 
#                      data.category = "Transcriptome Profiling",
#                      data.type = "Gene Expression Quantification", 
#                      workflow.type = "HTSeq - FPKM")
# download data
#GDCdownload(query_fpkm)
# read downloaded data
#fpkm <- GDCprepare(query_fpkm)
# read downloaded data as data frame
#fpkmDF <- GDCprepare(query_fpkm, summarizedExperiment = FALSE)
# save imported object to file
#save(fpkm, file="GDCdata/fpkm.RData")
#save(fpkmDF, file="GDCdata/fpkmDF.RData")

# load saved data
load("GDCdata/fpkm.RData")
load("GDCdata/fpkmDF.RData")
fpkm

# show metadata
colData(fpkm)

# clinical metadata included
subtype_Race
subtype_PSA_preop
subtype_Tumor_cellularity_pathology
subtype_Reviewed_Gleason
subtype_Reviewed_Gleason_category
subtype_Reviewed_Gleason_sum
subtype_Clinical_Gleason_sum

# subset by metadata categories
fpkm[, fpkm$race == "white"]
fpkm[, fpkm$subtype_Race == "WHITE"]

# extract expression data (with seqnames, e.g., ENSG00000000003)
assays(fpkm)$"HTSeq - FPKM"
# show ranges
rowRanges(fpkm)
# show gene names
rowRanges(fpkm)$external_gene_name
# show ensembl_gene_id and external_gene_name
rowData(fpkm)

# TFAM: ENSG00000108064 (ENSG00000108064.9 in DF)

tfam <- fpkmDF %>%
  filter(X1 == "ENSG00000108064.9")
names <- colnames(fpkmDF)



# SPANXB1: ENSG00000227234 (ENSG00000227234.1 in DF)

# FPKM-UQ expression data (includes clinical metadata)
query_fpkmUQ <- GDCquery(project = "TCGA-PRAD", 
                         data.category = "Transcriptome Profiling",
                         data.type = "Gene Expression Quantification",
                         workflow.type = "HTSeq - FPKM-UQ")



## extra
legacy <- GDCquery(project="TCGA-PRAD",
                   data.category="Gene expression",
                   data.type = "Gene Expression Quantification", 
                   platform = "Illumina HiSeq", 
                   file.type  = "normalized_results",
                   experimental.strategy = "RNA-Seq",
                   legacy = TRUE)