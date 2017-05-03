# downloading and parsing data from TCGA 

# install TCGA bioconductor tools
source("https://bioconductor.org/biocLite.R")
biocLite("TCGAbiolinks")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)

# view data available for prostate cancer (TCGA-PRAD)
TCGAbiolinks:::getProjectSummary("TCGA-BRCA")

## clinical data 
# download and parse clinical data directly from TCGA
clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")

# inspecting variables of interest
str(clinical) # 500 total records
table(clinical$race) # XX white, X black, X asian, X not reported
table(clinical$vital_status) # X alive, X dead
table(clinical$morphology) 
clinical$days_to_death
clinical$bcr_patient_barcode # patient 

## FPKM expression data normalized for upper quartile (includes clinical metadata)
## aligned against Hg38

# data preparation and download (if fpkmUQBca.RData not in GDCdata/)
# identify desired data
query_fpkm <- GDCquery(project = "TCGA-BRCA", 
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification", 
                       workflow.type = "HTSeq - FPKM-UQ")
# download data
GDCdownload(query_fpkm)
# read downloaded data
fpkm <- GDCprepare(query_fpkm)

# save imported object to file
save(fpkm, file="GDCdata/fpkmUQBca.RData")

# load saved data
load("GDCdata/fpkmUQBca.RData")

# clinical metadata included (all samples, not individual)
table(fpkm$bcr_patient_barcode) # patient identifier
table(fpkm$race) # X white, X black, X asian, X not reported
table(fpkm$vital_status) # X alive, X dead
fpkm$days_to_death
table(fpkm$morphology) 
table(fpkm$shortLetterCode) # X; codes from fpkm$definition

# inspect object
assayNames(fpkm)
head(assay(fpkm), 1)
colSums(assay(fpkm))
rowRanges(fpkm)
colData(fpkm) # metadata
colnames(colData(fpkm)) # just metadata column names
# show gene names
rowRanges(fpkm)$external_gene_name

# extract genes of interest
# create object of ensembl_gene_id and external_gene_name
genes <-rowData(fpkm)
# find SPANX
spanx <- grep("spanx", genes$external_gene_name, ignore.case = TRUE)
genes[spanx, ]
# SPANXA1 ENSG00000198021
# SPANXA2 ENSG00000203926
# SPANXA2-OT1 ENSG00000277215
# SPANXB1 ENSG00000227234
# SPANXC ENSG00000198573
# SPANXD ENSG00000196406

##  assemble dataset for genes of interest and metadata
fpkmDat <- as.data.frame(t(assays(fpkm)[[1]])) # extract expression data
colnames(fpkmDat) # print gene names
rownames(fpkmDat) # show sample names
# extract gene data for target genes
fpkmGene <- fpkmDat %>%
  select(ENSG00000198021, ENSG00000203926, ENSG00000277215, ENSG00000227234, ENSG00000198573, ENSG00000196406)
# extract metadata
metaDat <-as.data.frame(colData(fpkm))
# bind metadata to gene expression data
fpkmGene <- cbind(fpkmGene, metaDat)
# create object of gene names in order
geneNames <- c("SPANXA1", "SPANXA2", "SPANXA2-OT1", "SPANXB1", "SPANXC", "SPANXD")
# create object of metadata names
metaNames <- colnames(colData(fpkm))
# replace column names
colnames(fpkmGene) <- c(geneNames, metaNames)

## clean data
# remove troublesome metadata
fpkmGene <- select(fpkmGene, -treatments)
# save aggregated data to file
write.table(fpkmGene, "targetGeneBca.csv")
