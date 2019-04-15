# downloading and parsing gene expression data from TCGA-HNSC (head and neck squamous cell carcinoma)

#install.packages("devtools")
#library(devtools)
#devtools::install_github("BioinformaticsFMRP/TCGAbiolinks")
library(TCGAbiolinks)
#source("https://bioconductor.org/biocLite.R")
#BiocInstaller::biocLite("SummarizedExperiment")
library(SummarizedExperiment)
#install.packages("dplyr")
library(dplyr)

# view data available for HNSCC
TCGAbiolinks:::getProjectSummary("TCGA-HNSC")

#### Gene expression ####

## FPKM expression data normalized for upper quartile (includes clinical metadata)
## aligned against Hg38

# data preparation and download (if fpkmUQBca.RData not in GDCdata/)
# identify desired data
query_fpkm <- GDCquery(project = "TCGA-HNSC", 
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification", 
                       workflow.type = "HTSeq - FPKM-UQ")
# download data
GDCdownload(query_fpkm)
# read downloaded data
fpkm <- GDCprepare(query_fpkm)

# save imported object to file
save(fpkm, file="GDCdata/fpkmUQHNSC.RData")

# load saved data
load("GDCdata/fpkmUQHNSC.RData")

# clinical metadata included (all samples, not individual)
table(fpkm$bcr_patient_barcode) # patient identifier
table(fpkm$race) # 468 white, 49 black, 11 asian, 16 not reported
table(fpkm$vital_status) # 294 alive, 252 dead
summary(fpkm$days_to_death) #296 NA
table(fpkm$shortLetterCode) # NT 44, TM 2, TP 500; codes from fpkm$definition

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
# find SDCBP (MDA-9/Syntenin)
mda9 <- grep("SDCBP", genes$external_gene_name, ignore.case = TRUE)
genes[mda9, ]
# SDCBP (MDA-9) ENSG00000137575
# find SIRPA
sirpa <- grep("SIRPA", genes$external_gene_name, ignore.case = TRUE)
genes[sirpa, ]
# SIRPA ENSG00000198053

##  assemble dataset for genes of interest and metadata
fpkmDat <- as.data.frame(t(assays(fpkm)[[1]])) # extract expression data
colnames(fpkmDat) # print gene names
rownames(fpkmDat) # show sample names
# extract gene data for target genes
fpkmGene <- fpkmDat %>%
  select(ENSG00000198021, ENSG00000203926, ENSG00000277215, ENSG00000227234, ENSG00000198573, ENSG00000196406, ENSG00000136238, ENSG00000146648, ENSG00000224057, ENSG00000107295, ENSG00000137575, ENSG00000198053)
# extract metadata
metaDat <-as.data.frame(colData(fpkm))
# bind metadata to gene expression data
fpkmGene <- cbind(fpkmGene, metaDat)
# create object of gene names in order
geneNames <- c("SPANXA1", "SPANXA2", "SPANXA2-OT1", "SPANXB1", "SPANXC", "SPANXD", "RAC1", "EGFR", "EGFR-AS1", "SH3GL2", "MDA9", "SIRPA")
# create object of metadata names
metaNames <- colnames(colData(fpkm))
# replace column names
colnames(fpkmGene) <- c(geneNames, metaNames)

## clean data
# remove troublesome metadata
fpkmGene <- select(fpkmGene, -primary_diagnosis, -treatments, -tissue_or_organ_of_origin, -site_of_resection_or_biopsy, -primary_site, -disease_type)
# save aggregated data to file
write.table(fpkmGene, "targetGeneHNSCC.csv")
