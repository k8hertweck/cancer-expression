# downloading and parsing data from TCGA 

# install TCGA bioconductor tools
source("https://bioconductor.org/biocLite.R")
biocLite("TCGAbiolinks")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)

# view data available for prostate cancer (TCGA-PRAD)
TCGAbiolinks:::getProjectSummary("TCGA-PRAD")

#### Clinical data #### 
# incomplete compared to metadata with Gene Expression?)
# download and parse clinical data directly from TCGA
clinical <- GDCquery_clinic(project = "TCGA-PRAD", type = "clinical")

# inspecting variables of interest
str(clinical) # 500 total records
table(clinical$race) # 147 white, 7 black, 2 asian, 344 not reported
table(clinical$vital_status) # 490 alive, 10 dead
table(clinical$morphology) 
clinical$days_to_death
clinical$bcr_patient_barcode # patient 

#### Gene expression ####

## FPKM expression data normalized for upper quartile (includes clinical metadata)
## aligned against Hg38

# data preparation and download (if fpkm.RData not in GDCdata/)
# identify desired data
query_fpkm <- GDCquery(project = "TCGA-PRAD", 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - FPKM-UQ")
# download data
GDCdownload(query_fpkm)
# read downloaded data
fpkm <- GDCprepare(query_fpkm)
# save imported object to file
save(fpkm, file="GDCdata/fpkmUQPca.RData")

# read downloaded data as data frame (adds numbers to gene names)
#fpkmDF <- GDCprepare(query_fpkm, summarizedExperiment = FALSE)
#fpkmDFnoUQ.RData and fpkmNoUQ.RData are from non-normalized data (not included in git repo, code not shown)

# load saved data
load("GDCdata/fpkmUQPca.RData")
fpkm <- fpkmUQ
rm(fpkmUQ)

# clinical metadata included (all samples, not individual)
table(fpkm$bcr_patient_barcode) # patient identifier
table(fpkm$race) # 190 white, 11 black, 2 asian, 348 not reported
table(fpkm$vital_status) # 541 alive, 10 dead
fpkm$days_to_death
table(fpkm$morphology) 
table(fpkm$shortLetterCode) # 498 TP (Primary solid Tumor), 52 NT (Solid Tissue Normal), 1 TM (Metastatic); codes from fpkm$definition
table(fpkm$subtype_Race)
fpkm$subtype_Tumor_cellularity_pathology
table(fpkm$subtype_Reviewed_Gleason)
table(fpkm$subtype_Reviewed_Gleason_category)
table(fpkm$subtype_Reviewed_Gleason_sum) # prefer to subtype_Clinical_Gleason_sum
table(fpkm$subtype_Clinical_Gleason_sum) 

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
# find TFAM
tfam <- grep("TFAM$", genes$external_gene_name, perl = TRUE)
genes[tfam, ]
# TFAM = ENSG00000108064

## assemble dataset for genes of interest and metadata
fpkmDat <- as.data.frame(t(assays(fpkm)[[1]])) # extract expression data
colnames(fpkmDat) # print gene names
rownames(fpkmDat) # show sample names
# extract gene data for target genes
fpkmGene <- fpkmDat %>%
  select(ENSG00000108064)
# extract metadata
metaDat <-as.data.frame(colData(fpkm))
# bind metadata to gene expression data
fpkmGene <- cbind(fpkmGene, metaDat)
# create object of gene names in order
geneNames <- "TFAM"
# create object of metadata names
metaNames <- colnames(colData(fpkm))
# replace column names
colnames(fpkmGene) <- c(geneNames, metaNames)

## clean data
# normalize subtype_Race colum
fpkmGene$subtype_Race <- gsub("WHITE", "white", fpkmGene$subtype_Race)
fpkmGene$subtype_Race <- gsub("BLACK_OR_AFRICAN_AMERICAN", "black or african american", fpkmGene$subtype_Race)
fpkmGene$subtype_Race <- gsub("BLACK OR AFRICAN AMERICAN", "black or african american", fpkmGene$subtype_Race)
fpkmGene$subtype_Race <- gsub("ASIAN", "asian", fpkmGene$subtype_Race)
# replace "not reported" in race column with NA
fpkmGene$race[fpkmGene$race == "not reported"] <- NA
# aggregate race columns 
fpkmGene$race <- ifelse(is.na(fpkmGene$race), fpkmGene$subtype_Race, fpkm$race)
# remove extraneous race column
fpkmGene <- select(fpkmGene, -subtype_Race)
# check race column stats
table(fpkmGene$race)
# force race and Gleason to factor
fpkmGene$race <- as.factor(fpkmGene$race)
fpkmGene$subtype_Reviewed_Gleason_sum <- as.factor(fpkmGene$subtype_Reviewed_Gleason_sum)
# complete metadata for: TCGA-HC-8260 (black), TCGA-HC-8262 (white)
fpkmGene$race <- ifelse(fpkmGene$bcr_patient_barcode == "TCGA-HC-8260", "black or african american", as.character(fpkmGene$race))
fpkmGene$race <- ifelse(fpkmGene$bcr_patient_barcode == "TCGA-HC-8262", "white", as.character(fpkmGene$race))
# remove troublesome metadata
fpkmGene <- select(fpkmGene, -treatments)
# save aggregated data to file
write.table(fpkmGene, "targetGenePca.csv")

#### Simple Nucelotide Variation ####