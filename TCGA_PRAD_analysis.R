# analyzing expression data from TCGA 

# install TCGA bioconductor tools
source("https://bioconductor.org/biocLite.R")
biocLite("TCGAbiolinks")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)

# view data available for prostate cancer (TCGA-PRAD)
TCGAbiolinks:::getProjectSummary("TCGA-PRAD")

## clinical data (incomplete compared to metadata with Gene Expression?)
# download and parse clinical data directly from TCGA
clinical <- GDCquery_clinic(project = "TCGA-PRAD", type = "clinical")

# inspecting variables of interest
str(clinical) # 500 total records
table(clinical$race) # 147 white, 7 black, 2 asian, 344 not reported
table(clinical$vital_status) # 490 alive, 10 dead
table(clinical$morphology) 
clinical$days_to_death
clinical$bcr_patient_barcode # patient 

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
# read downloaded data as data frame (adds numbers to gene names)
#fpkmDF <- GDCprepare(query_fpkm, summarizedExperiment = FALSE)
# save imported object to file
#save(fpkm, file="GDCdata/fpkm.RData")
#save(fpkmDF, file="GDCdata/fpkmDF.RData")

# load saved data
load("GDCdata/fpkm.RData")
load("GDCdata/fpkmDF.RData")
fpkm

# clinical metadata included (all samples, not individual)
table(fpkm$bcr_patient_barcode) # patient identifier
table(fpkm$race) # 190 white, 11 black, 2 asian, 348 not reported
table(fpkm$vital_status) # 541 alive, 10 dead
fpkm$days_to_death
table(fpkm$morphology) 
table(fpkm$shortLetterCode) # TP (Primary solid Tumor), NT (Solid Tissue Normal), TM (Metastatic); codes from fpkm$definition
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
# show gene names
rowRanges(fpkm)$external_gene_name
# show ensembl_gene_id and external_gene_name
rowData(fpkm)

##  assemble dataset for genes of interest and metadata
fpkmDat <- as.data.frame(t(assays(fpkm)[[1]])) # extract expression data
colnames(fpkmDat) # print gene names
rownames(fpkmDat) # show sample names
# extract gene data for target genes
# TFAM: ENSG00000108064 
# SPANXB1: ENSG00000227234
fpkmGene <- fpkmDat %>%
  select(ENSG00000108064, ENSG00000227234)
# find average sample-wide expression
colMeans(fpkmGene)
# bind metadata to gene expression data
fpkmGene <- cbind(fpkmGene, fpkm$bcr_patient_barcode, fpkm$race, fpkm$vital_status, fpkm$days_to_death, fpkm$morphology, fpkm$shortLetterCode, fpkm$subtype_Race, fpkm$subtype_Tumor_cellularity_pathology, fpkm$subtype_Reviewed_Gleason, fpkm$subtype_Reviewed_Gleason_category, fpkm$subtype_Reviewed_Gleason_sum)
# replace column names
colnames(fpkmGene) <- c("TFAM", "SPANXB1", "bcr_patient_barcode", "race", "vital_status", "days_to_death", "morphology", "shortLetterCode", "subtype_Race", "subtype_Tumor_cellularity_pathology", "subtype_Reviewed_Gleason", "subtype_Reviewed_Gleason_category", "subtype_Reviewed_Gleason_sum")

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
# force race and Gleason to factor
fpkmGene$race <- as.factor(fpkmGene$race)
fpkmGene$subtype_Reviewed_Gleason_sum <- as.factor(fpkmGene$subtype_Reviewed_Gleason_sum)

# incomplete metadata for: TCGA-HC-8260 (black), TCGA-HC-8262 (white)
# add extra metadata for genetic variants

# save aggregated data to file
write.csv(fpkmGene, "targetGeneData.csv")

# read in saved data
fpkmGene <- read.csv("targetGeneData.csv")
fpkmGene <- select(fpkmGene, -subtype_Race)
# force race and Gleason to factor
fpkmGene$race <- as.factor(fpkmGene$race)
fpkmGene$subtype_Reviewed_Gleason_sum <- as.factor(fpkmGene$subtype_Reviewed_Gleason_sum)  

## visualizing data distribution
hist(fpkmGene$TFAM)
hist(fpkmGene$SPANXB1)
plot(fpkmGene$race)
plot(fpkmGene$morphology)
plot(fpkmGene$shortLetterCode)
plot(fpkmGene$subtype_Tumor_cellularity_pathology)
plot(fpkmGene$subtype_Reviewed_Gleason)
plot(fpkmGene$subtype_Reviewed_Gleason_category)
plot(fpkmGene$subtype_Reviewed_Gleason_sum)

## statistical testing

# does gene expression differ between normal and prostate?
normVcancer <- filter(fpkmGene, shortLetterCode != "TM") # remove metastasis
# TFAM unpaired, all data
t.test(TFAM ~ shortLetterCode, data = normVcancer) # p=0.2911
ggplot(normVcancer, aes(shortLetterCode, TFAM)) + 
  geom_boxplot()
# SPANXB1 unpaired, all data
t.test(SPANXB1 ~ shortLetterCode, data = normVcancer) # p=0.06523
ggplot(normVcancer, aes(shortLetterCode, SPANXB1)) + 
  geom_boxplot()

# paired samples
normID <- normVcancer %>% # list normal samples
  filter(shortLetterCode == "NT") %>%
  select(bcr_patient_barcode)
normVcancerPaired <- normVcancer %>%
  filter(bcr_patient_barcode %in%
           normID$bcr_patient_barcode)
# three samples (remove one): TCGA-HC-7740, TCGA-HC-8258

# does gene expression differ between white and AA?
whiteVaa <- fpkmGene %>%
  filter(!is.na(fpkmGene$race)) %>% # remove missing data
  filter(race != "asian") # remove asian
# TFAM and race
t.test(TFAM ~ race, data=whiteVaa) # p=0.7693
ggplot(whiteVaa, aes(race, TFAM)) + 
  geom_boxplot()
# SPANXB1 and race
t.test(SPANXB1 ~ race, data=whiteVaa) # p=0.8854
ggplot(whiteVaa, aes(race, SPANXB1)) + 
  geom_boxplot()

# is gene expression higher in deceased individuals?
# TFAM and vital status
t.test(TFAM ~ vital_status, data=fpkmGene) # p=0.1704
ggplot(fpkmGene, aes(vital_status, TFAM)) + 
  geom_boxplot()
# SPANXB1 and vital status
t.test(SPANXB1 ~ vital_status, data=fpkmGene) # p=0.4386
ggplot(fpkmGene, aes(vital_status, SPANXB1)) + 
  geom_boxplot()

# is gene expression related to tumor progression?
gleason <- fpkmGene %>%
  filter(!is.na(subtype_Reviewed_Gleason_sum))
# TFAM anova for Gleason
summary(aov(TFAM ~ subtype_Reviewed_Gleason_sum, dat=gleason)) #p=2.49e-05
ggplot(gleason, aes(subtype_Reviewed_Gleason_sum, TFAM)) + 
  geom_boxplot()
ggsave("figures/TFAM_GleasonSum.pdf")
summary(aov(TFAM ~ subtype_Reviewed_Gleason, dat=gleason)) #p=0.00199
ggplot(gleason, aes(subtype_Reviewed_Gleason, TFAM)) + 
  geom_boxplot()
ggsave("figures/TFAM_Gleason.pdf")
summary(aov(TFAM ~ morphology, dat=gleason)) #p=2.95e-08
ggplot(gleason, aes(morphology, TFAM)) + 
  geom_boxplot()
summary(aov(TFAM ~ subtype_Tumor_cellularity_pathology, dat=gleason)) #p=0.00876
ggplot(gleason, aes(subtype_Tumor_cellularity_pathology, TFAM)) + 
  geom_boxplot()
# SPANXB1 anova Gleason sum
summary(aov(SPANXB1 ~ subtype_Reviewed_Gleason_sum, dat=gleason)) #p=0.735
ggplot(gleason, aes(subtype_Reviewed_Gleason_sum, SPANXB1)) + 
  geom_boxplot()
summary(aov(SPANXB1 ~ subtype_Reviewed_Gleason, dat=gleason)) #p=0.878
ggplot(gleason, aes(subtype_Reviewed_Gleason, SPANXB1)) + 
  geom_boxplot()
summary(aov(SPANXB1 ~ morphology, dat=gleason)) #p=2.09e-07
ggplot(gleason, aes(morphology, SPANXB1)) + 
  geom_boxplot() # 8255/3 only has two data points
summary(aov(SPANXB1 ~ subtype_Tumor_cellularity_pathology, dat=gleason)) #p=0.591
ggplot(gleason, aes(subtype_Tumor_cellularity_pathology, SPANXB1)) + 
  geom_boxplot() 

# is gene expression related to genetic variants(especially SPANXB1)?

###########

# FPKM-UQ expression data (includes clinical metadata)
query_fpkmUQ <- GDCquery(project = "TCGA-PRAD", 
                         data.category = "Transcriptome Profiling",
                         data.type = "Gene Expression Quantification",
                         workflow.type = "HTSeq - FPKM-UQ")
# download data
#GDCdownload(query_fpkmUQ)
# read downloaded data
#fpkmUQ <- GDCprepare(query_fpkmUQ)
# save imported object to file
#save(fpkmUQ, file="GDCdata/fpkmUQ.RData")

# load saved data
load("GDCdata/fpkmUQ.RData")
fpkmUQ

# clinical metadata included (all samples, not individual)
fpkmUQ$bcr_patient_barcode # patient 
table(fpkmUQ$race) # 190 white, 11 black, 2 asian, 348 not reported
table(fpkmUQ$vital_status) # 541 alive, 10 dead
fpkmUQ$days_to_death
table(fpkmUQ$morphology) 
table(fpkmUQ$shortLetterCode) # TP (Primary solid Tumor), NT (Solid Tissue Normal), TM (Metastatic); codes from fpkmUQ$definition
table(fpkmUQ$subtype_Race)
fpkmUQ$subtype_Tumor_cellularity_pathology
table(fpkmUQ$subtype_Reviewed_Gleason)
table(fpkmUQ$subtype_Reviewed_Gleason_category)
table(fpkmUQ$subtype_Reviewed_Gleason_sum) # prefer to subtype_Clinical_Gleason_sum
table(fpkmUQ$subtype_Clinical_Gleason_sum) 

# inspect object
assayNames(fpkmUQ)
head(assay(fpkmUQ), 1)
colSums(assay(fpkmUQ))
rowRanges(fpkmUQ)
colData(fpkmUQ) # metadata
# show gene names
rowRanges(fpkmUQ)$external_gene_name
# show ensembl_gene_id and external_gene_name
rowData(fpkmUQ)

# extract gene data for all 
fpkmUQDat <- as.data.frame(t(assays(fpkmUQ)[[1]]))
colnames(fpkmUQDat) # print gene names
rownames(fpkmUQDat) # show sample names
# extract gene data for target genes
# TFAM: ENSG00000108064 
# SPANXB1: ENSG00000227234
fpkmUQGene <- fpkmUQDat %>%
  select(ENSG00000108064, ENSG00000227234)
# find average sample-wide expression
colMeans(fpkmUQGene)
# bind metadata to gene expression data
fpkmUQGene <- cbind(fpkmUQGene, fpkmUQ$bcr_patient_barcode, fpkmUQ$race, fpkmUQ$vital_status, fpkmUQ$days_to_death, fpkmUQ$morphology, fpkmUQ$shortLetterCode, fpkmUQ$subtype_Race, fpkmUQ$subtype_Tumor_cellularity_pathology, fpkmUQ$subtype_Reviewed_Gleason, fpkmUQ$subtype_Reviewed_Gleason_category, fpkmUQ$subtype_Reviewed_Gleason_sum)
# replace column names
colnames(fpkmUQGene) <- c("TFAM", "SPANXB1", "bcr_patient_barcode", "race", "vital_status", "days_to_death", "morphology", "shortLetterCode", "subtype_Race", "subtype_Tumor_cellularity_pathology", "subtype_Reviewed_Gleason", "subtype_Reviewed_Gleason_category", "subtype_Reviewed_Gleason_sum")
