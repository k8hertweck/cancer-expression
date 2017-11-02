# downloading and parsing mutation data from TCGA-PRAD (prostate cancer)

#install.packages("devtools")
#library(devtools)
#devtools::install_github("BioinformaticsFMRP/TCGAbiolinks")
library(TCGAbiolinks)
#source("https://bioconductor.org/biocLite.R")
biocLite("maftools") # requires R 3.4 or higher
library(maftools)
library(dplyr)
library(stringr)

# view data available for prostate cancer (TCGA-PRAD)
TCGAbiolinks:::getProjectSummary("TCGA-PRAD")

#### Clinical data #### 
# download and parse clinical data 
clinical <- GDCquery_clinic(project = "TCGA-PRAD", type = "clinical")

# inspecting variables of interest
str(clinical) # 500 total records
table(clinical$race) # 147 white, 7 black, 2 asian, 344 not reported
table(clinical$vital_status) # 490 alive, 10 dead
table(clinical$morphology) 
clinical$days_to_death
clinical$bcr_patient_barcode # patient 

#### Simple Nucelotide Variation ####

# data preparation and download (if .RData not in GDCdata/)
# identify samples and download
query_maf <- GDCquery_Maf("PRAD", pipelines = "muse")

# create maf object (without clinical data)
maf <- read.maf(query_maf, useAll = FALSE)
# visual summary of data
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
# draw oncoplot
oncoplot(maf = maf, top = 10, removeNonMutated = TRUE)
# assess transitions and transversions
maf_titv <- titv(maf = maf, plot = FALSE, useSyn = TRUE)
# plot transitions and transversions
plotTiTv(res = maf_titv)

# extract barcodes from maf query
barcodes <- sort(query_maf$Tumor_Sample_Barcode) %>%
  unique
str_trunc(barcodes, 12, side = "right", ellipsis = "")

# extract barcodes from clinical data
sub_id <- sort(clinical$submitter_id) %>%
  unique()
clinical <- rename(clinical, Tumor_Sample_Barcode = submitter_id)

# create maf object with clinical data attached
maf <- read.maf(query_maf, clinicalData = clinical, useAll = FALSE)
