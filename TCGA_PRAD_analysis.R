# analyzing expression data from TCGA 

# install TCGA bioconductor tools
source("https://bioconductor.org/biocLite.R")
biocLite("TCGAbiolinks")
library(TCGAbiolinks)

# view data available for prostate cancer (TCGA-PRAD)
TCGAbiolinks:::getProjectSummary("TCGA-PRAD")

## clinical data
# download and parse clinical data directly from TCGA
clinical <- GDCquery_clinic(project = "TCGA-PRAD", type = "clinical")

# inspecting variables of interest
str(clinical) # 500 total records
table(clinical$race) # 147 white, 7 black, 2 asian, 344 not reported

## expression data
GDCquery(project = "TCGA-PRAD", data.category = "Transcriptome Profiling")
