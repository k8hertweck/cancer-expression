#### analyzing expression data from TCGA ####

# load packages
library(dplyr)
library(ggplot2)

# define exclusion subsetting
`%ni%` <- Negate(`%in%`) 

# read in saved data
fpkmGene <- read.table("targetGeneHNSCC.csv")
# see all metadata
colnames(fpkmGene)
# view untransformed distribution
hist(fpkmGene$MDA9) # very slightly left skewed
hist(fpkmGene$SIRPA) # slightly left skewed
# save untransformed data
fpkmGeneNolog <- fpkmGene

# log transform genes
fpkmGene[1:12] <- fpkmGene[1:12] + 1 # add one pseudo count to all counts to remove zeros
fpkmGene[1:12] <- log2(fpkmGene[1:12]) # apply log2 transformation
hist(fpkmGene$SIRPA) 
hist(fpkmGene$MDA9) 

#### MDA-9/Syntenin ####

# high expression in AA is correlated with poor outcome when compared to non-hispanic whites


# high expression in HIS is correlated with poor outcome when compared to non-hispanic whites




#### SIRPA ####

# high expression in AA is correlated with poor outcome when compared to non-hispanic whites
