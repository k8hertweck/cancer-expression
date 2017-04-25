# analyzing expression data from TCGA 

# load packages
library(dplyr)
library(ggplot2)

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
