# analyzing expression data from TCGA 

# load packages
library(dplyr)
library(ggplot2)

# read in saved data
fpkmGene <- read.table("targetGeneData.csv")
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

# does gene expression differ between normal and prostate? (unpaired)
normVcancer <- filter(fpkmGene, shortLetterCode != "TM") # remove metastasis
# TFAM unpaired, all data
t.test(TFAM ~ shortLetterCode, data = normVcancer) # p=0.2191
table(normVcancer$shortLetterCode) # 52 normal, 498 tumor
ggplot(normVcancer, aes(shortLetterCode, TFAM)) + 
  geom_boxplot() +
  ylab("TFAM expression") +
  xlab("tissue type (unpaired)") +
  scale_x_discrete(labels=c("NT" = "normal", "TP" = "tumor")) +
  geom_boxplot() +
  theme_bw() + 
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
ggsave("figures/TFAMunpaired.jpg")
# SPANXB1 unpaired, all data
t.test(SPANXB1 ~ shortLetterCode, data = normVcancer) # p=0.0524
SB1 <- normVcancer %>%
  filter(SPANXB1 > 0) # 19 samples with SPANXB1 expression, all tumor
ggplot(normVcancer, aes(shortLetterCode, SPANXB1)) + 
  geom_boxplot() +
  ylab("SPANXB1 expression") +
  xlab("tissue type (unpaired)") +
  scale_x_discrete(labels=c("NT" = "normal", "TP" = "tumor")) +
  geom_boxplot() +
  theme_bw() + 
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
ggsave("figures/SPANXB1unpaired.jpg")

# does gene expression differ between normal and prostate? (paired)
normID <- normVcancer %>% # list normal samples
  filter(shortLetterCode == "NT") %>%
  select(bcr_patient_barcode)
normVcancerPaired <- normVcancer %>%
  filter(bcr_patient_barcode %in%
           normID$bcr_patient_barcode)
# three samples (remove one): TCGA-HC-7740, TCGA-HC-8258
filter(normVcancerPaired, patient == "TCGA-HC-7740")
grep("TCGA-HC-7740", normVcancerPaired$patient) # 5 (tumor) 11 (tumor missing metadata) 89 (normal)
filter(normVcancerPaired, patient == "TCGA-HC-8258")
grep("TCGA-HC-8258", normVcancerPaired$patient) # 10 (tumor) 60 (tumor missing metadata) 91 (normal)
# remove extra samples
normVcancerPaired <- normVcancerPaired[-c(11, 60),]
# summarize sample counts
table(normVcancerPaired$shortLetterCode) # 52 NT, 52 TP
# TFAM paired
t.test(TFAM ~ shortLetterCode, data = normVcancerPaired, paired=TRUE) # p=0.0009166
t.test(TFAM ~ shortLetterCode, data = normVcancerPaired) # p=0.0008068
ggplot(normVcancerPaired, aes(shortLetterCode, TFAM)) + 
  ylab("TFAM expression") +
  xlab("tissue type (paired samples)") +
  scale_x_discrete(labels=c("NT" = "normal", "TP" = "tumor")) +
  geom_boxplot() +
  theme_bw() + 
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
ggsave("figures/TFAMpaired.jpg")
# SPANXB1 paired
t.test(SPANXB1 ~ shortLetterCode, data = normVcancerPaired, paired=TRUE) # p=0.04633
t.test(SPANXB1 ~ shortLetterCode, data = normVcancerPaired) # p=0.04633
normVcancerPaired %>%
  filter(SPANXB1 > 0) # only four samples with expression
ggplot(normVcancerPaired, aes(shortLetterCode, SPANXB1)) + 
  ylab("SPANXB1 expression") +
  xlab("tissue type (paired samples)") +
  scale_x_discrete(labels=c("NT" = "normal", "TP" = "tumor")) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
ggsave("figures/SPANXB1paired.jpg")

# add SPANXA2.OT1

# does gene expression differ between white and AA?
whiteVaa <- fpkmGene %>%
  filter(!is.na(fpkmGene$race)) %>% # remove missing data
  filter(race != "asian") %>% # remove asian
  filter(shortLetterCode == "TP") # only tumor
table(whiteVaa$race) # 23 AA, 175 white
# TFAM and race
t.test(TFAM ~ race, data=whiteVaa) # p=0.3886
ggplot(whiteVaa, aes(race, TFAM)) + 
  ylab("TFAM expression") +
  xlab("race") +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
ggsave("figures/TFAMrace.jpg")
# SPANXB1 and race
t.test(SPANXB1 ~ race, data=whiteVaa) # p=0.8875
ggplot(whiteVaa, aes(race, SPANXB1)) + 
  ylab("SPANXB1 expression") +
  xlab("race") +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
ggsave("figures/SPANXB1race.jpg")

# is gene expression higher in deceased individuals?
vital <- fpkmGene %>%
  filter(shortLetterCode == "TP")
table(vital$vital_status) # 488 alive, 10 dead (all dead are tumor)
# TFAM and vital status 
t.test(TFAM ~ vital_status, data=vital) # p=0.2257
ggplot(vital, aes(vital_status, TFAM)) + 
  ylab("TFAM expression") +
  xlab("vital status") +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
ggsave("figures/TFAMvital.jpg")
# SPANXB1 and vital status
t.test(SPANXB1 ~ vital_status, data=vital) # p=0.4691
ggplot(vital, aes(vital_status, SPANXB1)) + 
  ylab("SPANXB1 expression") +
  xlab("vital status") +
  geom_boxplot() +
  theme_bw()+
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
ggsave("figures/SPANXB1vital.jpg")

# is gene expression related to tumor progression?
gleason <- fpkmGene %>%
  filter(!is.na(subtype_Reviewed_Gleason_sum))
# summarize counts
table(gleason$subtype_Reviewed_Gleason_sum) # 6=62, 7=172, 8=44, 9=41, 10=1 person
table(gleason$subtype_Reviewed_Gleason)
table(gleason$morphology)
# TFAM 
summary(aov(TFAM ~ subtype_Reviewed_Gleason_sum, dat=gleason)) #p=0.000164
ggplot(gleason, aes(subtype_Reviewed_Gleason_sum, TFAM)) + 
  ylab("TFAM expression") +
  xlab("Gleason score") +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
ggsave("figures/TFAM_GleasonSum.jpg")
summary(aov(TFAM ~ subtype_Reviewed_Gleason, dat=gleason)) #p=0.0105
ggplot(gleason, aes(subtype_Reviewed_Gleason, TFAM)) + 
  ylab("TFAM expression") +
  xlab("Gleason score") +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
ggsave("figures/TFAM_Gleason.jpg")
summary(aov(TFAM ~ morphology, dat=gleason)) #p=1.09e-07
ggplot(gleason, aes(morphology, TFAM)) + 
  geom_boxplot()
summary(aov(TFAM ~ subtype_Tumor_cellularity_pathology, dat=gleason)) #p=0.00391
ggplot(gleason, aes(subtype_Tumor_cellularity_pathology, TFAM)) + 
  geom_boxplot()
summary(aov(TFAM ~ subtype_Subtype, dat=gleason)) #p=0.818
ggplot(gleason, aes(subtype_Subtype, TFAM)) + 
  geom_boxplot()
summary(aov(TFAM ~ subtype_Residual_tumor, dat=gleason)) #p=0.313; should remove NA
ggplot(gleason, aes(subtype_Residual_tumor, TFAM)) + 
  geom_boxplot()
summary(aov(TFAM ~ subtype_SCNA_cluster, dat=gleason)) #p=0.000859
ggplot(gleason, aes(subtype_SCNA_cluster, TFAM)) + 
  geom_boxplot()

# SPANXB1 
summary(aov(SPANXB1 ~ subtype_Reviewed_Gleason_sum, dat=gleason)) #p=0.716
ggplot(gleason, aes(subtype_Reviewed_Gleason_sum, SPANXB1)) + 
  geom_boxplot()
summary(aov(SPANXB1 ~ subtype_Reviewed_Gleason, dat=gleason)) #p=0.863
ggplot(gleason, aes(subtype_Reviewed_Gleason, SPANXB1)) + 
  geom_boxplot()
summary(aov(SPANXB1 ~ morphology, dat=gleason)) #p=3.54e-07
ggplot(gleason, aes(morphology, SPANXB1)) + 
  geom_boxplot() # 8255/3 only has two data points
summary(aov(SPANXB1 ~ subtype_Tumor_cellularity_pathology, dat=gleason)) #p=0.591
ggplot(gleason, aes(subtype_Tumor_cellularity_pathology, SPANXB1)) + 
  geom_boxplot() 
summary(aov(SPANXB1 ~ subtype_Subtype, dat=gleason)) #p=0.851
ggplot(gleason, aes(subtype_Subtype, SPANXB1)) + 
  geom_boxplot()
summary(aov(SPANXB1 ~ subtype_Residual_tumor, dat=gleason)) #p=0.315; should remove NA
ggplot(gleason, aes(subtype_Residual_tumor, SPANXB1)) + 
  geom_boxplot()
summary(aov(SPANXB1 ~ subtype_SCNA_cluster, dat=gleason)) #p=0.237
ggplot(gleason, aes(subtype_SCNA_cluster, SPANXB1)) + 
  geom_boxplot()

# is gene expression related to other expression or genetic variants?

## TFAM
# ETV1, ETV4

## SPANXB1
# methylation
# ETV1, FLI1, FOXA1, MED12, ADH1. KMT2A, KMT2C
