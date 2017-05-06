# analyzing expression data from TCGA 

# load packages
library(dplyr)
library(ggplot2)

# define exclusion subsetting
`%ni%` <- Negate(`%in%`) 

# read in saved data
fpkmGene <- read.table("targetGeneBca.csv")
# see all metadata
colnames(fpkmGene)

## visualizing data distribution for variables of interest
# not variable or not reported: classification_of_tumor, last_known_disease_status, tumor_grade, progression_or_recurrence, disease_type
hist(fpkmGene$SPANXB1)
plot(fpkmGene$shortLetterCode) # same (abbreviations): plot(fpkmGene$definition)
table(fpkmGene$shortLetterCode) 
# 113 NT= normal tissue, 1102 TP= primary tumor, 7 TM= metastatic 
table(fpkmGene$tumor_stage) # same: table(fpkmGene$subtype_Converted.Stage)
table(fpkmGene$subtype_AJCC.Stage)
table(fpkmGene$vital_status) #table(fpkmGene$subtype_Vital.Status) has missing data
table(fpkmGene$morphology)
table(fpkmGene$subtype_Metastasis) # same: table(fpkmGene$subtype_Metastasis.Coded)

# triple negative indicators
table(fpkmGene$subtype_ER.Status)
table(fpkmGene$subtype_PR.Status)
table(fpkmGene$subtype_HER2.Final.Status)

## initial data subsets
# assess number of samples with SPANXB1 expression
SB1 <- fpkmGene %>%
  filter(SPANXB1 > 0) # 544 samples with SPANXB1 expression
meta <- fpkmGene %>%
  filter(shortLetterCode == "TM")
norm <- fpkmGene %>%
  filter(shortLetterCode == "NT")
tum <- fpkmGene %>%
  filter(shortLetterCode == "TP")

## Q1 Compare spanxb1 expression between normal and BC patients
# all "normal" samples are paired with a Bca sample
normVcancer <- rbind(norm, tum)
# unpaired, all data
t.test(SPANXB1 ~ shortLetterCode, data = normVcancer) # p=0.008919
table(normVcancer$shortLetterCode) # 113 normal, 1102 tumor
ggplot(normVcancer, aes(shortLetterCode, SPANXB1)) + 
  geom_boxplot() +
  ylab("SPANXB1 expression") +
  xlab("tissue type (unpaired)") +
  scale_x_discrete(labels=c("NT" = "normal", "TP" = "tumor")) +
  geom_boxplot() +
  theme_bw() + 
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
#ggsave("figures/SPANXB1unpaired.jpg")
# paired
# create paired sample dataset
normVcancerPaired <- normVcancer %>%
  filter(bcr_patient_barcode %in%
           norm$bcr_patient_barcode)
# find patients with more than 1 tumor sample to pair
normVcancerPaired %>% 
  group_by(patient) %>%
  tally() %>%
  filter(n > 2)
# find patients with only 1 normal sample
normVcancerPaired %>%
  group_by(patient) %>%
  tally() %>%
  filter(n == 1)
# only one sample: TCGA-BH-A0BS 
filter(normVcancerPaired, patient == "TCGA-BH-A0BS")
grep("TCGA-BH-A0BS", normVcancerPaired$patient) # remove 67 (no tumor to pair) 
# three or four samples: TCGA-A7-A0DB, TCGA-A7-A0DC, TCGA-A7-A13E
filter(normVcancerPaired, patient == "TCGA-A7-A0DB")
grep("TCGA-A7-A0DB", normVcancerPaired$patient) # remove 182 and 216 (extra tumor)
filter(normVcancerPaired, patient == "TCGA-A7-A0DC")
grep("TCGA-A7-A0DC", normVcancerPaired$patient) # remove 145 (tumor missing metadata)
filter(normVcancerPaired, patient == "TCGA-A7-A13E")
grep("TCGA-A7-A13E", normVcancerPaired$patient) # remove 221 and 226 (extra tumor)
# remove extra samples
normVcancerPaired <- normVcancerPaired[-c(67, 182, 216, 145, 221, 226),]
# summarize sample counts
table(normVcancerPaired$shortLetterCode) # 112 NT, 112 TP
# perform t test
t.test(SPANXB1 ~ shortLetterCode, data = normVcancerPaired) # 0.0006777
ggplot(normVcancerPaired, aes(shortLetterCode, SPANXB1)) + 
  ylab("SPANXB1 expression") +
  xlab("tissue type (paired samples)") +
  scale_x_discrete(labels=c("NT" = "normal", "TP" = "tumor")) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
#ggsave("figures/SPANXB1paired.jpg")

## Q2 Compare spanxb1 expression between metastatic vs. non metastatic BC patients
# combine tumor and metastasis samples
tumVmetaAll <- rbind(tum, meta)
table(tumVmetaAll$shortLetterCode) # 7 TM 1102 TP
# perform t test (unpaired data)
t.test(SPANXB1 ~ definition, data = tumVmetaAll) # 0.01223
ggplot(tumVmetaAll, aes(definition, SPANXB1)) + 
  ylab("SPANXB1 expression") +
  xlab("tissue type") +
  scale_x_discrete(labels=c("TM" = "metastatic", "TP" = "primary tumor")) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
#ggsave("figures/SPANXB1unpairedMetastasis.jpg")
# remove TP for metastatic patients 
TPnoMeta <- tum %>%
  filter(bcr_patient_barcode %ni%
           meta$bcr_patient_barcode)
TPnoMeta <- rbind(TPnoMeta, meta)
# summarize sample counts
table(TPnoMeta$shortLetterCode) # TM 7, TP 1095
# perform t test (unpaired data, with TP from metastasis patients removed)
t.test(SPANXB1 ~ shortLetterCode, data = TPnoMeta) # 0.01221
ggplot(TPnoMeta, aes(definition, SPANXB1)) + 
  ylab("SPANXB1 expression") +
  xlab("tissue type") +
  scale_x_discrete(labels=c("TM" = "metastatic", "TP" = "primary tumor")) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
#ggsave("figures/SPANXB1unpairedMetastasisFiltered.jpg")
# paired
# create paired metastatic/nonmetastatic sample dataset
tumMet <- tum %>%
  filter(bcr_patient_barcode %in%
           meta$bcr_patient_barcode)
tumVmetPaired <- rbind(tumMet, meta)
# summarize sample counts
table(tumVmetPaired$shortLetterCode) # 7 pairs
# perform t test
t.test(SPANXB1 ~ shortLetterCode, data = tumVmetPaired) # 0.3132
ggplot(tumVmetPaired, aes(shortLetterCode, SPANXB1)) + 
  ylab("SPANXB1 expression") +
  xlab("tissue type (paired samples)") +
  scale_x_discrete(labels=c("TM" = "metastatic", "TP" = "tumor")) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
#ggsave("figures/SPANXB1pairedMetastasis.jpg")

## Q3 Compare spanxb1 expression with clinical stages of BC patients
fpkmGene
table(fpkmGene$subtype_AJCC.Stage)
table(fpkmGene$tumor_stage)

## Q4 Compare spanxb1 expression between ER/PR/HER2 positive vs. ER/PR/HER2 negative (a.k.a.TNBC) patients
# extract triple negative
TNBC <- fpkmGene %>% 
  filter(subtype_ER.Status == "Negative" & subtype_PR.Status == "Negative" & subtype_HER2.Final.Status == "Negative")
# extract normal
norm <- fpkmGene %>%
  filter(shortLetterCode == "NT")
# combine triple negative and normal
TNBCnorm <- rbind(TNBC, norm) # 113 normal, 118 TNBC
t.test(SPANXB1 ~ shortLetterCode, data=TNBCnorm) # 0.003443
ggplot(TNBCnorm, aes(shortLetterCode, SPANXB1)) + 
  ylab("SPANXB1 expression") +
  xlab("tissue type") +
  scale_x_discrete(labels=c("NT" = "normal", "TP" = "TNBC")) +
  geom_boxplot() +
  theme_bw()+
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
#ggsave("figures/SPANXB1TNBC.jpg")

## Q5 Compare spanxb1 expression between metastatic vs. non metastatic TNBC patients

## Q6 Compare spanxb1 expression with survival outcome of TNBC patients
# is gene expression higher in deceased individuals?
vital <- fpkmGene %>%
  filter(shortLetterCode == "TP") %>%
  filter(!is.na(vital_status))
table(vital$vital_status) # 947 alive, 154 dead 
# SPANXB1 and vital status
t.test(SPANXB1 ~ vital_status, data=vital) # p=0.2003
ggplot(vital, aes(vital_status, SPANXB1)) + 
  ylab("SPANXB1 expression") +
  xlab("vital status") +
  geom_boxplot() +
  theme_bw()+
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
#ggsave("figures/SPANXB1vital.jpg")

## Q7 Compare spanxb1 expression with survival outcome of ER/PR/HER2 positive patients alone and in combination

## Q8 Compare RAC1/SPANXB1 expression together in normal vs. TNBC

## Q9 Compare RAC1/SPANXB1 expression in metastatic vs. not met TNBC

## Q10 Compare RAC1/SPANXB1 expression with survival of TNBC

## use Wilcoxin signed-rank or Mann-Whitney for small sample sizes?

## correct for multiple comparisons?

## transform data, or otherwise manage counts of zero?
# log transform data
# rlog from DESeq2

## permutation tests (requires at least two groups of six samples)
# calculate test statistic for dataset
# permute labels on samples at random 1000 times
# recalculate test statistic for each of 1000 permuted datasets
# calculate percentage of permuted test statistics exceed original
