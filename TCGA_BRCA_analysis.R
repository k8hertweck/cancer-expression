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
# view untransformed distribution
hist(fpkmGene$SPANXB1) # very left skewed 
hist(fpkmGene$RAC1) # slightly left skewed

# log transform gene expression data
fpkmGene[1:7] <- fpkmGene[1:7] + 1 # add one pseudo count to all counts to remove zeros
fpkmGene[1:7] <- log2(fpkmGene[1:7]) # apply log2 transformation

## visualizing data distribution for variables of interest
# not variable or not reported: classification_of_tumor, last_known_disease_status, tumor_grade, progression_or_recurrence, disease_type
hist(fpkmGene$SPANXB1) # left skewed
hist(fpkmGene$RAC1) # fairly normal, slightly left skewed
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
t.test(SPANXB1 ~ shortLetterCode, data = normVcancer) # p=2.2e-16
table(normVcancer$shortLetterCode) # 113 normal, 1102 tumor
ggplot(normVcancer, aes(shortLetterCode, SPANXB1)) + 
  geom_boxplot() +
  ylab("log2 SPANXB1 expression") +
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
t.test(SPANXB1 ~ shortLetterCode, data = normVcancerPaired) # 9.122e-09
ggplot(normVcancerPaired, aes(shortLetterCode, SPANXB1)) + 
  ylab("log2 SPANXB1 expression") +
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
t.test(SPANXB1 ~ definition, data = tumVmetaAll) # 0.7794
ggplot(tumVmetaAll, aes(definition, SPANXB1)) + 
  ylab("log2 SPANXB1 expression") +
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
t.test(SPANXB1 ~ shortLetterCode, data = TPnoMeta) # p=0.7739
ggplot(TPnoMeta, aes(definition, SPANXB1)) + 
  ylab("log2 SPANXB1 expression") +
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
t.test(SPANXB1 ~ shortLetterCode, data = tumVmetPaired) # 0.4917
ggplot(tumVmetPaired, aes(shortLetterCode, SPANXB1)) + 
  ylab("log2 SPANXB1 expression") +
  xlab("tissue type (paired samples)") +
  scale_x_discrete(labels=c("TM" = "metastatic", "TP" = "tumor")) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
#ggsave("figures/SPANXB1pairedMetastasis.jpg")

## Q3 Compare spanxb1 expression with clinical stages of BC patients
table(fpkmGene$subtype_AJCC.Stage) # missing data are "[Not Available]"
table(fpkmGene$tumor_stage) # missing data are "not reported"
# compare AJCC stages
AJCC <- tum %>%
  filter(subtype_AJCC.Stage != "[Not Available]") %>%
  filter(subtype_AJCC.Stage != "Stage X")
table(AJCC$subtype_AJCC.Stage)
# perform ANOVA
summary(aov(SPANXB1 ~ subtype_AJCC.Stage, data = AJCC)) # 0.907
ggplot(AJCC, aes(subtype_AJCC.Stage, SPANXB1)) + 
  ylab("log2 SPANXB1 expression") +
  xlab("stage") +
  scale_x_discrete(labels=c("Stage I" = "I", "Stage IA" = "IA", "Stage IB" = "IB", "Stage II" = "II", "Stage IIA" = "IIA", "Stage IIB" = "IIB", "Stage III" = "III", "Stage IIIA" = "IIIA", "Stage IIIB" = "IIIB", "Stage IIIC" = "IIIC", "Stage IV" = "IV")) +
  geom_boxplot() +
  theme_bw() 
#ggsave("figures/SPANXB1.AJCC.jpg")
# compare tumor stages
stage <- tum %>%
  filter(tumor_stage != "not reported") %>%
  filter(tumor_stage != "stage x")
table(stage$tumor_stage)
# perform ANOVA
summary(aov(SPANXB1 ~ tumor_stage, data = stage)) #p=0.514
ggplot(stage, aes(tumor_stage, SPANXB1)) + 
  ylab("log2 SPANXB1 expression") +
  xlab("stage") +
  scale_x_discrete(labels=c("stage i" = "I", "stage ia" = "IA", "stage ib" = "IB", "stage ii" = "II", "stage iia" = "IIA", "stage iib" = "IIB", "stage iii" = "III", "stage iiia" = "IIIA", "stage iiib" = "IIIB", "stage iiic" = "IIIC", "stage iv" = "IV")) +
  geom_boxplot() +
  theme_bw() 
#ggsave("figures/SPANXB1.stage.jpg")

## Q4 Compare spanxb1 expression between ER/PR/HER2 positive vs. ER/PR/HER2 negative (a.k.a.TNBC) patients
# extract triple negative from tumor samples
TNBCneg <- tum %>% 
  filter(subtype_ER.Status == "Negative" & subtype_PR.Status == "Negative" & subtype_HER2.Final.Status == "Negative") %>%
  mutate(triple = "negative")
# extract triple positive from tumor samples
TNBCpos <- tum %>% 
  filter(subtype_ER.Status == "Positive" & subtype_PR.Status == "Positive" & subtype_HER2.Final.Status == "Positive") %>%
  mutate(triple = "positive")
# combine triple negative and all positive
TNBCboth <- rbind(TNBCneg, TNBCpos) # 118 TNBC, 59 all positive
# perform t test 
t.test(SPANXB1 ~ triple, data=TNBCboth) # 0.06836
ggplot(TNBCboth, aes(triple, SPANXB1)) + 
  ylab("log2 SPANXB1 expression") +
  xlab("ER/PR/HER2 status") +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
#ggsave("figures/SPANXB1TNBC.jpg")

## Q5 Compare spanxb1 expression between metastatic vs. non metastatic TNBC patients
table(TNBCneg$shortLetterCode) #no TNBC neg or pos are metastatic

## Q6 Compare spanxb1 expression with survival outcome of TNBC patients
table(TNBCneg$vital_status) # 100 alive, 18 dead 
# SPANXB1 and vital status
t.test(SPANXB1 ~ vital_status, data=TNBCneg) # 0.1257
ggplot(TNBCneg, aes(vital_status, SPANXB1)) + 
  ylab("log2 SPANXB1 expression") +
  xlab("vital status") +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
#ggsave("figures/SPANXB1.TNBCneg.vital.jpg")
# SPANXB1 and time to death

## Q7 Compare spanxb1 expression with survival outcome of ER/PR/HER2 positive patients alone and in combination (e.g., co-expression of genes)
table(TNBCpos$vital_status) # 52 alive, 7 dead 
# SPANXB1 and vital status
t.test(SPANXB1 ~ vital_status, data=TNBCpos) # 0.9511
ggplot(TNBCpos, aes(vital_status, SPANXB1)) + 
  ylab("log2 SPANXB1 expression") +
  xlab("vital status") +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
#ggsave("figures/SPANXB1.TNBCpos.vital.jpg")

## Q8 Compare RAC1/SPANXB1 expression together in normal vs. TNBC
# aggregate data from norm and TNBC
normTNBC <- rbind(norm, TNBCneg[,-92]) # removes column triple
# plot both together
ggplot(normTNBC, aes(SPANXB1, RAC1, col=shortLetterCode)) +
  geom_point() +
  ylab("log2 RAC1 expression") +
  xlab("log2 SPANXB1 expression") +
  theme_bw() +
  theme(legend.position="none") + # blue=TNBC, red=normal
  geom_smooth(data=subset(normTNBC, shortLetterCode == "TP"), method = "lm", se = FALSE) +
  geom_smooth(data=subset(normTNBC, shortLetterCode == "NT"), method = "lm", se = FALSE)
#ggsave("figures/SPANXB1.RAC1.TNBC.jpg")
# filter for only expression > 0 for both
bothExp <- normTNBC %>%
  filter(SPANXB1 > 0.000000) %>%
  filter(RAC1 > 0.000000)
ggplot(bothExp, aes(SPANXB1, RAC1, col=shortLetterCode)) +
  geom_point() +
  ylab("log2 RAC1 expression") +
  xlab("log2 SPANXB1 expression") +
  theme_bw() +
  theme(legend.position="none") + # blue=TNBC, red=normal
  geom_smooth(data=subset(bothExp, shortLetterCode == "TP"), method = "lm", se = FALSE) +
  geom_smooth(data=subset(bothExp, shortLetterCode == "NT"), method = "lm", se = FALSE)
#ggsave("figures/SPANXB1.RAC1.TNBC.filtered.jpg")
# logistic regression (unequal groupings)?

## Q9 Compare RAC1/SPANXB1 expression in metastatic vs. not met TNBC
table(TNBCneg$shortLetterCode) #no TNBC neg are metastatic

## Q10 Compare RAC1/SPANXB1 expression with survival of TNBC
TNBCneg
ggplot(TNBCneg, aes(SPANXB1, RAC1, col=vital_status)) +
  geom_point() +
  theme_bw()
# logistic regression (unequal groupings)?
q10model1 <- glm(vital_status ~ SPANXB1, data = TNBCneg, family = binomial())
confint(q10model1, parm = "SPANXB1")
exp(coef(q10model1)["SPANXB1"])
exp(confint(q10model1, parm = "SPANXB1"))
summary(q10model1)
q10model2 <- glm(vital_status ~ SPANXB1 + RAC1, data = TNBCneg, family = binomial())
summary(q10model2)
anova(q10model1, q10model2, test = "Chisq")

## use Wilcoxin signed-rank or Mann-Whitney for small sample sizes?

## correct for multiple comparisons?

# rlog from DESeq2 as alternative transformation method

## permutation tests (requires at least two groups of six samples)
# calculate test statistic for dataset
# permute labels on samples at random 1000 times
# recalculate test statistic for each of 1000 permuted datasets
# calculate percentage of permuted test statistics exceed original
