#### analyzing expression data from TCGA ####

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
hist(fpkmGene$EGFR) # very left skewed
hist(fpkmGene$EGFR.AS1) # very left skewed
hist(fpkmGene$SH3GL2) # very left skewed
# save untransformed data
fpkmGeneNolog <- fpkmGene

# log transform gene expression data
fpkmGene[1:10] <- fpkmGene[1:10] + 1 # add one pseudo count to all counts to remove zeros
fpkmGene[1:10] <- log2(fpkmGene[1:10]) # apply log2 transformation

# add SH3GL2/SPANXB1 expression ratio
fpkmGene <- fpkmGene %>%
  mutate(SH3GL2_SPANXB1_ratio = SH3GL2 / SPANXB1)

## visualizing data distribution for variables of interest
# not variable or not reported: classification_of_tumor, last_known_disease_status, tumor_grade, progression_or_recurrence, disease_type
hist(fpkmGene$SPANXB1) # left skewed (many zeros)
hist(fpkmGene$RAC1) # fairly normal, slightly left skewed
hist(fpkmGene$EGFR) # normal
hist(fpkmGene$EGFR.AS1) # left skewed (many zeros)
hist(fpkmGene$SH3GL2) # left skewed (many zeros)
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

#### Q1 Compare spanxb1 expression between normal and BC patients ####
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
  theme(axis.text=element_text(size=16), axis.title = element_text(size=16))
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
# sort by barcode
normVcancerPaired <- normVcancerPaired[order(normVcancerPaired$patient),] 
# summarize sample counts
table(normVcancerPaired$shortLetterCode) # 112 NT, 112 TP
# perform t test (paired)
t.test(SPANXB1 ~ shortLetterCode, paired=TRUE, data = normVcancerPaired) # 4.177e-09
ggplot(normVcancerPaired, aes(shortLetterCode, SPANXB1)) + 
  ylab("log2 SPANXB1 expression") +
  xlab("tissue type (paired samples)") +
  scale_x_discrete(labels=c("NT" = "normal", "TP" = "tumor")) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=16), axis.title = element_text(size=16))
#ggsave("figures/SPANXB1paired.jpg")

#### Q2 Compare spanxb1 expression between metastatic vs. non metastatic BC patients ####
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
# sort for paired samples
tumVmetPaired <- tumVmetPaired[order(tumVmetPaired$patient)]
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

#### Q3 Compare spanxb1 expression with clinical stages of BC patients ####
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
  xlab(" AJCC stage") +
  scale_x_discrete(labels=c("Stage I" = "I", "Stage IA" = "IA", "Stage IB" = "IB", "Stage II" = "II", "Stage IIA" = "IIA", "Stage IIB" = "IIB", "Stage III" = "III", "Stage IIIA" = "IIIA", "Stage IIIB" = "IIIB", "Stage IIIC" = "IIIC", "Stage IV" = "IV")) +
  geom_boxplot() +
  theme_bw() 
#ggsave("figures/SPANXB1.AJCC.jpg")
# AJCC stages with zeros removed
AJCC_filt <- AJCC %>%
  filter(SPANXB1 > 0)
table(AJCC_filt$subtype_AJCC.Stage)
summary(aov(SPANXB1 ~ subtype_AJCC.Stage, data = AJCC_filt)) # 0.0346
ggplot(AJCC_filt, aes(subtype_AJCC.Stage, SPANXB1)) + 
  ylab("log2 SPANXB1 expression") +
  xlab("AJCC stage") +
  scale_x_discrete(labels=c("Stage I" = "I", "Stage IA" = "IA", "Stage IB" = "IB", "Stage II" = "II", "Stage IIA" = "IIA", "Stage IIB" = "IIB", "Stage III" = "III", "Stage IIIA" = "IIIA", "Stage IIIB" = "IIIB", "Stage IIIC" = "IIIC", "Stage IV" = "IV")) +
  geom_boxplot() +
  theme_bw() 
#ggsave("figures/SPANXB1.AJCC_filt.jpg")
# compare tumor stages
stage <- tum %>%
  filter(tumor_stage != "not reported") %>%
  filter(tumor_stage != "stage x")
table(stage$tumor_stage)
# perform ANOVA
summary(aov(SPANXB1 ~ tumor_stage, data = stage)) #p=0.514
ggplot(stage, aes(tumor_stage, SPANXB1)) + 
  ylab("log2 SPANXB1 expression") +
  xlab("tumor stage") +
  scale_x_discrete(labels=c("stage i" = "I", "stage ia" = "IA", "stage ib" = "IB", "stage ii" = "II", "stage iia" = "IIA", "stage iib" = "IIB", "stage iii" = "III", "stage iiia" = "IIIA", "stage iiib" = "IIIB", "stage iiic" = "IIIC", "stage iv" = "IV")) +
  geom_boxplot() +
  theme_bw() 
#ggsave("figures/SPANXB1.stage.jpg")

#### Q4 Compare spanxb1 expression between normal and ER/PR/HER2 negative (a.k.a.TNBC) patients ####
# extract triple negative from tumor samples
TNBCneg <- tum %>% 
  filter(subtype_ER.Status == "Negative" & subtype_PR.Status == "Negative" & subtype_HER2.Final.Status == "Negative") %>%
  mutate(triple = "TNBC")
# convert norm to include TNBC status
normT <- norm %>%
  mutate(triple = "normal")
# combine triple negative and normal
normTNBC <- rbind(normT, TNBCneg) # 118 TNBC, 113 normal
# perform t test 
t.test(SPANXB1 ~ triple, data=normTNBC) # 4.395e-07
ggplot(normTNBC, aes(triple, SPANXB1)) + 
  ylab("log2 SPANXB1 expression") +
  xlab("tissue type (unpaired)") +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=16), axis.title = element_text(size=16))
#ggsave("figures/SPANXB1TNBCnormal.jpg")

#### Q4.5 Compare spanxb1 expression between ER/PR/HER2 positive vs. ER/PR/HER2 negative (a.k.a.TNBC) patients ####
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

#### Q5 Compare spanxb1 expression between metastatic vs. non metastatic TNBC patients ####
table(TNBCneg$shortLetterCode) #no TNBC neg or pos are metastatic

#### Compare spanxb1 expression with survival outcome ####
tumVital <- tum %>%
  filter(vital_status != "NA")
table(tumVital$vital_status) # 100 alive, 18 dead 
# SPANXB1 and vital status
t.test(SPANXB1 ~ vital_status, data=tumVital) # 0.5902
ggplot(tumVital, aes(vital_status, SPANXB1)) + 
  ylab("log2 SPANXB1 expression") +
  xlab("vital status") +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
#ggsave("figures/SPANXB1.vital.jpg")

# read in curated survival data
surv_dat <- read.csv("brca_survival.csv")
# compare samples in each set
length(tum$patient) #1102
length(unique(tum$patient)) #1091
length(surv_dat$bcr_patient_barcode) #1097
length(unique(surv_dat$bcr_patient_barcode)) #1097
# filter out patients in tum that aren't in survival data
surv <- merge(tum, surv_dat, by.x="patient", by.y="bcr_patient_barcode") %>%
  filter(DFI.time.cr != "NA")
range(surv$DFI.time.cr)

# plot all SPANXB1 and DFI
ggplot(surv, aes(SPANXB1, DFI.time.cr)) +
  geom_point() +
  ylab("disease free interval (days)") +
  xlab("log2 SPANXB1 expression") +
  theme_bw() +
  theme(legend.position="none") + 
  geom_smooth(method = "lm", se = FALSE)
#ggsave("figures/SPANXB1.survival.jpg")
# model with DFI
surv.mod.D <- lm(SPANXB1 ~ DFI.time.cr, data=surv_ex)
summary(surv.mod.D) # p=0.1759, R2=0.003988, no relationship

# plot DFI by SPANXB1 only in expressed
surv_ex <- filter(surv, SPANXB1 > 0)
ggplot(surv_ex, aes(SPANXB1, DFI.time.cr)) +
  geom_point() +
  ylab("disease free interval (days)") +
  xlab("log2 SPANXB1 expression") +
  theme_bw() +
  theme(legend.position="none") + 
  geom_smooth(method = "lm", se = FALSE)
#ggsave("figures/SPANXB1.ex.survival.jpg")
# model with DFI
surv.mod.D.ex <- lm(SPANXB1 ~ DFI.time.cr, data=surv_ex)
summary(surv.mod.D.ex) # p=0.1759, R2=0.003988, no relationship

#### Q6 Compare spanxb1 expression with survival outcome of TNBC patients ####
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

#### Q7 Compare spanxb1 expression with survival outcome of ER/PR/HER2 positive patients alone and in combination (e.g., co-expression of genes) ####
table(TNBCpos$vital_status) # 52 alive, 7 dead 
# SPANXB1 and vital status of triple positive patients
t.test(SPANXB1 ~ vital_status, data=TNBCpos) # 0.9511
ggplot(TNBCpos, aes(vital_status, SPANXB1)) + 
  ylab("log2 SPANXB1 expression") +
  xlab("vital status") +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
#ggsave("figures/SPANXB1.TNBCpos.vital.jpg")
# SPANXB1, ER+ (subtype_ER.Status) and vital status
table(fpkmGene$subtype_ER.Status) # 175 Negative, 591 Positive, Indeterminate, Not Performed, Performed but Not Available
SP.ER <- fpkmGene %>% 
  filter(subtype_ER.Status == "Positive" | subtype_ER.Status == "Negative")
t.test(SPANXB1 ~ subtype_ER.Status, data=SP.ER) # 0.004189
ggplot(SP.ER, aes(subtype_ER.Status, SPANXB1)) + 
  ylab("log2 SPANXB1 expression") +
  xlab("ER status") +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
#ggsave("figures/SPANXB1.ERstatus.jpg")
# logistic regression # 766 individuals
SP.ER.mod <- glm(vital_status ~ SPANXB1 + subtype_ER.Status, data = SP.ER, family = "binomial")
summary(SP.ER.mod) # nope
# odds ratio and 95% CI
exp(cbind(OR = coef(SP.ER.mod), confint(SP.ER.mod)))
# SPANXB1, PR+ (subtype_PR.Status) and vital status
table(fpkmGene$subtype_PR.Status) # 247 Negative, 516 Positive
SP.PR <- fpkmGene %>% 
  filter(subtype_PR.Status == "Positive" | subtype_PR.Status == "Negative")
t.test(SPANXB1 ~ subtype_PR.Status, data=SP.PR) # 0.2112
ggplot(SP.PR, aes(subtype_PR.Status, SPANXB1)) + 
  ylab("log2 SPANXB1 expression") +
  xlab("PR status") +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
#ggsave("figures/SPANXB1.PRstatus.jpg")
# logistic regression # 763 individuals
SP.PR.mod <- glm(vital_status ~ SPANXB1 + subtype_PR.Status, data = SP.PR, family = "binomial")
summary(SP.PR.mod) # nope
# odds ratio and 95% CI
exp(cbind(OR = coef(SP.PR.mod), confint(SP.PR.mod)))
# SPANXB1, HER2+ (subtype_HER2.Final.Status) and vital status
table(fpkmGene$subtype_HER2.Final.Status) # 642 Negative, 111 Positive
SP.HER <- fpkmGene %>% 
  filter(subtype_HER2.Final.Status == "Positive" | subtype_HER2.Final.Status == "Negative")
t.test(SPANXB1 ~ subtype_HER2.Final.Status, data=SP.HER) # 0.3968
ggplot(SP.HER, aes(subtype_HER2.Final.Status, SPANXB1)) + 
  ylab("log2 SPANXB1 expression") +
  xlab("HER2 status") +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
#ggsave("figures/SPANXB1.HERstatus.jpg")
# logistic regression # 753 individuals
SP.HER.mod <- glm(vital_status ~ SPANXB1 + subtype_HER2.Final.Status, data = SP.HER, family = "binomial")
summary(SP.HER.mod) # nope
# odds ratio and 95% CI
exp(cbind(OR = coef(SP.HER.mod), confint(SP.HER.mod)))
# SPANXB1 and ER/PR/HER2+ # 715 individuals
SP.PR.ER.HER <- fpkmGene %>% 
  filter(subtype_PR.Status == "Positive" | subtype_PR.Status == "Negative") %>%
  filter(subtype_ER.Status == "Positive" | subtype_ER.Status == "Negative") %>%
  filter(subtype_HER2.Final.Status == "Positive" | subtype_HER2.Final.Status == "Negative")
# logistic regression
SP.PR.ER.HER.mod <- glm(vital_status ~ SPANXB1 + subtype_PR.Status + subtype_ER.Status + subtype_HER2.Final.Status, data = SP.PR.ER.HER, family = "binomial")
summary(SP.PR.ER.HER.mod) # nope
# odds ratio and 95% CI
exp(cbind(OR = coef(SP.HER.mod), confint(SP.HER.mod)))

#### Q8 Compare RAC1/SPANXB1 expression together in normal vs. TNBC ####
# aggregate data from norm and TNBC
normTNBC <- rbind(norm, TNBCneg[,-95]) # removes column triple
table(normTNBC$shortLetterCode) # 113 normal, 118 TNBC (unpaired)
# linear regression (only normal)
norm.mod <- lm(SPANXB1 ~ RAC1, data=norm)
summary(norm.mod) # p=0.697, R2=0.001371
# linear regression (only TNBC)
TNBCmod <- lm(SPANXB1 ~ RAC1, data=TNBCneg)
summary(TNBCmod) # p=0.1707, R2=0.01612
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
# linear regression (only normal)
normExp <- norm %>%
  filter(SPANXB1 > 0.000000) %>%
  filter(RAC1 > 0.000000)
normExp.mod <- lm(SPANXB1 ~ RAC1, data=normExp)
summary(normExp.mod) # p=0.7759, R2=0.005982
# linear regression (only TNBC)
TNBCexp <- TNBCneg %>%
  filter(SPANXB1 > 0.000000) %>%
  filter(RAC1 > 0.000000)
TNBCexp.mod <- lm(SPANXB1 ~ RAC1, data=TNBCexp)
summary(TNBCexp.mod) # p=0.07648, R2=0.06265
# filter for only expression > 0 for both
bothExp <- normTNBC %>%
  filter(SPANXB1 > 0.000000) %>%
  filter(RAC1 > 0.000000)
# plot both together
ggplot(bothExp, aes(SPANXB1, RAC1, col=shortLetterCode)) +
  geom_point() +
  ylab("log2 RAC1 expression") +
  xlab("log2 SPANXB1 expression") +
  theme_bw() +
  theme(legend.position="none") + # blue=TNBC, red=normal
  geom_smooth(data=subset(bothExp, shortLetterCode == "TP"), method = "lm", se = FALSE) +
  geom_smooth(data=subset(bothExp, shortLetterCode == "NT"), method = "lm", se = FALSE)
#ggsave("figures/SPANXB1.RAC1.TNBC.filtered.jpg")

#### Q9 Compare RAC1/SPANXB1 expression in metastatic vs. not met TNBC ####
table(TNBCneg$shortLetterCode) #no TNBC neg are metastatic

#### Q10 Compare RAC1/SPANXB1 expression with survival of TNBC ####
table(TNBCneg$vital_status) # 100 alive, 18 dead
# linear regression (only live)
TNBClive <- TNBCneg %>%
  filter(vital_status == "alive")
TNBClive.mod <- lm(SPANXB1 ~ RAC1, data=TNBClive)
summary(TNBClive.mod) # p=0.2782, R2=0.01199
# linear regression (only dead)
TNBCdead <- TNBCneg %>%
  filter(vital_status == "dead")
TNBCdead.mod <- lm(SPANXB1 ~ RAC1, data=TNBCdead)
summary(TNBCdead.mod) # p=0.4423, R2=0.03735
# plot both together
ggplot(TNBCneg, aes(SPANXB1, RAC1, col=vital_status)) +
  geom_point() +
  ylab("log2 RAC1 expression") +
  xlab("log2 SPANXB1 expression") +
  theme_bw() +
  theme(legend.position="none") + # blue=dead, red=alive
  geom_smooth(data=subset(bothExp, vital_status == "dead"), method = "lm", se = FALSE) +
  geom_smooth(data=subset(bothExp, vital_status == "alive"), method = "lm", se = FALSE)
#ggsave("figures/SPANXB1.RAC1.vital.jpg")
# logistic regression 
SP.RAC.mod <- glm(vital_status ~ SPANXB1 + RAC1, data = TNBCneg, family = "binomial")
summary(SP.RAC.mod) # nope
# odds ratio and 95% CI
exp(cbind(OR = coef(SP.RAC.mod), confint(SP.RAC.mod)))

#### Q11 Compare EGFR/SPANXB1 expression together in normal vs. TNBC ####
table(normTNBC$shortLetterCode) # 113 normal, 118 TNBC (unpaired)
# linear regression (only normal)
norm.mod <- lm(SPANXB1 ~ EGFR, data=norm)
summary(norm.mod) # p=0.7676, R2=0.0007897
# linear regression (only TNBC)
TNBCmod <- lm(SPANXB1 ~ EGFR, data=TNBCneg)
summary(TNBCmod) # p=0.528, R2=0.003441
# plot both together
ggplot(normTNBC, aes(SPANXB1, EGFR, col=shortLetterCode)) +
  geom_point() +
  ylab("log2 EGFR expression") +
  xlab("log2 SPANXB1 expression") +
  theme_bw() +
  theme(legend.position="none") + # blue=TNBC, red=normal
  geom_smooth(data=subset(normTNBC, shortLetterCode == "TP"), method = "lm", se = FALSE) +
  geom_smooth(data=subset(normTNBC, shortLetterCode == "NT"), method = "lm", se = FALSE)
#ggsave("figures/SPANXB1.EGFR.jpg")

#### Q12 Compare EGFR/SPANXB1 expression with survival of TNBC ####
# logistic regression 
SP.EG.mod <- glm(vital_status ~ SPANXB1 + EGFR, data = TNBCneg, family = "binomial")
summary(SP.EG.mod) # nope
# odds ratio and 95% CI
exp(cbind(OR = coef(SP.EG.mod), confint(SP.EG.mod)))

#### Compare SH3GL2/SPANXB1 expression together in normal vs. TNBC ####
hist(norm$SPANXB1)
hist(norm$SH3GL2)
# coexpression in normal (113 samples)
norm.mod <- lm(SPANXB1 ~ SH3GL2, data=norm)
summary(norm.mod) # p=0.5118, R2=0.003887, not coexpressed in normal
# coexpression in tumor (1102 samples)
tum.mod <- lm(SPANXB1 ~ SH3GL2, data=tum)
summary(tum.mod) # p=0.005938, R2=0.00686, coexpressed in tumor
# coexpression in TNBC (118 samples)
TNBC.mod <- lm(SPANXB1 ~ SH3GL2, data=TNBCneg)
summary(TNBC.mod) # p=0.004334, R2=0.06803, coexpressed in TNBC
# plot coexpression in normal vs tumor
ggplot(normVcancer, aes(SPANXB1, SH3GL2, col=shortLetterCode)) +
  geom_point() +
  ylab("log2 SH3GL2 expression") +
  xlab("log2 SPANXB1 expression") +
  theme_bw() +
  theme(legend.position="none") + # blue=tumor, red=normal
  geom_smooth(data=subset(normTNBC, shortLetterCode == "TP"), method = "lm", se = FALSE) +
  geom_smooth(data=subset(normTNBC, shortLetterCode == "NT"), method = "lm", se = FALSE)
#ggsave("figures/SPANXB1.SH3GL2.tumor.jpg")
# plot coexpression in normal vs TNBC
ggplot(normTNBC, aes(SPANXB1, SH3GL2, col=triple)) +
  geom_point() +
  ylab("log2 SH3GL2 expression") +
  xlab("log2 SPANXB1 expression") +
  theme_bw() +
  theme(legend.position="none") + # blue=TNBC, red=normal
  geom_smooth(data=subset(normTNBC, triple == "TNBC"), method = "lm", se = FALSE) +
  geom_smooth(data=subset(normTNBC, triple == "normal"), method = "lm", se = FALSE)
#ggsave("figures/SPANXB1.SH3GL2.TNBC.jpg")
# plot coexpression in tumor only
tum_above <- tum %>%
  filter(SH3GL2 > 0 & SPANXB1 > 0)
ggplot(tum_above, aes(SPANXB1, SH3GL2)) +
  geom_point() +
  ylab("log2 SH3GL2 expression") +
  xlab("log2 SPANXB1 expression") +
  theme_bw() +
  geom_smooth(method = "lm", se = FALSE)
# plot coexpression in TNBC only 
TNBC_above <- TNBCneg %>%
  filter(SH3GL2 > 0 & SPANXB1 > 0)
ggplot(TNBC_above, aes(SPANXB1, SH3GL2)) +
  geom_point() +
  ylab("log2 SH3GL2 expression") +
  xlab("log2 SPANXB1 expression") +
  theme_bw() +
  geom_smooth(method = "lm", se = FALSE)
# expression ratios: normal v tumor
t.test(SH3GL2_SPANXB1_ratio ~ shortLetterCode, data = normVcancer) 
ggplot(normVcancer, aes(SH3GL2_SPANXB1_ratio, fill = shortLetterCode)) +
  geom_density(alpha = 0.5) +
  xlab("SH3GL2:SPANXB1 expression ratio") + # red = normal, blue = tumor
  scale_fill_discrete(guide=FALSE) +
  theme_bw()
#ggsave("figures/expressionRatioNormVtum.jpg")
# expression ratios: TNBC
ggplot(TNBCboth, aes(SH3GL2_SPANXB1_ratio, fill = triple)) +
  geom_density(alpha = 0.5) +
  xlab("SH3GL2:SPANXB1 expression ratio") + # red = TNBC, blue = all positive
  scale_fill_discrete(guide=FALSE) +
  theme_bw()
#ggsave("figures/expressionRatioTNBC.jpg")

#### Compare SH3GL2/SPANXB1 expression together in primary tumor vs. metastatic ####
# coexpression in metastatic (7 samples)
meta.mod <- lm(SPANXB1 ~ SH3GL2, data=meta)
summary(meta.mod) # p=0.2932, R2=0.2162, not coexpressed in metastatic
ggplot(tumVmetaAll, aes(SPANXB1, SH3GL2, col=shortLetterCode)) + 
  geom_point() +
  ylab("log2 SH3GL2 expression") +
  xlab("log2 SPANXB1 expression") +
  theme_bw() +
  theme(legend.position="none") + # blue=primary tumor, red=metastatic
  geom_smooth(data=subset(normTNBC, shortLetterCode == "TP"), method = "lm", se = FALSE) +
  geom_smooth(data=subset(normTNBC, shortLetterCode == "TM"), method = "lm", se = FALSE)
#ggsave("figures/SPANXB1unpairedMetastasis.jpg")

#### Compare SH3GL2/SPANXB1 expression with survival ####
TNBCboth$triple <- as.factor(TNBCboth$triple)
# triple status given SPANXB1 and SH3GL2
SP.SH.TNBC.mod <- glm(triple ~ SPANXB1 + SH3GL2, data=TNBCboth, family="binomial")
summary(SP.SH.TNBC.mod) # nope
# vital status given SPANXB1 and SH3GL2
SP.SH.vital.mod <- glm(vital_status ~ SPANXB1 + SH3GL2, data = normTNBC, family = "binomial")
summary(SP.SH.vital.mod) # nope
# plot both genes, color by survival
ggplot(normTNBC, aes(SPANXB1, SH3GL2, col=shortLetterCode)) +
  geom_point() +
  ylab("log2 SH3GL2 expression") +
  xlab("log2 SPANXB1 expression") +
  theme_bw() +
  theme(legend.position="none") + # blue=TNBC, red=normal
  geom_smooth(data=subset(normTNBC, shortLetterCode == "TP"), method = "lm", se = FALSE) +
  geom_smooth(data=subset(normTNBC, shortLetterCode == "NT"), method = "lm", se = FALSE)
#ggsave("figures/SPANXB1.SH3GL2.jpg")
# plot both genes together if overexpressed
SP.SH.pos <- rbind(norm, tum) %>%
  filter(SPANXB1 > 1 & SH3GL2 > 1) #488 samples
table(SP.SH.pos$shortLetterCode) # 16 normal, 472 tumor
ggplot(SP.SH.pos, aes(SPANXB1, SH3GL2, col=shortLetterCode)) +
  geom_point() +
  ylab("log2 SH3GL2 expression") +
  xlab("log2 SPANXB1 expression") +
  theme_bw() +
  theme(legend.position="none") + # blue=TNBC, red=normal
  geom_smooth(data=subset(SP.SH.pos, shortLetterCode == "TP"), method = "lm", se = FALSE) +
  geom_smooth(data=subset(SP.SH.pos, shortLetterCode == "NT"), method = "lm", se = FALSE)
#ggsave("figures/SPANXB1.SH3GL2.positive.jpg")
# vital status given SPANXB1 in all TNBC
t.test(SPANXB1 ~ vital_status, data = TNBCneg) # p=0.1257
wilcox.test(SPANXB1 ~ vital_status, data = TNBCneg) # p=0.1225
# vital status given SPANXB1 in only SPANXB1+ tumors
SP.TNBC <- TNBCneg %>%
  filter(SPANXB1 > 1)
t.test(SPANXB1 ~ vital_status, data = SP.TNBC) # p=0.9059
wilcox.test(SPANXB1 ~ vital_status, data = SP.TNBC) # p=0.955
# vital status given SPANXB1 and SH3GL2 in all TNBC
SP.SH.mod <- glm(vital_status ~ SPANXB1 + SH3GL2, data = TNBCneg, family = "binomial")
summary(SP.SH.mod) # nope
# odds ratio and 95% CI
exp(cbind(OR = coef(SP.SH.mod), confint(SP.SH.mod)))
# vital status given SPANXB1 and SH3GL2 in only SPANXB1+ TNBC tumors
SP.SH.TNBC.mod <- glm(vital_status ~ SPANXB1 + SH3GL2, data = SP.TNBC, family = "binomial")
summary(SP.SH.TNBC.mod) # nope
# odds ratio and 95% CI
exp(cbind(OR = coef(SP.SH.mod), confint(SP.SH.mod)))
# vital status of SPANXB1+/- and SH3GL2+/- in only TNBC
SP.SH.table <- table(SP.TNBC$SH3GL2 >1, SP.TNBC$vital_status)
chisq.test(SP.SH.table, simulate.p.value = TRUE) # 0.5667
# vital status of SPANXB1+/- and SH3GL2+/- in all tumors
SP.SH.pos <- tum %>% 
  filter(SPANXB1 > 0 & SH3GL2 > 0) %>%
  mutate(SP.SH = "positive") # 472
SP.SH.neg <- tum %>% 
  filter(SPANXB1 == 0 | SH3GL2 == 0) %>%
  mutate(SP.SH = "negative")
SP.SH <- rbind(SP.SH.pos, SP.SH.neg)
SP.SH.tbl <- table(SP.SH$SP.SH, SP.SH$vital_status)
chisq.test(SP.SH.tbl, simulate.p.value = TRUE) # p=0.5467
# vital status of SPANXB1+/- and SH3GL2+/- in all tumors, strict comparison
SP.SH.part <- tum %>%
  filter(SPANXB1 ==0 & SH3GL2 ==0) %>%
  mutate(SP.SH="negative")
SP.SH2 <- rbind(SP.SH.pos, SP.SH.part)
SP.SH.tbl2 <- table(SP.SH$SP.SH, SP.SH$vital_status)
chisq.test(SP.SH.tbl2, simulate.p.value = TRUE) # p=0.5617

## use Wilcoxin signed-rank or Mann-Whitney for small sample sizes?

## web tool for Kaplan Meier plots: http://kmplot.com/analysis/index.php?p=service&cancer=breast#

## correct for multiple comparisons?

# rlog from DESeq2 as alternative transformation method

## permutation tests (requires at least two groups of six samples)
# calculate test statistic for dataset
# permute labels on samples at random 1000 times
# recalculate test statistic for each of 1000 permuted datasets
# calculate percentage of permuted test statistics exceed original

