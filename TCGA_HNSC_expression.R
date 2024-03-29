#### head and neck expression analysis ####

# load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, survminer, survival, viridis)

#### data pre-processing ####

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

# assessing clinical data
table(fpkmGene$shortLetterCode) # 44 NT, 2 TM, 500 TP
table(fpkmGene$race) # 49 AA, 468 white
table(fpkmGene$ethnicity) # 27 HIS, 480 not (some overlap with AA?)
table(fpkmGene$classification_of_tumor) # not reported
table(fpkmGene$tumor_stage) # not reported, stage i - ivc
table(fpkmGene$tumor_grade) # not reported
table(fpkmGene$progression_or_recurrence) # not reported
table(fpkmGene$shortLetterCode) # only 2 metastasis

# extract tumor data
tum <- fpkmGene %>%
  filter(shortLetterCode == "TP")
# extracting AA and white who are not HIS
BW_nonhis <- tum %>%
  filter(race == "white" | race == "black or african american") %>%
  filter(ethnicity == "not hispanic or latino")
# extracting whites with known HIS
his_nonhis <- tum %>%
  filter(race == "white") %>%
  filter(ethnicity == "hispanic or latino" | ethnicity == "not hispanic or latino") 

#### MDA-9/Syntenin ####

## high expression in AA is correlated with poor outcome when compared to non-hispanic whites
## compare vital status based on MDA9
t.test(MDA9 ~ vital_status, data=BW_nonhis) # 0.001108
ggplot(BW_nonhis, aes(vital_status, MDA9)) + 
  ylab("log2 MDA9 expression") +
  xlab("vital status") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.3) +
  theme_bw() 

## MDA-9 expression in AA/CA by tumor stage
BW_nonhis_levels <- BW_nonhis %>%
  filter(tumor_stage != "not reported") %>%
  droplevels()
race_summary <- BW_nonhis_levels %>%  
  group_by(race, tumor_stage, .drop=FALSE) %>%
  tally()
ggplot(BW_nonhis_levels, aes(x = tumor_stage, y = MDA9, color=race)) +
  ylab("log2 MDA9 expression") +
  xlab("tumor stage") +
  scale_x_discrete(labels=c("stage i" = "I", "stage ii" = "II", "stage iii" = "III", "stage iva" = "IVA", "stage ivb" = "IVB", "stage ivc" = "IVC")) +
  geom_boxplot() +
  geom_text(data=filter(race_summary, race=="black or african american"), 
            aes(tumor_stage, Inf, label = n), vjust=2, hjust=2) +
  geom_text(data=filter(race_summary, race=="white"), 
            aes(tumor_stage, Inf, label = n), vjust=2, hjust=-0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.title.y = element_text(size=12))
#ggsave("figures/AA_CA_tumor_stage.jpeg", width = 8)
table(BW_nonhis$tumor_stage, BW_nonhis$race)

## MDA-9 expression in CA/HIS by tumor stage
his_nonhis_levels <- his_nonhis %>%
  filter(tumor_stage != "not reported", tumor_stage != "stage ivc") %>%
  droplevels()
ethnicity_summary <- his_nonhis_levels %>% 
  group_by(ethnicity, tumor_stage, .drop=FALSE) %>%
  tally()
ggplot(his_nonhis_levels, aes(tumor_stage, MDA9, color=ethnicity)) +
  ylab("log2 MDA9 expression") +
  xlab("tumor stage") +
  scale_x_discrete(labels=c("stage i" = "I", "stage ii" = "II", "stage iii" = "III", "stage iva" = "IVA", "stage ivb" = "IVB")) +
  geom_boxplot() +
  geom_text(data=filter(ethnicity_summary, ethnicity=="hispanic or latino"), 
            aes(tumor_stage, Inf, label = n), vjust=2, hjust=2) +
  geom_text(data=filter(ethnicity_summary, ethnicity=="not hispanic or latino"), 
            aes(tumor_stage, Inf, label = n), vjust=2, hjust=-0.5) +
  theme_bw()  +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.title.y = element_text(size=12))  
#ggsave("figures/his_nonhis_tumor_stage.jpeg", width = 7)
table(his_nonhis$tumor_stage, his_nonhis$ethnicity)

# KM from TCGA biolinks
#TCGAbiolinks::TCGAanalyze_SurvivalKM()

## kaplan meier: days to last followup, comparing AA and CA
# remove missing data 
BW_nonhis_follow <- BW_nonhis %>%
  filter(!is.na(days_to_last_follow_up)) 
BW_nonhis_follow <- BW_nonhis_follow %>% 
  mutate(vital = (as.numeric(vital_status)) - 1)
# fit model
BW_nonhis_follow_fit <- survfit(Surv(days_to_last_follow_up, 
                    vital) ~ race, 
               data = BW_nonhis_follow)
ggsurvplot(BW_nonhis_follow_fit, data = BW_nonhis_follow, 
           risk.table = TRUE,
           pval = TRUE,
           pval.method = TRUE,
           pval.coord = c(1600, 0.8),
           pval.method.coord = c(1600, 0.9),
           conf.int = TRUE,
           xlim = c(0, 2000),
           break.time.by = 500,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = T,
           risk.table.y.text = FALSE)

## kaplan meier: days to death, comparing AA and CA (better than days to last followup)
# remove missing data 
BW_nonhis_death <- BW_nonhis %>%
  filter(!is.na(days_to_death)) 
BW_nonhis_death <- BW_nonhis_death %>% 
  mutate(vital = (as.numeric(vital_status)) - 1)
# fit model
BW_nonhis_death_fit <- survfit(Surv(days_to_death, 
                         vital) ~ race, 
                    data = BW_nonhis_death)
ggsurvplot(BW_nonhis_death_fit, 
           data = BW_nonhis_death, 
           risk.table = TRUE,
           pval = TRUE,
           pval.method = TRUE,
           pval.coord = c(2000, 0.8),
           pval.method.coord = c(2000, 0.9),
           conf.int = TRUE,
           xlim = c(0, 3000),
           break.time.by = 500,
           xlab = "days to death",
           ggtheme = theme_bw(),
           risk.table.y.text.col = T,
           risk.table.y.text = FALSE,
           legend.labs = c("AA", "CA"))
#ggsave("figures/AA_CA_days_to_death.jpeg") # width 700, height 525
# splitting into high and low gene expression
BW_nonhis_death <- BW_nonhis_death %>%
  mutate(mda9_hl = MDA9 > median(MDA9))
BW_nonhis_death$mda9_hl[BW_nonhis_death$mda9_hl == TRUE] <- "high"
BW_nonhis_death$mda9_hl[BW_nonhis_death$mda9_hl == FALSE] <- "low"
# high and low gene expression
BW_nonhis_death_mda9_fit <- survfit(Surv(days_to_death, vital) ~ mda9_hl, 
                                         data = BW_nonhis_death)
ggsurvplot(BW_nonhis_death_mda9_fit, 
           data = BW_nonhis_death, 
           risk.table = TRUE,
           pval = TRUE,
           pval.method = TRUE,
           pval.coord = c(2000, 0.8),
           pval.method.coord = c(2000, 0.9),
           conf.int = TRUE,
           xlim = c(0, 3000),
           break.time.by = 500,
           xlab = "days to death",
           ggtheme = theme_bw(),
           risk.table.y.text.col = TRUE,
           risk.table.y.text = FALSE,
           legend.labs = c("high", "low")
           )
#ggsave("figures/days_to_death_MDA9_AA_CA.jpg")
# AA vs CA, high and low gene expression
BW_nonhis_death_mda9_fit <- survfit(Surv(days_to_death, vital) ~ mda9_hl + race, 
                                    data = BW_nonhis_death)
ggsurvplot(BW_nonhis_death_mda9_fit, 
           data = BW_nonhis_death, 
           palette = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c"),
           risk.table = TRUE,
           pval = TRUE,
           pval.method = TRUE,
           pval.coord = c(2000, 0.8),
           pval.method.coord = c(2000, 0.9),
           conf.int = TRUE,
           xlim = c(0, 3000),
           break.time.by = 500,
           xlab = "days to death",
           ggtheme = theme_minimal(),
           risk.table.y.text.col = TRUE,
           risk.table.y.text = FALSE, 
           legend.labs = c("AA, high", "CA, high", "AA, low", "CA, low"))
#ggsave("figures/AA_CA_days_to_death_MDA9.jpg")

# high expression in HIS is correlated with poor outcome when compared to non-hispanic whites
# vital status
t.test(MDA9 ~ vital_status, data=his_nonhis) # 0.0008478
ggplot(his_nonhis, aes(vital_status, MDA9)) + 
  ylab("log2 MDA9 expression") +
  xlab("vital status") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.3) +
  theme_bw() 

## kaplan meier: days to last followup, comparing his and nonhis
# remove missing data 
his_nonhis_follow <- his_nonhis %>%
  filter(!is.na(days_to_last_follow_up)) 
his_nonhis_follow <- his_nonhis_follow %>% 
  mutate(vital = (as.numeric(vital_status)) - 1)
# fit model
his_nonhis_follow_fit <- survfit(Surv(days_to_last_follow_up, 
                                     vital) ~ ethnicity, 
                                data = his_nonhis_follow)
ggsurvplot(his_nonhis_follow_fit, data = his_nonhis_follow, 
           risk.table = TRUE,
           pval = TRUE,
           pval.method = TRUE,
           pval.coord = c(1600, 0.8),
           pval.method.coord = c(1600, 0.9),
           conf.int = TRUE,
           xlim = c(0, 2000),
           break.time.by = 500,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = T,
           risk.table.y.text = FALSE)

## kaplan meier: days to death, comparing his and nonhis (better than days to last followup)
# remove missing data 
his_nonhis_death <- his_nonhis %>%
  filter(!is.na(days_to_death)) 
his_nonhis_death <- his_nonhis_death %>% 
  mutate(vital = (as.numeric(vital_status)) - 1)
# fit model
his_nonhis_death_fit <- survfit(Surv(days_to_death, 
                                    vital) ~ ethnicity, 
                               data = his_nonhis_death)
ggsurvplot(his_nonhis_death_fit, 
           data = his_nonhis_death, 
           risk.table = TRUE,
           pval = TRUE,
           pval.method = TRUE,
           pval.coord = c(2000, 0.8),
           pval.method.coord = c(2000, 0.9),
           conf.int = TRUE,
           xlim = c(0, 3000),
           break.time.by = 500,
           xlab = "days to death",
           ggtheme = theme_bw(),
           risk.table.y.text.col = T,
           risk.table.y.text = FALSE,
           legend.labs = c("hispanic", "not hispanic"))
#ggsave("figures/his_nonhis_days_to_death.jpeg") # width 700, height 525
# splitting into high and low gene expression
his_nonhis_death <- his_nonhis_death %>%
  mutate(mda9_hl = MDA9 > median(MDA9))
his_nonhis_death$mda9_hl[his_nonhis_death$mda9_hl == TRUE] <- "high"
his_nonhis_death$mda9_hl[his_nonhis_death$mda9_hl == FALSE] <- "low"
# high and low gene expression
his_nonhis_death_mda9_fit <- survfit(Surv(days_to_death, vital) ~ mda9_hl, 
                                    data = his_nonhis_death)
ggsurvplot(his_nonhis_death_mda9_fit, 
           data = his_nonhis_death, 
           risk.table = TRUE,
           pval = TRUE,
           pval.method = TRUE,
           pval.coord = c(2000, 0.8),
           pval.method.coord = c(2000, 0.9),
           conf.int = TRUE,
           xlim = c(0, 3000),
           break.time.by = 500,
           xlab = "days to death",
           ggtheme = theme_bw(),
           risk.table.y.text.col = TRUE,
           risk.table.y.text = FALSE,
           legend.labs = c("high", "low")
)
#ggsave("figures/days_to_death_MDA9_hisnonhis.jpg")
# his vs non, high and low gene expression
his_nonhis_death_mda9_fit <- survfit(Surv(days_to_death, vital) ~ mda9_hl + ethnicity, 
                                    data = his_nonhis_death)
ggsurvplot(his_nonhis_death_mda9_fit, 
           data = his_nonhis_death, 
           palette = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c"),
           risk.table = TRUE,
           xlim = c(0, 3000),
           break.time.by = 500,
           xlab = "days to death",
           ggtheme = theme_minimal(),
           risk.table.y.text.col = TRUE,
           risk.table.y.text = FALSE, 
           legend.labs = c("hispanic, high", "not hispanic, high", "hispanic, low", "not hispanic, low"))
#ggsave("figures/his_nonhis_days_to_death_MDA9.jpg")

#### SIRPA: not needed right now ####

# high expression in AA is correlated with poor outcome when compared to non-hispanic whites
t.test(SIRPA ~ vital_status, data=BW_nonhis) # 0.3977
ggplot(BW_nonhis, aes(vital_status, SIRPA)) + 
  ylab("log2 MDA9 expression") +
  xlab("vital status") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.3) +
  theme_bw() 
