#### head and neck expression analysis ####

# load packages
library(dplyr)
library(ggplot2)
library(survminer)
library(survival)

# define exclusion subsetting
`%ni%` <- Negate(`%in%`) 

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

# extract tumor data
tum <- fpkmGene %>%
  filter(shortLetterCode == "TP")
# extracting whites with known HIS
his_nonhis <- fpkmGene %>%
  filter(race == "white") %>%
  filter(ethnicity == "hispanic or latino" | ethnicity == "not hispanic or latino") 
# extracting AA and white who are not HIS
BW_nonhis <- fpkmGene %>%
  filter(race == "white" | race == "black or african american") %>%
  filter(ethnicity == "not hispanic or latino") 

#### MDA-9/Syntenin ####

## high expression in AA is correlated with poor outcome when compared to non-hispanic whites
## compare vital status based on MDA9
t.test(MDA9 ~ vital_status, data=BW_nonhis) # 0.004855
ggplot(BW_nonhis, aes(vital_status, MDA9)) + 
  ylab("log2 MDA9 expression") +
  xlab("vital status") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.3) +
  theme_bw() 
#ggsave("figures/MDA9BW_nonhis_vital_.jpg")

# KM from TCGA biolinks
#TCGAbiolinks::TCGAanalyze_SurvivalKM()

## kaplan meier: days to last followup, comparing AA and CA
# remove missing data 
BW_nonhis_km <- BW_nonhis %>%
  filter(!is.na(days_to_last_follow_up)) 
BW_nonhis_km <- BW_nonhis_km %>% 
  mutate(vital = (as.numeric(vital_status)) - 1)
# fit model
mda9_fit_BW_nonhis_follow <- survfit(Surv(days_to_last_follow_up, 
                    vital) ~ race, 
               data = BW_nonhis_km)
ggsurvplot(mda9_fit_BW_nonhis_follow, data = BW_nonhis_km, 
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

## kaplan meier: days to death, comparing AA and CA
# remove missing data 
BW_nonhis_km <- BW_nonhis %>%
  filter(!is.na(days_to_death)) 
BW_nonhis_km <- BW_nonhis_km %>% 
  mutate(vital = (as.numeric(vital_status)) - 1)
# fit model
mda9_fit_BW_nonhis_death <- survfit(Surv(days_to_death, 
                         vital) ~ race, 
                    data = BW_nonhis_km)
ggsurvplot(mda9_fit_BW_nonhis_death, 
           data = BW_nonhis_km, 
           risk.table = TRUE,
           pval = TRUE,
           pval.method = TRUE,
           pval.coord = c(2000, 0.8),
           pval.method.coord = c(2000, 0.9),
           conf.int = TRUE,
           xlim = c(0, 3000),
           break.time.by = 500,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = T,
           risk.table.y.text = FALSE)
# splitting into high and low gene expression
BW_nonhis_km <- BW_nonhis_km %>%
  mutate(mda9_hl = MDA9 > median(MDA9))
BW_nonhis_km$mda9_hl[BW_nonhis_km$mda9_hl == TRUE] <- "high"
BW_nonhis_km$mda9_hl[BW_nonhis_km$mda9_hl == FALSE] <- "low"
# AA and CA, high and low
mda9_fit_BW_nonhis_death_gene <- survfit(Surv(days_to_death, vital) ~ mda9_hl, 
                                    data = BW_nonhis_km)
ggsurvplot(mda9_fit_BW_nonhis_death_gene, 
           data = BW_nonhis_km, 
           risk.table = TRUE,
           pval = TRUE,
           pval.method = TRUE,
           pval.coord = c(2000, 0.8),
           pval.method.coord = c(2000, 0.9),
           conf.int = TRUE,
           xlim = c(0, 3000),
           break.time.by = 500,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = T,
           risk.table.y.text = FALSE)
# AA, high and low (median from all gene expression)
AA_nonhis_km <- filter(BW_nonhis_km, race == "black or african american")
mda9_fit_AA_nonhis_death_gene <- survfit(Surv(days_to_death, vital) ~ mda9_hl, 
                                         data = AA_nonhis_km)
ggsurvplot(mda9_fit_AA_nonhis_death_gene, 
           data = AA_nonhis_km, 
           risk.table = TRUE,
           pval = TRUE,
           pval.method = TRUE,
           pval.coord = c(2000, 0.8),
           pval.method.coord = c(2000, 0.9),
           conf.int = TRUE,
           xlim = c(0, 3000),
           break.time.by = 500,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = T,
           risk.table.y.text = FALSE)


# high expression in HIS is correlated with poor outcome when compared to non-hispanic whites
# vital status
t.test(MDA9 ~ vital_status, data=his_nonhis) # 0.004231
ggplot(his_nonhis, aes(vital_status, MDA9)) + 
  ylab("log2 MDA9 expression") +
  xlab("vital status") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.3) +
  theme_bw() 
#ggsave("figures/MDA9his_nonhis_vital.jpg") 
# days to death

#### SIRPA: not needed right now ####

# high expression in AA is correlated with poor outcome when compared to non-hispanic whites
t.test(SIRPA ~ vital_status, data=BW_nonhis) # 0.3977
ggplot(BW_nonhis, aes(vital_status, SIRPA)) + 
  ylab("log2 MDA9 expression") +
  xlab("vital status") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.3) +
  theme_bw() 
