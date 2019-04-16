#### head and neck expression analysis ####

# load packages
library(dplyr)
library(ggplot2)

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

# high expression in AA is correlated with poor outcome when compared to non-hispanic whites
# vital status
t.test(MDA9 ~ vital_status, data=BW_nonhis) # 0.004855
ggplot(BW_nonhis, aes(vital_status, MDA9)) + 
  ylab("log2 MDA9 expression") +
  xlab("vital status") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.3) +
  theme_bw() 
#ggsave("figures/MDA9BW_nonhis_vital_.jpg")
# kaplan meier
https://genomicsclass.github.io/book/pages/tcga.html

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

#### SIRPA ####

# high expression in AA is correlated with poor outcome when compared to non-hispanic whites
t.test(SIRPA ~ vital_status, data=BW_nonhis) # 0.3977
ggplot(BW_nonhis, aes(vital_status, SIRPA)) + 
  ylab("log2 MDA9 expression") +
  xlab("vital status") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.3) +
  theme_bw() 