# analyzing expression data from TCGA 

# load packages
library(dplyr)
library(ggplot2)

# read in saved data
fpkmGene <- read.table("targetGeneBca.csv")
# force race to factor
fpkmGene$race <- as.factor(fpkmGene$race)

## visualizing data distribution
hist(fpkmGene$SPANXB1)
plot(fpkmGene$race)
plot(fpkmGene$morphology)
plot(fpkmGene$shortLetterCode)

## statistical testing

# does gene expression differ between normal and prostate? (unpaired)
normVcancer <- filter(fpkmGene, shortLetterCode != "TM") # remove metastasis
# unpaired, all data
t.test(SPANXB1 ~ shortLetterCode, data = normVcancer) # p=0.008919
SB1 <- normVcancer %>%
  filter(SPANXB1 > 0) # 541 samples with SPANXB1 expression
table(normVcancer$shortLetterCode) # 113 normal, 1102 tumor
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
# create paired sample dataset
normID <- normVcancer %>% # list normal samples
  filter(shortLetterCode == "NT") %>%
  select(bcr_patient_barcode)
normVcancerPaired <- normVcancer %>%
  filter(bcr_patient_barcode %in%
           normID$bcr_patient_barcode)
# find patients with more than 1 tumor sample to pair
normVcancerPaired %>% 
  group_by(patient) %>%
  tally() %>%
  filter(n > 2)
# find patients with more than 1 sample of any type
normVcancerPaired %>%
  group_by(patient) %>%
  tally() %>%
  filter(n == 1)
# only one sample: TCGA-BH-A0BS (coded NT/normal, but also stage iiia?)
filter(normVcancerPaired, patient == "TCGA-BH-A0BS")
grep("TCGA-BH-A0BS", normVcancerPaired$patient) # remove 145 
# three or four samples: TCGA-A7-A0DB, TCGA-A7-A0DC, TCGA-A7-A13E
filter(normVcancerPaired, patient == "TCGA-A7-A0DB")
grep("TCGA-A7-A0DB", normVcancerPaired$patient) # remove 123 and 201 (extra tumor)
filter(normVcancerPaired, patient == "TCGA-A7-A0DC")
grep("TCGA-A7-A0DC", normVcancerPaired$patient) # remove 60 (tumor missing metadata)
filter(normVcancerPaired, patient == "TCGA-A7-A13E")
grep("TCGA-A7-A13E", normVcancerPaired$patient) # remove 213 and 220 (extra tumor)

# remove extra samples
normVcancerPaired <- normVcancerPaired[-c(145, 123, 201, 60, 213, 220),]
# summarize sample counts
table(normVcancerPaired$shortLetterCode) # 112 NT, 112 TP
# SPANXB1 paired
t.test(SPANXB1 ~ shortLetterCode, data = normVcancerPaired, paired=TRUE) # p=0.0007278
t.test(SPANXB1 ~ shortLetterCode, data = normVcancerPaired) # 0.0006777
ggplot(normVcancerPaired, aes(shortLetterCode, SPANXB1)) + 
  ylab("SPANXB1 expression") +
  xlab("tissue type (paired samples)") +
  scale_x_discrete(labels=c("NT" = "normal", "TP" = "tumor")) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
ggsave("figures/SPANXB1paired.jpg")

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
ggsave("figures/SPANXB1vital.jpg")

# add other SPANX?