# survival plots

install.packages('survminer')
source("https://bioconductor.org/biocLite.R")
biocLite("RTCGA.clinical")

library(survminer)
library(RTCGA.clinical)
library(survival)

# load data
survivalTCGA(BRCA.clinical,
             extract.cols = "admin.disease_code") -> BRCA.survInfo

# fit model
fit <- survfit(Surv(times, patient.vital_status) ~ admin.disease_code,
               data = BRCA.survInfo)

# Visualize with survminer
ggsurvplot(fit, data = BRCA.survInfo, risk.table = TRUE)
