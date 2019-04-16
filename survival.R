# example kaplan-meier code
# https://genomicsclass.github.io/book/pages/tcga.html

library(survival)
library(survminer)
library(TCGAbiolinks)
library(dplyr)

# example clinical data
clin_all <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
# remove missing data for days to last followup
clin <- clin_all %>%
  filter(!is.na(days_to_last_follow_up))

# survival plot
ev <- 1*(clin$vital_status == 1)
fut <- as.numeric(clin$days_to_last_follow_up)
su <- Surv(fut, ev)
plot(survfit(su~t_stage, data=clin), lwd=2, lty=1:4, xlim=c(0,2000))
ntab <- table(clin$t_stage)
ns <- paste("[n=", ntab, "]", sep="")
legend(100, .4, lty=1:4, lwd=2, legend=paste(levels(clin$), ns))

# fit model
fit <- survfit(Surv(times, patient.vital_status) ~ admin.disease_code,
               data = BRCA.survInfo)
?survfit
# Visualize with survminer
ggsurvplot(fit, data = BRCA.survInfo, risk.table = TRUE)
?ggsurvplot
