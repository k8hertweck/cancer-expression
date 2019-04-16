# example kaplan-meier code
# https://genomicsclass.github.io/book/pages/tcga.html
# https://github.com/kassambara/survminer
# https://www.biostars.org/p/267016/
# https://www.biostars.org/p/300304/#300420 (cox models)

library(survival)
library(survminer)
library(TCGAbiolinks)
library(dplyr)

#### load example clinical data ####
clin_all <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
# remove missing data for days to last followup
clin <- clin_all %>%
  filter(!is.na(days_to_last_follow_up)) %>%
  filter(ethnicity != "not reported")
clin$days_to_last_follow_up <- as.numeric(clin$days_to_last_follow_up)
clin <- clin %>% 
  mutate(vital = (as.numeric(as.factor(vital_status))) - 1)

#### survival plot ####

# important functions
?Surv # Create a survival object, usually used as a response variable in a model formula.
?survfit # creates survival curves from a formula (e.g. the Kaplan-Meier)
?ggsurvplot

# survival plot
# create event, or status indicator: 0=alive, 1=dead
# create survival object
su <- Surv(clin$days_to_last_follow_up, clin$vital)
# plot survival curves 
plot(survfit(su~ethnicity, data=clin), lwd=2, lty=1:4)
ntab <- table(clin$ethnicity)
ns <- paste("[n=", ntab, "]", sep="")
legend(100, .6, lty=1:4, lwd=2, legend=paste(levels(clin$ethnicity), ns))

# fit model
fit <- survfit(Surv(days_to_last_follow_up, 
                    vital) ~ ethnicity, 
               data = clin)
# Visualize with survminer
ggsurvplot(fit, data = clin, risk.table = TRUE)
