#### Functions for TCGA data ####

# define exclusion subsetting
`%ni%` <- Negate(`%in%`) 

#### tumor: extract only tumor data ####
tumor <- function(dat){
  dat %>%
    filter(shortLetterCode == "TP")
}

#### AAvsCA: extracting AA and CA who are not HIS for tumor samples ####
AAvsCA <- function(dat){
  dat %>%
    tumor() %>%
    filter(race == "white" | race == "black or african american") %>%
    filter(ethnicity == "not hispanic or latino")
} 
# test
AAvsCA(tum)

#### HISvsNonHIS: HIS and nonHIS from CA for tumor samples ####
HISvsNonHIS <- function(dat){
  dat %>%
    tumor() %>%
    filter(race == "white") %>%
    filter(ethnicity == "hispanic or latino" | ethnicity == "not hispanic or latino")
}
# test
HISvsNonHIS(tum)
