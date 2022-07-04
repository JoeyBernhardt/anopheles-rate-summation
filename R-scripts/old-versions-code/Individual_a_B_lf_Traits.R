## Kerri Miazgowicz, University of Georgia
## October 1, 2018
## 
## Purpose: Generate indiviudal-specific values for lifespan, daily biting rate, and lifetime egg production for constant and fluctuating temperatures
## ID / temp / a / lf / B <- from master raw csv files
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) Calculate a / lf / B across all datasets
##           3) Export files to be used in TraitFitting.R
##        


#####################
###
#  1) Set-up,load packages, get data, etc.
###
######################
library(dplyr)

#Set working directory
mainDir = "C:/Users/Kerri/Desktop/Fluctuation_BayesianFits"
setwd(mainDir)

#Read in csv files for each dataset
constant.data <- read.csv("data/constant_master.csv")
constant.data[is.na(constant.data)] <- 0
fluc9.data <- read.csv("data/fluc_dtr9_master.csv")
fluc9.data[is.na(fluc9.data)] <- 0
fluc12.data <- read.csv("data/fluc_dtr12_master.csv")
fluc12.data[is.na(fluc12.data)] <- 0


#####################
###
#  2) Calculate a / lf / B for each individual across all datasets
###
######################

#create new dataframe for each temperature profile for the aggregated data

constant.individual.trait <- summarise(group_by(constant.data, Female,Treatment), lifespan = max(Day), bite.rate = sum(Feed)/ lifespan, lifetime.eggs = sum(Count))
fluc9.individual.trait <- summarise(group_by(fluc9.data, Female,Treatment), lifespan = max(Day), bite.rate = sum(Feed)/ lifespan, lifetime.eggs = sum(Count))
fluc12.individual.trait <- summarise(group_by(fluc12.data, Female,Treatment), lifespan = max(Day), bite.rate = sum(Feed)/ lifespan, lifetime.eggs = sum(Count))

#####################
###
#  3) Export files to be used in TraitFitting.R
###
######################

write.csv(constant.individual.trait, "constant.individual.trait.csv")
write.csv(fluc9.individual.trait, "fluc9.individual.trait.csv")
write.csv(fluc12.individual.trait, "fluc12.individual.trait.csv")

