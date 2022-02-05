## Kerri Miazgowicz, University of Georgia
## June 6, 2019
## 
## Purpose: Generate Population and replicate-specific values for r, Tc, and R0(mosqutio fitness) for constant and fluctuating temperatures
## 
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) Import block-specific Kaplan-Meier survival estimates for lx (we have right censored data)from 'Analysis_survival_fluc.rmd'
##           3) Calculate r / Tc / R0 across all datasets and for each 'Block' within a dataset
##           4) Export files to be used in TraitFitting.R
##                  i. 'constant.population.R0_Tc_r'
##                  ii. 'dtr9.population.R0_Tc_r'
##                  iii.  'dtr12.population.R0_Tc_r'        


#####################
###
#  1) Set-up,load packages, get data, etc.
###
######################
library(dplyr)

#Set working directory
mainDir = "C:/Users/Kerri/Desktop/RateSummationProject"
setwd(mainDir)

#Read in csv files for each dataset
constant.data <- read.csv("data/constant_master.csv")
constant.data[is.na(constant.data)] <- 0
constant.data$Program <- "constant"
constant.data$Censor <- ifelse(constant.data$Dead == 2, 1, 0) #Add a censor column for Sx calculation
constant.data$Dead[constant.data$Dead == 2] <- 0 #Change censored individuals to 0 to not mess with Sx calc

fluc9.data <- read.csv("data/fluc_dtr9_master.csv")
fluc9.data[is.na(fluc9.data)] <- 0
fluc9.data$Program <- "dtr9"
fluc9.data$Censor <- ifelse(fluc9.data$Dead == 2, 1, 0) 
fluc9.data$Dead[fluc9.data$Dead == 2] <- 0 #Change censored individuals to 0 to not mess with Sx calc

fluc12.data <- read.csv("data/fluc_dtr12_master.csv")
fluc12.data[is.na(fluc12.data)] <- 0
fluc12.data$Program <- "dtr12"
fluc12.data$Censor <- ifelse(fluc12.data$Dead == 2, 1, 0) 
fluc12.data$Dead[fluc12.data$Dead == 2] <- 0 #Change censored individuals to 0 to not mess with Sx calc

all.data <- rbind(constant.data,fluc9.data, fluc12.data)
#####################
###
#  2) Import block-specific Kaplan-Meier survival estimates for lx (we have right censored data)from 'Analysis_survival_fluc.rmd'
###
######################

#Place all datasets of kaplan-meier estimates into a list
###################
dir <- "data/"
file_raw <- list.files(path=dir,pattern="*values.csv") #double check there are no other datasets with this ending
file_list <- paste(dir, file_raw, sep="")
data.list <- list()

for (i in 1:length(file_list))
{
  file <- file_list[i]  #df for each csv file of kaplan-meier estimates
  dataset <- read.table(file, header=T, sep=",")
    #change names of headers to day and survival
    names(dataset) <- c("day", "survival")
    #add program name
    program_ID <- sub("*.csv", "", file_raw[i])
    dataset$program_ID <- program_ID
    dataset$day <- as.numeric(dataset$day)
    data.list[[i]] <- dataset
  }
  
updated.data.list<- list()
#Fill in day gaps
for(y in 1:length(data.list)){
  df <- data.list[[y]]
   for(z in 1:max(df$day)){
    if(z == df$day[z+1]){
      
    }else{
      new.row <- data.frame(day = z,survival = as.numeric(df$survival[z]), program_ID =df$program_ID[z])
      df <- rbind(df, new.row)
      df <- df[order(df$day),]
    } #end else statement
    updated.data.list[[y]] <- df
  }#end z for loop
}#end y for loop

#Now survival for all updated.data.list corresponds to lx
#instead of a list; place into an extended dataframe ('Treatment', 'Block','Survival', 'Day','Program')
for(x in 1:length(updated.data.list)){
  df<- updated.data.list[[x]]
  names(df)<- c('Day', "survival", "ID")
  df$Treatment <- as.integer(substr(df$ID,5,6)) #string subset from "ID"
  df$Block <- substr(df$ID, regexpr("B", df$ID) + 1,regexpr("B", df$ID)+ 1 )
  #Note: block from survival doesnt match raw block # because of how the surival estimates were named in the code....
  if(substr(df$ID, regexpr("R", df$ID) + 1,regexpr("R", df$ID)+ 1 )== "0"){df$Program <- "constant"}
  if(substr(df$ID, regexpr("R", df$ID) + 1,regexpr("R", df$ID)+ 1 )== "9"){df$Program <- "dtr9"}
  if(substr(df$ID, regexpr("R", df$ID) + 1,regexpr("R", df$ID)+ 1 )== "1"){df$Program <- "dtr12"}
  if(x ==1){
    all.KM.survival <- df
  }else{
    all.KM.survival <- rbind(all.KM.survival, df)
    
  } #ends else
}#ends x loop
all.KM.survival$Block <- as.integer(all.KM.survival$Block)
#subset all.KM.survival to each 'Program'
constant.KM.survival <- subset(all.KM.survival, Program == "constant")
dtr9.KM.survival <- subset(all.KM.survival, Program == "dtr9")
dtr12.KM.survival <- subset(all.KM.survival, Program == "dtr12")


#check that the max number of days corresponds
#update this lx to the constant.rep.daily summaries and dtr9; and dtr12 datasets below

#####################
###
#  3) Calculate r / Tc / R for each Block across all datasets
###
######################
#Note: mx = the estimated amount of female offspring; assumes a 1:1 female to male ratio of progeny
#Censorship is not accounted for properly -> use km survival estimates instead for lx 
#This lx column will be added and called lx.KM

#To avoid a loop to tally all calcs since I have censored individuals; I use cumsum() over grouped variables
censor.day.counts <- mutate(group_by(constant.data, Day, Block, Treatment), total.day.censor= cumsum(Censor))%>% ungroup
constant.data$T.censor <- censor.day.counts$total.day.censor

#Create new dataframe for each temperature|Day|Block profile for the aggregated data
constant.rep.daily <- summarise(group_by(constant.data, Block, Treatment, Day), Sx = n()-sum(Dead)-sum(Censor), dead = sum(Dead),censor = sum(Censor), lx = Sx / (30-max(T.censor)), daily.eggs = sum(Count), mx = 0.5*daily.eggs, lxmx = lx *mx)%>% ungroup

#Create new dataframe fro each temperature|Day|Block profile but includes the lx from the Kaplan-Meier survival since I have censored individuals
####
#step 1: add right lx values to the Temp|Block combo
#step 2: merge new lx values with constant.rep.daily
constant.rep.daily <- merge(constant.rep.daily, constant.KM.survival, by = c("Treatment", "Block","Day"))
#step 3: calculate mxlx and other values
constant.rep.daily$mxsurvival <- constant.rep.daily$mx * constant.rep.daily$survival
####
constant.rep.daily$mxsurvivalx <- constant.rep.daily$mxsurvival * constant.rep.daily$Day
constant.rep.sums <- summarise(group_by(constant.rep.daily, Block, Treatment), R0 = sum(mxsurvival), Tc = sum(mxsurvivalx)/sum(mxsurvival), est.r = log(R0)/Tc)%>%ungroup              
constant.rep.sums$Program <- "constant"

#I do NOT solve for r by iteration with the Euler equation as I censor individuals at some treatments and thus survival!=0 ; therefore cumsum() !=1

#repeat above with the fluctuation data
dtr9.rep.daily <- summarise(group_by(fluc9.data, Block, Treatment, Day), daily.eggs = sum(Count), mx = 0.5*daily.eggs)%>% ungroup
dtr9.rep.daily <- merge(dtr9.rep.daily, dtr9.KM.survival, by = c("Treatment", "Block","Day"))
dtr9.rep.daily$mxsurvival <- dtr9.rep.daily$mx * dtr9.rep.daily$survival
dtr9.rep.daily$mxsurvivalx <- dtr9.rep.daily$mxsurvival * dtr9.rep.daily$Day
dtr9.rep.sums <- summarise(group_by(dtr9.rep.daily, Block, Treatment), R0 = sum(mxsurvival), Tc = sum(mxsurvivalx)/sum(mxsurvival), est.r = log(R0)/Tc)%>%ungroup              
dtr9.rep.sums$Program <- "dtr9"

#repeat above with the fluctuation data
dtr12.rep.daily <- summarise(group_by(fluc9.data, Block, Treatment, Day), daily.eggs = sum(Count), mx = 0.5*daily.eggs)%>% ungroup
dtr12.rep.daily <- merge(dtr12.rep.daily, dtr12.KM.survival, by = c("Treatment", "Block","Day"))
dtr12.rep.daily$mxsurvival <- dtr12.rep.daily$mx * dtr12.rep.daily$survival
dtr12.rep.daily$mxsurvivalx <- dtr12.rep.daily$mxsurvival * dtr12.rep.daily$Day
dtr12.rep.sums <- summarise(group_by(dtr12.rep.daily, Block, Treatment), R0 = sum(mxsurvival), Tc = sum(mxsurvivalx)/sum(mxsurvival), est.r = log(R0)/Tc)%>%ungroup              
dtr12.rep.sums$Program <- "dtr12"

#Plot of all population metrics
all.pop.metrics <- rbind(constant.rep.sums, dtr9.rep.sums, dtr12.rep.sums)

library(ggplot2)
ggplot(all.pop.metrics, aes(x = Treatment, y = R0, group = Program))+
  geom_point(aes(colour = Program))

ggplot(all.pop.metrics, aes(x = Treatment, y = Tc, group = Program))+
  geom_point(aes(colour = Program))

ggplot(all.pop.metrics, aes(x = Treatment, y = est.r, group = Program))+
  geom_point(aes(colour = Program))


#####################
###
#  3) Export files to be used in TraitFitting.R or the file name we end up using
###
######################

write.csv(constant.rep.sums, "data/constant.population.R0_Tc_r.csv")
write.csv(dtr9.rep.sums, "data/dtr9.population.R0_Tc_r.csv")
write.csv(dtr12.rep.sums, "data/dtr12.population.R0_Tc_r.csv")

