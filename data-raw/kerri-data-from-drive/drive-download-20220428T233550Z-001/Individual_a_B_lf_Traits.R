## Kerri Miazgowicz, University of Georgia
## October 1, 2018
## Stats included October 1, 2019
## Rerun 4/22/2020
## 
## Purpose: Generate indiviudal-specific values for lifespan, daily biting rate, and lifetime egg production for constant and fluctuating temperatures
## ID / temp / a / lf / B <- from master raw csv files
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) Calculate a / lf / B across all datasets
##           3) Export files to be used in TraitFitting.R -> provided to Joey and Marta for fitting in code file ''.
##           4) Perform GLMM models on a / lf/ B across all datasets 
##        


#####################
###
#  1) Set-up,load packages, get data, etc.
###
######################
library(dplyr)

#Set working directory
mainDir = "C:/Users/Kerri/Desktop/Chapter2 InProgress"
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

constant.individual.trait <- constant.data %>%
                             dplyr::group_by(Female,Treatment,Block)%>%
                             dplyr::summarise(lifespan = max(Day), bite.rate = sum(Feed)/ lifespan, lifetime.eggs = sum(Count))
                             
fluc9.individual.trait <- fluc9.data %>%
                          dplyr::group_by(Female,Treatment, Block)%>%
                          dplyr::summarise(lifespan = max(Day), bite.rate = sum(Feed)/ lifespan, lifetime.eggs = sum(Count))


fluc12.individual.trait <- fluc12.data%>% 
                            dplyr::group_by(Female,Treatment, Block)%>%
                            dplyr::summarise(lifespan = max(Day), bite.rate = sum(Feed)/ lifespan, lifetime.eggs = sum(Count))




#####################
###
#  3) Export files to be used in TraitFitting.R
###
######################

#write.csv(constant.individual.trait, "data/constant.individual.trait.csv", row.names = F)
#write.csv(fluc9.individual.trait, "data/fluc9.individual.trait.csv", row.names = F)
#write.csv(fluc12.individual.trait, "data/fluc12.individual.trait.csv", row.names =F)

#####################
###
#  4) Perform GLMM models on a / lf/ B across all datasets 
### Note: Need to keep block variable for random effects; won't be able to include donor or day due to summary data structure
#####################
#

#Repeat above summary, but include Block for analysis
#create new dataframe for each temperature profile for the aggregated data

constant.individual.trait <- as.data.frame(summarise(group_by(constant.data, Female,Treatment, Block), lifespan = max(Day), bite.rate = sum(Feed)/ lifespan, lifetime.eggs = sum(Count)))
fluc9.individual.trait <- as.data.frame(summarise(group_by(fluc9.data, Female,Treatment, Block), lifespan = max(Day), bite.rate = sum(Feed)/ lifespan, lifetime.eggs = sum(Count)))
fluc12.individual.trait <- as.data.frame(summarise(group_by(fluc12.data, Female,Treatment, Block), lifespan = max(Day), bite.rate = sum(Feed)/ lifespan, lifetime.eggs = sum(Count)))

#replace "block" and "individual" values within each dataset with '.DTR0, 9, or 12' so that when I merge all data together those stay seperate; will become string value
#convert block to a string variable and then add a character segment to it
constant.individual.trait$Block <- paste0(as.character(constant.individual.trait$Block),".dtr0")
constant.individual.trait$Female <- paste0(as.character(constant.individual.trait$Female),".dtr0")
constant.individual.trait$Program <- as.character("dtr0")
fluc9.individual.trait$Block <- paste0(as.character(fluc9.individual.trait$Block), ".dtr9")
fluc9.individual.trait$Female <- paste0(as.character(fluc9.individual.trait$Female), ".dtr9")
fluc9.individual.trait$Program <- as.character("dtr9")
fluc12.individual.trait$Block <- paste0(as.character(fluc12.individual.trait$Block), ".dtr12")
fluc12.individual.trait$Female <- paste0(as.character(fluc12.individual.trait$Female), ".dtr12")
fluc12.individual.trait$Program <- as.character("dtr12")

#merge all data into one dataframe
all.individual.trait <- rbind(constant.individual.trait, fluc9.individual.trait, fluc12.individual.trait)

#perform GLMM analysis
#Main effects = Treatment(scaled & centered), Program (categorical)
#Random = Block (categorical), Female(categorical)
#Response variables = lifespan(continuous/count), bite.rate(rate), lifetime.eggs(continuous/count)

#Following Courtney's previous code[Tesla temperature], making Temperature a continuous variable
#feeddata$Treatment <- as.factor(feeddata$Treatment)
all.individual.trait$Female <- as.factor(all.individual.trait$Female)
all.individual.trait$Block <- as.factor(all.individual.trait$Block)
all.individual.trait$Program <- as.factor(all.individual.trait$Program)

#Centering and scaling continuous variables of interest (Treatment)
all.individual.trait$Tscale <- scale(all.individual.trait$Treatment)

library(fitdistrplus)
library(lme4)
library(car)



#Looking at my data distribution
hist(all.individual.trait$lifespan, breaks=50, col="red") # right-skewed
hist(all.individual.trait$bite.rate, breaks=50, col="red") #normal between 0 and 1
hist(all.individual.trait$lifetime.eggs, breaks=50, col="red") #0-inflated nbinomial
#descdist(all.individual.trait$lifetime.eggs, discrete = FALSE) # This suggests to use a uniform distribution

###LIFESPAN####
fit.gamma <- fitdist(all.individual.trait$lifespan, distr = "gamma", method = "mme") #gamma fits fine
plot(fit.gamma)

lifespan.temp.model <- glm(lifespan ~ Tscale*Program, family = Gamma(link ='inverse'), data = all.individual.trait)
Anova(lifespan.temp.model)
summary(lifespan.temp.model)

#model does not converge if I include any random factors
lf.m2 <- glmer(lifespan ~ Tscale*Program + (1|Block), family = Gamma(link = "inverse"), data = all.individual.trait)

###################3

###LIFETIME EGGS####
fit.nbinom <- fitdist(all.individual.trait$lifetime.eggs, distr = "nbinom", method = "mme") #gamma fits fine
plot(fit.nbinom)

summary(m1 <- glm.nb(lifetime.eggs ~ Tscale*Program, data = all.individual.trait))

lifeeggs.temp.model <- glm(lifetime.eggs ~ Tscale*Program, family = poisson(link ='log'), data = all.individual.trait)
Anova(lifeeggs.temp.model)
summary(lifeeggs.temp.model)

lf.m2 <- glmer(lifetime.eggs ~ Tscale*Program + (1|Block), family = poisson(link = "log"), data = all.individual.trait)
AIC(lf.m2)
Anova(lf.m2)
summary(lifeeggs.temp.model)

library(pscl)
zero.infl1 <- zeroinfl(lifetime.eggs ~ Tscale*Program, data = all.individual.trait) #runs zero-inflated pois
zero.infl2 = zeroinfl(lifetime.eggs ~ Tscale*Program, data = all.individual.trait, dist = "negbin") #runs zero-inflated nb

summary(zero.infl1) #zero inflated seems to be working well
summary(zero.infl2) #zero inflated seems to be working well
Anova(zero.infl2)
AIC(zero.infl2)
#####################

###Biting rate########
fit.gamma2 <- fitdist(all.individual.trait$bite.rate, distr = "gamma", method = "mme")
#to create no non-postive values
all.individual.trait$bite.rate.t <- all.individual.trait$bite.rate + 0.0001
plot(fit.gamma2)
biterate.temp.model <- glm(bite.rate.t ~ Tscale*Program, family =Gamma(link = 'inverse'), data = all.individual.trait)
Anova(biterate.temp.model)
summary(biterate.temp.model)
