###Code written by Kerri Miazgowicz July 30 2020
################################# Contents
# 1. Set-up + data vis 
# 2. NLS fits
################################# 1. Set-up + data vis

###### Load packages
library(plyr)
library(dplyr)
library(tidyverse)

#Set working directory
mainDir = "C:/Users/Kerri/Desktop/Chapter2 NLS"
setwd(mainDir)

######  Load raw trait data
alldata <- data.frame(read.csv("data-raw/fluc.program.trait.means.csv", colClasses = c("DTR" ="character")))

#change default plot specifications
par(mar=c(1,1,1,1))
dev.off()

#example for how to seperate out the data for each DTR
constant.data <- subset(alldata, DTR == "constant")
dtr9.data<- subset(alldata, DTR == "+/- 4.5")
dtr12.data<- subset(alldata, DTR == "+/- 6")


#Can I pre-define each function and call it in the nls fitting - YES
#Define each function of interest
QUAD <- function(cf.q,cf.T0,cf.Tm, Treatment){
  return(-1 * cf.q * (Treatment - cf.T0) * (Treatment - cf.Tm))
}
BRIERE <- function(cf.q, cf.T0,cf.Tm,Treatment){
  return( cf.q * Treatment * (Treatment - cf.T0) * sqrt((cf.Tm - Treatment)*(cf.Tm > Treatment)))
}
NORBERG <- function(cf.a, cf.b, cf.w, cf.z, Treatment){
  return( cf.a * exp(cf.b * Treatment) * (1-((Treatment-cf.z)/(cf.w/2))^2) )
}
MOD_SS <- function(cf.c, cf.Topt, cf.ED,TempK){ #Needs to be fit with temperature in Kelvin
  return(cf.c * 10^11 * exp(-0.6/(8.62*10^-5*TempK)) / (1 + exp(-1/(8.62*10^-5*TempK) * (cf.ED - TempK*(cf.ED/cf.Topt + 8.62*10^-5*log(0.6/(cf.ED-0.6)) ) ) ) ) )
}

#call nls function over each dataset

#add in error around the parameter estimates
lf.quad <- nls(lifespan ~ QUAD(cf.q, cf.T0, cf.Tm, Treatment),
               data= constant.data,
               start = list(cf.q=0.15, cf.T0=4.7, cf.Tm=36.5))
plot(constant.data$Treatment, constant.data$lifespan, xlim= c(0,50))
lines(seq(0,50,0.1),predict(lf.quad,data.frame(Treatment= seq(0,50,0.1))))

lf.briere <- nls(lifespan ~ BRIERE(cf.q, cf.T0,cf.Tm,Treatment),
               data= constant.data,
               start = list(cf.q=0.03, cf.T0=3.8, cf.Tm=29.6))
plot(constant.data$Treatment, constant.data$lifespan, xlim= c(0,50))
lines(seq(0,50,0.1),predict(lf.briere,data.frame(Treatment= seq(0,50,0.1))))

lf.norberg <- nls(lifespan ~ NORBERG(cf.a, cf.b, cf.w, cf.z, Treatment),
                  data = constant.data,
                  start = list(cf.a=17, cf.b=0.03, cf.w=28, cf.z=21))
plot(constant.data$Treatment, constant.data$lifespan, xlim= c(0,50))
lines(seq(0,50,0.1),predict(lf.norberg,data.frame(Treatment= seq(0,50,0.1))))

#Fit this in temp Kelvin
lf.modSS <- nls(lifespan ~ MOD_SS(cf.c, cf.Topt, cf.ED, TempK),
                data = constant.data,
                start = list(cf.c=10 , cf.Topt=294, cf.ED = 2))
plot(constant.data$Treatment, constant.data$lifespan, xlim= c(0,50))
lines(seq(0,50,0.1),predict(lf.modSS,data.frame(TempK= seq(273.15,273.15+50,0.1))))

###FIT ANY ADDITIONAL FUNCTIONS HERE

#extracting model fit parameters
lf.quad$m$getPars() #cf.q- 0.1186192, cf.T0- 1.5491849, cf.Tm- 37.571815

#determining the best fitting model
#Joey/Marta - How do we extract the AIC values from these fits to determine the 'best' functional form for each trait | dtr
#Joey- How do I get the parameter error interval???
