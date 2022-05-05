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
alldata <- data.frame(read.csv("fluc.program.trait.means.csv", colClasses = c("DTR" ="character")))

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

#add in error around the parameter estimates later
EFD.quad <- nls(EFD ~ QUAD(cf.q, cf.T0, cf.Tm, Treatment),
               data= constant.data,
               start = list(cf.q=0.15, cf.T0=4.7, cf.Tm=36.5))
plot(constant.data$Treatment, constant.data$EFD, xlim= c(0,50),ylim = c(-1,20))
lines(seq(0,50,0.1),predict(EFD.quad,data.frame(Treatment= seq(0,50,0.1))))

EFD.briere <- nls(EFD ~ BRIERE(cf.q, cf.T0,cf.Tm,Treatment),
                 data= constant.data,
                 start = list(cf.q=0.03, cf.T0=3.8, cf.Tm=29.6))
plot(constant.data$Treatment, constant.data$EFD, xlim= c(0,50),ylim = c(-1,20))
lines(seq(0,50,0.1),predict(EFD.briere,data.frame(Treatment= seq(0,50,0.1))))

EFD.norberg <- nls(EFD ~ NORBERG(cf.a, cf.b, cf.w, cf.z, Treatment),
                  data = constant.data,
                  start = list(cf.a=.5, cf.b=-0.1, cf.w=80, cf.z=60))
plot(constant.data$Treatment, constant.data$EFD, xlim= c(0,50),ylim = c(-1,20))
lines(seq(0,50,0.1),predict(EFD.norberg,data.frame(Treatment= seq(0,50,0.1))))

#Fit this in temp Kelvin
EFD.modSS <- nls(EFD ~ MOD_SS(cf.c, cf.Topt, cf.ED, TempK),
                data = constant.data,
                start = list(cf.c=10 , cf.Topt=280, cf.ED = 2))
plot(constant.data$Treatment, constant.data$EFD, xlim= c(0,50),ylim = c(-1,20))
lines(seq(0,50,0.1),predict(EFD.modSS,data.frame(TempK= seq(273.15,273.15+50,0.1))))

###FIT ANY ADDITIONAL FUNCTIONS HERE

#determining the best fitting model
AIC(br.quad) #-58.837
AIC(br.briere) #-54.71379
AIC(br.modSS) #71.08535
AIC(br.norberg) #-57.02561

#Norberg likely will behave the best in rate summation- proceed with bootstrapping on that one

#extracting model fit parameters
br.briere$m$getPars() #cf.q= 0.0001299665 cf.T0= -2.6744415063 cf.T0= 43.3665526574 
library(nlstools)

#95% confidence interval on model parameters
br.constant.briere.boot <- nlsBoot(br.briere, niter = 1000) #bootstap

plot(br.constant.briere.boot) #generic plot

#shows the spread of the parameter values
plot(br.constant.briere.boot, type = c("boxplot"),
     mfr  = c(ceiling(sqrt(ncol(br.constant.briere.boot$coefboot))),
              ceiling(sqrt(ncol(br.constant.briere.boot$coefboot)))),
     ask=FALSE,)

summary(br.constant.briere.boot) #Provides estimates and error; also provides the median and 95% confidence interval

#Plot the summary metrics for the norberg bootstrap from the 'summary()' using median here
#Median of bootstrap estimates and percentile confidence intervals
#           Median          2.5%        97.5%
#cf.q   0.0001305231  7.937853e-05 1.799421e-04
#cf.T0 -2.6434635250 -1.465596e+01 3.553586e+00
#cf.Tm 43.3584191747  4.092300e+01 4.855576e+01

#See if there is a way to code that directly later; for now create table to be used during plotting
boot.curves<- data.frame("median" = c(0.0001305231,-2.6434635250,43.3584191747),"low.2.5" =c(7.937853e-05,-1.465596e+01,4.092300e+01), "upper.97.5"=c(1.799421e-04,3.553586e+00,4.855576e+01))
row.names(boot.curves)<- c("cf.q","cf.T0","cf.Tm")

#derive curve outputs using the parameter values-> will need these for use later anyways in RS
temp.gradient <- seq(0,50,0.1)
briere.outputs.medians <- BRIERE(cf.q = boot.curves$median[1], cf.T0 =boot.curves$median[2] , cf.Tm = boot.curves$median[3], Treatment = temp.gradient )
briere.outputs.lowerCI <- BRIERE(cf.q = boot.curves$low.2.5[1], cf.T0 =boot.curves$low.2.5[2] , cf.Tm = boot.curves$low.2.5[3], Treatment = temp.gradient )
briere.outputs.upperCI <- BRIERE(cf.q = boot.curves$upper.97.5[1], cf.T0 =boot.curves$upper.97.5[2] , cf.Tm = boot.curves$upper.97.5[3], Treatment = temp.gradient )

#Plot-
plot(constant.data$Treatment, constant.data$bite.rate, xlim= c(0,50), ylim=c(-0.1,1))
lines(seq(0,50,0.1),predict(br.briere,data.frame(Treatment= seq(0,50,0.1))), col="red")
abline(h=0, col ="black", lty=2)
#add boot data
lines(temp.gradient,briere.outputs.medians, col ="blue")
lines(temp.gradient,briere.outputs.lowerCI, col= "blue", lty =2)
lines(temp.gradient,briere.outputs.upperCI, col="blue",lty=2)

