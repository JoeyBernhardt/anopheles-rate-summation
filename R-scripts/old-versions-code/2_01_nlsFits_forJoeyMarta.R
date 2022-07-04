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

#call nls function over each dataset
lf.quad <- nls(lifespan ~(-1 * cf.q * (Treatment - cf.T0) * (Treatment - cf.Tm)),
               data= constant.data,
               start = list(cf.q=0.15, cf.T0=4.7, cf.Tm=36.5))
plot(constant.data$Treatment, constant.data$lifespan)
lines(seq(0,50,0.1),predict(lf.quad,data.frame(Treatment= seq(0,50,0.1))))

#add in error around the parameter estimates
lf.quad2 <- nls(lifespan ~ QUAD(cf.q, cf.T0, cf.Tm, Treatment),
               data= constant.data,
               start = list(cf.q=0.15, cf.T0=4.7, cf.Tm=36.5))
plot(constant.data$Treatment, constant.data$lifespan)
lines(seq(0,50,0.1),predict(lf.quad2,data.frame(Treatment= seq(0,50,0.1))))

#extracting model fit parameters
lf.quad$m$getPars() #cf.q- 0.1186192, cf.T0- 1.5491849, cf.Tm- 37.571815
#Joey- How do I get the parameter error interval???
