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

summary(lf.norberg)
AIC(lf.norberg)

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
lf.quad$m$getPars() #cf.q= 0.1186192, cf.T0= 1.5491849, cf.Tm= 37.571815

library(nlstools)
#95% confidence interval on model parameters
confint(lf.quad, data =constant.data ) # cf.q = (0.05066386 to 0.1865744);; #cf.T0 (-18.23065968 to 7.5007900); #cf.Tm(36.03025195 to 40.3914812 )


#determining the best fitting model
#Joey/Marta - How do we extract the AIC values from these fits to determine the 'best' functional form for each trait | dtr

#Bootstrapping each model fit for error
#   i.Take a sampling of parameter values using a normal distribution;
#   ii. Generate new curves with the sampled parameter values

n.iterations <- 7500

#generate normal distribution- standard deviation from the 95% CI interval?
#I think the 95% CI is 2 standard deviations away from the mean
sigma <- (40.3914812- 37.57)/2
plot(seq(0,50,0.1),dnorm(seq(0,50,0.1),37.57, sigma)) #Probability function
plot(rnorm(7500,37.57, sigma)) #samples randomly from the specified normal distribution
hist(rnorm(7500,37.57, sigma)) #shows the frequency distribution (new random set from previous code line)

#Generate a code flow through that for the 'best-fitting' model;
# a. stores the mean parameter values
# b. stores the 2.5% lower bound parameter value
# c. stores the 97.5% upper bound parameter value
# d. calculates and stores sigma (from upper bound): (mean - error interval)/2
# d. uses df of a-c to generate 7500 new parameter values; and adds an addition column ('out') for the equation outputs from those parameter values
# e. graph each function output

quad.params <- data.frame("mean"=lf.quad$m$getPars(),"2.5" = confint(lf.quad, data =constant.data )[,1], "97.5" = confint(lf.quad, data =constant.data )[,2] )
colnames(quad.params)<- c("mean", "2.5", "97.5")
quad.params$sigma <- quad.params$`97.5`- quad.params$mean

#Need to do this for each parameter -(have each row be the estimated parameter value; each column as one of the parameters)
#    EXAMPLE: quad.params.boot.list <- rnorm(7500, quad.params$mean, quad.params$sigma)
#Column name becomes previous row name  #still can't figure out how to have variable output be used directly as character during assigment
for( i in 1:nrow(quad.params)){
  quad.params.boot.df[,i] <- data.frame("out" =rnorm(7500, quad.params$mean[i], quad.params$sigma[i]))
  }
colnames(quad.params.boot.df)<- rownames(quad.params)
#remake as a matrix for easier math- if need; if not remove!

#Calculate function output with new parameter values for each iteration (aka row)
#essentially apply QUAD function to each row across our temp gradient to generate 7500 new curves
temp.gradient <- seq(0,50,0.1)

#Below needs to be converted to an apply statement
#should be a [1:501][1:7500] matrix output of the function evaluated across all temp.gradient values for each set of parameters

#initializing the matrix
quad.boot.outputs <- matrix(nrow = 501, ncol = 7500)

#Its using the same parameter values for each call-> the function isnt taking in different values...
#worked for the first column.....(now maybe for everything....) 
for(l in 1: length(temp.gradient)){
  quad.boot.outputs[l,] <- QUAD(cf.q = test[,1], cf.T0 = test[,2], cf.Tm = test[,3], Treatment = temp.gradient[l])
#plot all of the new functions with the mean and error functions
  print(l)
  } # end  loop (columns)

plot.new()
plot(temp.gradient, quad.boot.outputs[,1], ylim = c(0, 60), type = 'n')
for(j in 1:ncol(quad.boot.outputs)){
  lines(temp.gradient, quad.boot.outputs[,j], col=alpha(rgb(0,0,0), 0.1))
}
#show the mean on top
quad.outputs.means <- QUAD(cf.q = quad.params$mean[1], cf.T0 =quad.params$mean[2] , cf.Tm = quad.params$mean[3], Treatment = temp.gradient )
lines(temp.gradient, quad.outputs.means, col = "red")

#above is for visualization

#calculate the 95% HPD of the parameter values and mean/median from the bootstrapping and use those to show the error around the curve

