## Erin Mordecai, Stanford University
## August 8, 2018
## Modified by Kerri Miazgowicz on April 13,2020 to generate curve outputs for EIP

## Purpose: Plot trait thermal response comparisons
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) Plot thermal performance curves
##           3) Plot Tmin, Topt, and Tmax mean and 95% CI

##########
###### 1. Set up workspace, load packages, get data, etc.
##########

# Set working directory
mainDir = "C:/Users/Kerri/Desktop/Chapter1 Submission"
setwd(mainDir)

# Load libraties for plotting traits
library(plyr) # Slices and dices data
library(plotrix) # For standard error function
library(coda) #for trait plots

# Load PDR fits
#1. Shapiro PDR
load("saved posteriors inf/shapiro_PDR_briere_inf.Rdata")
EIP.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred #actually fits 1/EIP (aka is PDR)
remove(model.out) #make sure there won't be any errors

#2. Johnson PDR
load("Johnson et al posteriors/PDR.preds.Rdata")
PDR.preds <- PDR.preds

# Load temperature sequence used to calculate trajectories
load("saved posteriors/temps.Rdata")
N.Temp.xs <-length(Temp.xs)

##########
# Create a viewable output file of the TPC generated from EIP (uniform) to be used to generate the TPC for gamma (the porportion of mosquitoes which survive through the latency period)
# Need to use the clacPostQuants function on the EIP.preds
########

##########
###### calcPostQuants - Function to calculate quantiles of derived quantity posteriors (for plotting)
##########

# Arguments:
#   input       data frame with posterior distributions of trait predicted over a temperature gradient (i.e., 'trait.preds')
#   temp.list   the temperature gradient itself (i.e., 'Temp.xs')

calcPostQuants = function(input, temp.list) {
  
  # Get length of gradient
  N.grad.xs <- length(temp.list)
  
  # Create output dataframe
  output.df <- data.frame("mean" = numeric(N.grad.xs), "median" = numeric(N.grad.xs), "lowerCI" = numeric(N.grad.xs), "upperCI" = numeric(N.grad.xs))
  
  # Calculate mean & quantiles
  for(i in 1:N.grad.xs){
    output.df$mean[i] <- mean(input[ ,i])
    output.df$lowerCI[i] <- quantile(input[ ,i], 0.025)
    output.df$upperCI[i] <- quantile(input[ ,i], 0.975)
    output.df$median[i] <- quantile(input[ ,i], 0.5)
  }
  
  output.df # return output
  
}

###########################
##########Part 1 Shapiro.uniform based EIP
########################
Shapiro.EIP.out <- calcPostQuants(EIP.preds, Temp.xs)
#Add in Temp.xs column to the dataframe
Shapiro.EIP.out$temp <- Temp.xs

#Extract out the EIP values for the constant temperature conditions used in my experiment (@temp = 16,20,24,28,32,36) *Remember EIP is fit to 1/EIP(50)
EIP.inf.values <- c(1/Shapiro.EIP.out[Shapiro.EIP.out$temp == 16,1], 1/Shapiro.EIP.out[Shapiro.EIP.out$temp == 20,1],1/Shapiro.EIP.out[Shapiro.EIP.out$temp == 24,1], 1/Shapiro.EIP.out[Shapiro.EIP.out$temp == 28,1],1/Shapiro.EIP.out[Shapiro.EIP.out$temp == 32,1],1/Shapiro.EIP.out[Shapiro.EIP.out$temp == 36,1])
EIP.inf.temp <- c(16,20,24,28,32,36)
EIP.inf.curve.values <- data.frame(EIP = EIP.inf.values, temp = EIP.inf.temp)
EIP.inf.curve.values
#Create an output dataframe with these values
#write.csv(EIP.inf.curve.values,"data/EIP.inf.curve.output.csv", row.names = FALSE)

############################
#######Part 2 Johnson.uniform based EIP
###########################
Johnson.EIP.out <- calcPostQuants(PDR.preds, Temp.xs)
#Add in Temp.xs column to the dataframe
Johnson.EIP.out$temp <- Temp.xs
#Convert values from rate to duration (time, days)
Johnson.EIP.out$EIP <- 1/ Johnson.EIP.out$mean

#Extract out the EIP vlaues for the constant temperature conditions used in my experiment
J.EIP.uniform.values <- c(Johnson.EIP.out[Johnson.EIP.out$temp==16,6], Johnson.EIP.out[Johnson.EIP.out$temp==20,6] ,Johnson.EIP.out[Johnson.EIP.out$temp==24,6],Johnson.EIP.out[Johnson.EIP.out$temp==28,6], Johnson.EIP.out[Johnson.EIP.out$temp==32,6],Johnson.EIP.out[Johnson.EIP.out$temp==36,6])
J.EIP.uniform.temp <- c(16,20,24,28,32,36)
J.EIP.uni.curve.values <- data.frame(EIP = J.EIP.uniform.values, temp = J.EIP.uniform.temp)

#Create an output dataframe with these values
#write.csv(J.EIP.uni.curve.values, "data/J.EIP.uni.curve.values.csv", row.names = FALSE)
