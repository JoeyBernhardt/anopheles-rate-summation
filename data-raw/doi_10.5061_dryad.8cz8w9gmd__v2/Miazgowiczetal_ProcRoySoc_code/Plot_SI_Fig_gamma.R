## Erin Mordecai, Stanford University
## August 8, 2018
## Updated by Kerri Miazgowicz on November 19th, 2019 to include the truncated model fits.
## Updated by KM on April 13, 2020 to revise for resubmission
## Generate the plot specifically comparing gamma across datasets
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

# Check whether there's a folder in the directory for saving plots
# If not, create one
subDir = "plots"
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)

# Load libraties for plotting traits
library(plyr) # Slices and dices data
library(dplyr)
library(plotrix) # For standard error function
library(coda) #for trait plotting

#Load trait data the curve was fit over
data.gamma <- read.csv("data/Miaz_gamma.values.csv")

# Load temperature sequence used to calculate trajectories
load("saved posteriors/temps.Rdata")
N.Temp.xs <-length(Temp.xs)

############ Load traits fits I need
#1.An.stephensi lifetime model; gamma [Miaz gompertz, shapiro EIP]
#2.An.stephensi estimated model; e^() [Miaz exp mu, shapiro EIP]
#3.Multi-species estimated model; e^() [Johnson exp mu, Johnson EIP]

#Gamma direct fit
#no informative priors can be added to this one
load("saved posteriors/gamma_quad_withT_all_uniform.Rdata")
gamma.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
gamma.params <- model.out$BUGSoutput$summary[1:5,1:9]



####Make a list of all of the trait function parameters and export them for quick reference
param.compile.fun = function(txt) {
  out = list()
  for (i in 1:length(txt)){
    df = eval(parse(text = txt[i]))
    out[[i]]<- df
  }
  return(out)
}
#looking at thermal limits for estimated.lifeeggs
all.param.list = list(c("gamma.params"))
trait.params = lapply(all.param.list, param.compile.fun)[[1]]
names(trait.params) <- all.param.list[[1]]
trait.params

#Export this for quick reference
#sink("output_final_gamma_TPC_params.txt")
#print(trait.params)
#sink()


#Miaz exp mu
load("saved posteriors inf/mu_exponential_quadratic_inf.Rdata")
mu.exp.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred  # actually fits for 1/mu

#Shapiro traits -> for determining the TPC thresholds
load("saved posteriors inf/shapiro_PDR_briere_inf.Rdata")
EIP.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred  # actually fits for 1/EIP50

#Johnson traits
# Fits from Johnson et al 2015
fits = paste("Johnson et al posteriors/", c(
  "PDR.preds.Rdata",
  "mu.preds.Rdata"), sep = "")
for (i in 1:length (fits)) load(fits[i])

#############
####### 2. Plot thermal performance curves
#############


###########
######
# Compare different ways of estimating trait TPCs
######
###########


##########
###### 2. calcPostQuants - Function to calculate quantiles of derived quantity posteriors (for plotting)
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


##########
###### 4. Function to calculate distributions of T0, Tmax, and peak R0
##########

# Arguments:
#   input       data frame with posterior distributions of R0 (or other quantity) predicted over a temperature gradient
#   temp.list   the temperature gradient itself (i.e., 'Temp.xs')

calcThreshPeakDists = function(input, temp.list) {
  
  # Create output dataframe
  output.df <- data.frame("peak" = numeric(nrow(input)), "T0" = numeric(nrow(input)), "Tmax" = numeric(nrow(input)))
  
  for (i in 1:nrow(input)) { # loop through each row of the input (MCMC step)
    
    output.df$peak[i] <- temp.list[which.max(input[i, ])] # Calculate index of R0 peak & store corresponding temperature
    
    # Makes lists of T0 and Tmax
    index.list <- which(input[i, ] > 0) # Create vector of list of indices where R0 > 1
    length.index.list <- length(index.list)
    output.df$T0[i] <- temp.list[index.list[1] - 1] # Store T0 (index prior to first value in index.list)
    output.df$Tmax[i] <- temp.list[index.list[length.index.list] + 1] # Store Tmax (index after to last value in index.list)
  }
  
  output.df # return
  
}

#######
#To be able to properly compare gammas (direct vs. indirect)
# I need to create a function that generates this new curve and the associated CI - similar to how it is done with R0
# direct gamma = gamma used in An.stephensi lifetime model
# indirect gamma = calculated gamma from e^()

### function to define indirect gamma for An.stephensi estimated model and Multi-species estimated model:
# constant to keep lifespan from being numerically zero
# assume minimum survival time is half an hour
ec.lf = 1/48
# constant to keep PDR from being numerically zero
# assume the maximum EIP is 100 days
ec.pdr = 1/100

indirect.Miaz.gamma.fun <- function(lf,PDR){  #takes in lf as lifespan
  value <-exp(-(1/(lf+ec.lf))*(1/(PDR+ec.pdr)))
  return(value)
}

indirect.Johnson.gamma.fun <- function(mu,PDR){  #takes in mu as 1/lf
  lifespan = 1/mu
  value <- exp(-(1/(lifespan+ec.lf))*(1/(PDR+ec.pdr)))
  return(value)
}


######  Calculate gamma posteriors and quantiles
gamma.out = calcPostQuants(gamma.preds, Temp.xs)
gamma.out$temp <- Temp.xs
gamma.TPdists <- calcThreshPeakDists(gamma.preds,Temp.xs)

# Calculate indirect gamma full posteriors and quantiles
#An.stephensi estimated model
estimated.gamma <- indirect.Miaz.gamma.fun(mu.exp.preds, EIP.preds) #mu.exp.preds as lifespan, #EIP.preds as PDR
estimated.gamma.out = calcPostQuants(estimated.gamma, Temp.xs)
estimated.gamma.out$temp = Temp.xs
estimated.gamma.TPdists = calcThreshPeakDists(estimated.gamma, Temp.xs)

#Multi-species estimated model
J.est.gamma <- indirect.Johnson.gamma.fun(mu.preds, PDR.preds) #mu.preds as mu, PDR.preds as PDR
J.est.gamma.out <- calcPostQuants(J.est.gamma, Temp.xs)
J.est.gamma.out$temp<- Temp.xs
J.est.gamma.TPdists = calcThreshPeakDists(J.est.gamma, Temp.xs)
##############
#######################################################################
##############Plotting functions
######################################################################
###### Plot TPCs with the raw data
trait.plots = function(traitname, fitname, df, labs, pch = 1, cols = c("black"),source = "inf", add ="false"){
  # Specify the temperature sequences and trait to plot
  par(mar=c(5,6,5,1))
  temp = df$temp
  traitplot = df[ , which(colnames(df)==traitname)]
  ifelse(source == "inf",load(paste("saved posteriors inf/", fitname, sep = "")),load(paste("saved posteriors/", fitname, sep = "")))
  model.fit = model.out
  
  # plot the data
  if(add== "false"){
    plot(traitplot ~ jitter(temp, 0.5), xlim = c(0, 45), ylim = c(0, max(traitplot, na.rm = T)*1.02), 
         ylab = labs, xlab = expression(paste("Temperature (",degree,"C)")), cex.lab = 1.4, cex.axis = 1.4, pch = pch)
  }else{
    points(traitplot ~ temp, data = df, pch = pch, col = cols)
  }
  # plot the TPC mean and 95% credible interval
  lines(model.fit$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = cols)
  lines(model.fit$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = cols)
  lines(model.fit$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd =2, col = cols)
}

# Alternative plotting function to display CI as dashed lines for the redefined data
plot.trait2 = function(temp = Temp.xs, df, modname = "", add = "false", cols = 1){
  if (add == "false") {
    par(mar=c(5,6,5,1))
    plot(temp, df$upperCI, type = "l", lty = 2, lwd = 1, xlab = expression(paste("Temperature (",degree,"C)")), ylab = modename, cex.axis = 1.2, cex.lab = 1.2, main = modname, cex.main = 1.5, yaxt = 'n', col = cols)
  } else {
    lines(temp, df$upperCI, lty = 2, lwd = 1, col = cols)
  }
  lines(temp, df$mean, lty = 1, lwd = 2, col = cols)
  lines(temp, df$lowerCI, lty = 2, lwd = 1, col = cols)
}


##########Supplemental Figure 4 - gamma overlay#############

pdf("plots/SI_gamma_overlay_.pdf",  width = 7, height = 5)
trait.plots("gamma", "gamma_quad_withT_all_uniform.Rdata", data.gamma, "Gamma",cols="black", source = "uni",pch =16)
# since we have to combine fits from multiple traits
plot.trait2(df = estimated.gamma.out, add = "true", cols = "dodgerblue")
plot.trait2(df = J.est.gamma.out, add = "true", cols = "purple")
legend('topleft', legend = c(expression(italic(An.~stephensi)~lifetime),expression(italic(An.~stephensi)~estimated), "Multi-species estimated"), lty = c(1,1,1), lwd = 2, col = c("black", "dodgerblue","purple"), bty = 'n', cex=0.8, pt.cex = 1)
dev.off()

########################################
#######################################3
#########################################



############
####### 3. Plot Tmin, Topt, and Tmax mean and 95% CI
###########

# Function to summarize the median and 95% CI
summary.fn = function(df){
  med = median(df, na.rm = T)
  int = c(HPDinterval(mcmc(df)))
  out = c(med, int)
  names(out) = c("median", "lowerCI", "upperCI")
  out
}

# Function to apply the summary function by column to the list of trait Tmin, Topt, and Tmax estimates
meta.summary.TPC = function(txt) {
  out = list()
  for (i in 1:length(txt)){
    df = eval(parse(text = txt[i]))
    tmp = calcThreshPeakDists(df, Temp.xs)
    out[[i]] = apply(tmp, 2, summary.fn)
  }
  out
}

#Deterimine T-thresholds for indirect gamma(was product of two other traits, evaluate like R0)
txt = "estimated.gamma.TPdists"
i=1
meta.summary.fun = function(txt) {
  out = list()
  for (i in 1:length(txt)){
    df = eval(parse(text = txt[i]))
    out[[i]] = apply(df, 2, summary.fn)
  }
  out
}

# Create a list of traits you'd like to compare
traits = list("gamma.preds")
trait.summaries = unlist(lapply(traits, meta.summary.TPC), recursive = F)
names(trait.summaries) = traits
trait.summaries
#Export this list for quick reference later on

#sink("output_final_gamma_thresholds.txt")
#print(trait.summaries)
#sink()

##How do I find the TPC for these? Theres something wrong with this flow here
##Nas in the TPdists df.

#looking at thermal limits for estimated.lifeeggs
#trait.subset = list(c("estimated.gamma.TPdists"))
#trait.subset.summaries = lapply(trait.subset, meta.summary.fun)[[1]]
#names(trait.subset.summaries) = trait.subset[[1]]
#trait.subset.summaries
#export this as well for quick reference later
#sink("output_final_est.gamma_thresholds.txt")
#print(trait.subset.summaries)
#sink()



