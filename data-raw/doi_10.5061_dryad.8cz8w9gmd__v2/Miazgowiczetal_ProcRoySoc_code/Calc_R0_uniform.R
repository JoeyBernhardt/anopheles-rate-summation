## Marta Shocket, Stanford University
## Updated by Erin Mordecai, August 6, 2018
## Updated by Kerri Miazgowicz between August 6, 2018 and August 5, 2019
## Updated by Kerri Miazgowicz on November 20,2019 to include truncated models
## Updated by KM on April 13,2020 to include the updated models over the data means
## Updated by KM on April 14,202 to report the DIC value of each model
##
## Purpose: Calculate temperature-dependent R0 from traits fit Use Bayesian Inference (JAGS)
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) Function to calculate quantiles of derived quantity posteriors (for plotting)
##           3) Functions defining R0 equations
##           4) Calculate R0 posteriors and quantiles
##           5) Function to calculate distributions of T0, Tmax, and peak R0
##           6) Calculate distributions of T0, Tmax, and peak R0


##########
###### 1. Set up workspace, load packages, get data, etc.
##########

# Set working directory
mainDir = "C:/Users/Kerri/Desktop/Chapter 1 Submission"
setwd(mainDir)

# Load libraties for fitting traits
library('R2jags') # Fits Bayesian models
library('mcmcplots') # Diagnostic plots for fits
library('MASS') # Fits distributions for informative priors
library('coda') # Calculate HPD intervals
library('dplyr') # Slice and dice data frames

##########Load uniform truncated normal trait fits for the functional forms not used in R0 to determine the DIC values
load("saved posteriors/gamma_briere_withT_all_uniform.Rdata")
gamma.briere.DIC <- model.out$BUGSoutput$DIC   #block-specific gompertz survivorship (Miazgowicz)with Shapiro EIP50
load("saved posteriors/krijn_pea_briere_withT_uniform.Rdata")
peaK.briere.DIC <- model.out$BUGSoutput$DIC
load("saved posteriors/krijn_MDR_quad_withT_uniform.Rdata")
MDRK.quad.DIC <- model.out$BUGSoutput$DIC
load("saved posteriors/lifeeggs_quad_withT_means_uniform.Rdata")
lifetime.eggs.quad.DIC <-  model.out$BUGSoutput$DIC
load("saved posteriors/shapiro_PDR_quad_withT_uniform.Rdata")
EIP.quad.DIC<- model.out$BUGSoutput$DIC
load("saved posteriors/shapiro_bc_briere_withT_uniform.Rdata")
bc..briere.DIC <- model.out$BUGSoutput$DIC





############ Load uniform truncated normal traits fits
# Fits from Miazgowicz and Shapiro data
load("saved posteriors/lifespan_quad_withT_means_uniform.Rdata")
lifespan.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
lifespan.DIC <- model.out$BUGSoutput$DIC


load("saved posteriors/bite.rate_briere_withT_means_uniform.Rdata")
bite.rate.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
bite.rate.DIC <-  model.out$BUGSoutput$DIC

load("saved posteriors/lifeeggs_briere_withT_means_uniform.Rdata")
lifetime.eggs.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
lifetime.eggs.DIC <-  model.out$BUGSoutput$DIC


load("saved posteriors/est.EFD_briere_withT_means_uniform.Rdata")
EFD.1stGC.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
EFD.1stGC.DIC <- model.out$BUGSoutput$DIC

load("saved posteriors/est.bite_briere_withT_means_uniform.Rdata")
est.bite.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
est.bite.DIC <- model.out$BUGSoutput$DIC

load("saved posteriors/est.lifespan_quad_withT_all_uniform.Rdata")
mu.exp.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred  # actually fits for 1/mu
mu.exp.DIC <- model.out$BUGSoutput$DIC

load("saved posteriors/shapiro_PDR_briere_withT_uniform.Rdata")
EIP.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred  # actually fits for 1/EIP50
EIP.DIC<- model.out$BUGSoutput$DIC

load("saved posteriors/shapiro_bc_quad_withT_uniform.Rdata")
bc.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
bc.DIC <- model.out$BUGSoutput$DIC

#Note: to prevent overlap in preds variable names I add a "K" after pea and MDR
load("saved posteriors/krijn_pea_quad_withT_uniform.Rdata")
peaK.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
peaK.DIC <- model.out$BUGSoutput$DIC

load("saved posteriors/krijn_MDR_briere_withT_uniform.Rdata")
MDRK.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
MDRK.DIC <- model.out$BUGSoutput$DIC

#Gamma substitutions
load("saved posteriors/gamma_quad_withT_all_uniform.Rdata")
gamma.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred    #block-specific gompertz survivorship (Miazgowicz)with Shapiro EIP50

# Fits from Johnson et al 2015
fits = paste("Johnson et al posteriors/", c(
  "a.preds.Rdata",
  "PDR.preds.Rdata",
  "MDR.preds.Rdata",
  "EFD.preds.Rdata",
  "pEA.preds.Rdata",
  "bcj.preds.Rdata",
  "mu.preds.Rdata"), sep = "")
for (i in 1:length (fits)) load(fits[i])

# Temperature sequences that was used to make all calculations
load("saved posteriors/temps.Rdata")
N.Temp.xs <-length(Temp.xs)

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
###### 3. Functions defining R0 equations
##########

# Arguments: dataframes with posterior distributions of each trait predicted over a gradient (i.e., 'trait.preds')

# constant to keep lifespan from being numerically zero
# assume minimum survival time is half an hour
ec.lf = 1/48

# constant to keep PDR from being numerically zero
# assume the maximum EIP is 100 days
ec.pdr = 1/100

###  Three R0 functions:
# 1a. M depends on eggs per female per day;
R0.est = function(a, bc, lf, PDR, EFD, pEA, MDR){
  M = EFD * pEA * MDR * (lf+ec.lf)^2
  (a^2 * bc * exp(-(1/(lf+ec.lf))*(1/(PDR+ec.pdr))) * M * (lf+ec.lf) )^0.5
}

# 2. M depends on lifetime fecundity; gamma TPC is substituted for exp^(-1/(lf*PDR))expression
R0.full = function(a, bc, lf, y, B, pEA, MDR){
  M = B * pEA * MDR * (lf+ec.lf)
  (a^2 * bc * y * M * (lf+ec.lf) )^0.5
}

# 3. Fit from mu instead of lifespan (Used for the Johnson Model)
R0.mu = function(a, bc, mu, PDR, EFD, pEA, MDR){
  lf = 1/mu
  M = EFD * pEA * MDR * (lf+ec.lf)^2
  (a^2 * bc * exp(-(1/(lf+ec.lf))*(1/(PDR+ec.pdr))) * M * (lf+ec.lf) )^0.5
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
    index.list <- which(input[i, ] > 0) # Create vector of list of indices where relative R0 > 0
    length.index.list <- length(index.list)
    output.df$T0[i] <- temp.list[index.list[1] - 1] # Store T0 (index prior to first value in index.list)
    output.df$Tmax[i] <- temp.list[index.list[length.index.list] + 1] # Store Tmax (index after to last value in index.list)
  }
  
  output.df # return
  
}

##########
###### 5. Calculate R0 posteriors and quantiles
##########

### Several R0 calculations to make and compare
# Calculate full posteriors and quantiles for each model formulation

# 5. Full Johnson et al. model
R0.5 = R0.mu(a.preds, bcj.preds, mu.preds, PDR.preds, EFD.preds, pEA.preds, MDR.preds)
R0.5.uni.out = calcPostQuants(R0.5, Temp.xs)
R0.5.TPdists = calcThreshPeakDists(R0.5, Temp.xs)

#6 An. lifetime model: Miazgowicz gomp survival; lifetime a,B,lf; Paaijmans et al. 2013 MDR and pea; Shapiro 2017 bc and EIP50
R0.6 = R0.full(bite.rate.preds, bc.preds, lifespan.preds, gamma.preds, lifetime.eggs.preds, peaK.preds, MDRK.preds)
R0.6.uni.out = calcPostQuants(R0.6, Temp.xs)
R0.6.TPdists = calcThreshPeakDists(R0.6, Temp.xs)

#7 An. estimated model: Miazgowicz exp survival; estimated a*,EFD*, lf*; Paaijmans et al. 2013 MDR and pea
R0.7 = R0.est(est.bite.preds, bc.preds, mu.exp.preds, EIP.preds, EFD.1stGC.preds, peaK.preds, MDRK.preds)
R0.7.uni.out = calcPostQuants(R0.7, Temp.xs)
R0.7.TPdists = calcThreshPeakDists(R0.7, Temp.xs)

##########
###### 6. Plot the resulting models
##########

# Alternative function that uses shaded polygons instead of dashed lines
plot.R0.poly = function(temp = Temp.xs, df, modname = "", add = "false", cols = "black"){
  if (add == "false") {
    par(mar=c(5,5,5,1))
    plot(temp, df$upperCI/max(df$upperCI, na.rm = T), lty = 1, xlab = expression(paste("Temperature (",degree,"C)")), ylab = expression(R[0]), cex.axis = 1.2, cex.lab = 1.2, main = modname, cex.main = 1.5, yaxt = 'n', col = "white")
    polygon(c(temp, rev(temp)), c(df$upperCI/max(df$upperCI, na.rm = T), rev(df$lowerCI/max(df$upperCI, na.rm = T))), col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  } else {
    polygon(c(temp, rev(temp)), c(df$upperCI/max(df$upperCI, na.rm = T), rev(df$lowerCI/max(df$upperCI, na.rm = T))), col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  }
  
  lines(temp, df$mean/max(df$upperCI, na.rm = T), lty = 1, lwd = 2, col = cols)
}

# Plot all the model means and credible intervals in different panels
#####################
pdf("plots/R0_all_uniform_vert.pdf", width = 4, height = 20)
par(mfrow = c(9, 1))
plot.R0.poly(df = R0.5.uni.out, modname = "Johnson\n estimated traits", cols = "purple")
plot.R0.poly(df = R0.6.uni.out, modname = "An. stephensi \n lifetime traits", cols = "black")
plot.R0.poly(df = R0.7.uni.out, modname = "An. stephensi \n estimated traits", cols = "dodgerblue")
par(mfrow = c(1,1))
dev.off()
######################

###Plot only the three main together
pdf("plots/R0_together_uniform.pdf", width = 7, height = 5)
plot.R0.poly(df = R0.5.uni.out, modname = "Johnson\n estimated traits", cols = "purple")
plot.R0.poly(df = R0.6.uni.out, modname = "An. stephensi \n lifetime traits", cols = "black",add = "true")
plot.R0.poly(df = R0.7.uni.out, modname = "An. stephensi \n estimated traits", cols = "dodgerblue", add = "true")
legend('topleft', legend = c("Multi-species estimated", "An.stephensi lifetime", "An.stephensi estimated"), lty = 1, col = c("purple","black","dodgerblue"), bty = 'n')
dev.off()

#####################

############
####### 7. Compare peak, T0, and Tmax across models
###########
summary.fn = function(df){
  med = median(df)
  int = c(HPDinterval(mcmc(df)))
  out = c(med, int)
  names(out) = c("median", "lowerCI", "upperCI")
  out
}

# apply the summary function by column to the list of R0 models
R0s = list(paste("R0.", seq(5,7,1), ".TPdists", sep = ""))
meta.summary.fun = function(txt) {
  out = list()
  for (i in 1:length(txt)){
  df = eval(parse(text = txt[i]))
  out[[i]] = apply(df, 2, summary.fn)
  }
  out
}

R0.summaries = lapply(R0s, meta.summary.fun)[[1]]
names(R0.summaries) = R0s[[1]]

# Create a plotting function with these summaries as inputs

plot.summaries = function(sums, names, xlims = c(10, 40), cols = c(1,1,1)){
  par(las = 1, mar = c(5, 10, 1, 1))
  nsums = length(sums)
  plot(0, -1, xlim = xlims, ylim = c(0, nsums+1), xlab = expression(paste("Temperature (",degree,"C)")), yaxt = 'n', ylab = "", cex = 2, bty = 'n')
  for (i in 1:nsums){
    tmp = sums[[i]]
    points(tmp[1,], rep((nsums+1-i), 3), pch = 16, cex = 1.8, col = cols)
    lines(tmp[2:3,1], rep((nsums+1-i), 2), lwd = 3, col = cols[1])
    lines(tmp[2:3,2], rep((nsums+1-i), 2), lwd = 3, col = cols[2])
    lines(tmp[2:3,3], rep((nsums+1-i), 2), lwd = 3, col = cols[3])
  }
  axis(2, at = rev(seq(1, nsums)), labels = names, tick = F)
}

# Plot the summary stats for each model
##################
pdf("plots/R0_uni_TPC_summary.pdf")
plot.summaries(R0.summaries, c("Multi-species estimated","An.stephensi lifetime","An.stephensi estimated"))
dev.off()
##################

