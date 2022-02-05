## Marta Shocket, Stanford University
## Updated by Erin Mordecai, August 6, 2018
## Updated by Kerri Miazgowicz between August 6, 2018 and August 5, 2019
## Updated by Kerri Miazgowicz on November 20,2019 to included truncated models
## Updated by KM on April 13,2020 to included revised models over data means
##
## Purpose: Calculate temperature-dependent R0 from traits fit Use Bayesian Inference (JAGS) - Traits fit with data-informed priors
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
mainDir = "C:/Users/Kerri/Desktop/Chapter1 Submission"
setwd(mainDir)

# Load libraties for fitting traits
library('R2jags') # Fits Bayesian models
library('mcmcplots') # Diagnostic plots for fits
library('MASS') # Fits distributions for informative priors
library('coda') # Calculate HPD intervals
library('dplyr') # Slice and dice data frames

############ Load traits fits
# Fits from Miazgowicz and Shapiro data - Informative priors
load("saved posteriors inf/lifespan_quadratic_inf.Rdata")
lifespan.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred

load("saved posteriors inf/bite.rate_briere_inf.Rdata")
bite.rate.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred

#Note: due to note having an appropiate trait to fit informative priors to we use the uniformative fits
load("saved posteriors/lifeeggs_briere_withT_means_uniform.Rdata")
lifetime.eggs.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred

#Note: due to the high amount of uncertainty on the inf fit; we are using the uni for this trait
load("saved posteriors/est.EFD_briere_withT_means_uniform.Rdata")
EFD.1stGC.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred

load("saved posteriors inf/est.bite_briere_inf.Rdata")
est.bite.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred

load("saved posteriors inf/mu_exponential_quadratic_inf.Rdata")
mu.exp.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred  # actually fits for 1/mu

load("saved posteriors inf/shapiro_PDR_briere_inf.Rdata")
EIP.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred  # actually fits for 1/EIP50

load("saved posteriors inf/shapiro_bc_quadratic_inf.Rdata")
bc.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred

#no informative priors can be added to this one
load("saved posteriors/gamma_quad_withT_all_uniform.Rdata")
gamma.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred

#Note: to prevent overlap in preds variable names I add a "K" after pea and MDR
load("saved posteriors inf/peaK_quadratic_inf.Rdata")
peaK.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
load("saved posteriors inf/MDRK_briere_inf.Rdata")
MDRK.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred


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

### Three R0 functions:
# 1. M depends on eggs per female per day
R0.est = function(a, bc, lf, PDR, EFD, pEA, MDR){
  M = EFD * pEA * MDR * (lf+ec.lf)^2
  (a^2 * bc * exp(-(1/(lf+ec.lf))*(1/(PDR+ec.pdr))) * M * (lf+ec.lf) )^0.5
}

# 2. M depends on lifetime fecundity
R0.full = function(a, bc, lf, y, B, pEA, MDR){
  M = B * pEA * MDR* (lf+ec.lf)
  (a^2 * bc * y * M * (lf+ec.lf) )^0.5
}

# 3. Fit from mu instead of lifespan
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

# 5. Full Johnson et al. model or Multi-species estimated model
R0.5 = R0.mu(a.preds, bcj.preds, mu.preds, PDR.preds, EFD.preds, pEA.preds, MDR.preds)
R0.5.out = calcPostQuants(R0.5, Temp.xs)
R0.5.TPdists = calcThreshPeakDists(R0.5, Temp.xs)

#6 with Paaijmans et al. 2013 MDR and pea;shapior bc and EIP50  -> MSP model lifetime or An. stephensi lifetime model
R0.6 = R0.full(bite.rate.preds, bc.preds, lifespan.preds, gamma.preds, lifetime.eggs.preds, peaK.preds, MDRK.preds)
R0.6.out = calcPostQuants(R0.6, Temp.xs)
R0.6.TPdists = calcThreshPeakDists(R0.6, Temp.xs)

#7 with Paaijmans et al. 2013 MDR and pea;shapiro bc and EIP50 -> MSP model estimated or An. stephensi estimated model
R0.7 = R0.est(est.bite.preds, bc.preds, mu.exp.preds, EIP.preds, EFD.1stGC.preds, peaK.preds, MDRK.preds)
R0.7.out = calcPostQuants(R0.7, Temp.xs)
R0.7.TPdists = calcThreshPeakDists(R0.7, Temp.xs)


##########
###### 6. Plot the resulting models
##########

# Alternative function that uses shaded polygons instead of dashed lines
plot.R0.poly = function(temp = Temp.xs, df, modname = "", add = "false", cols = "black", ltys = 1){
  if (add == "false") {
    par(mar=c(5,5,5,1))
    plot(temp, df$upperCI/max(df$upperCI, na.rm = T), lty = ltys, xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, main = modname, cex.main = 1.4, ylab = expression("Relative R"[0]), yaxt = 'n', col = "white")
    polygon(c(temp, rev(temp)), c(df$upperCI/max(df$upperCI, na.rm = T), rev(df$lowerCI/max(df$upperCI, na.rm = T))), col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  } else {
    polygon(c(temp, rev(temp)), c(df$upperCI/max(df$upperCI, na.rm = T), rev(df$lowerCI/max(df$upperCI, na.rm = T))), col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  }
  
  lines(temp, df$mean/max(df$upperCI, na.rm = T), lty = ltys, lwd = 2, col = cols)
}


# Plot all the model means and credible intervals in different panels
#####################
pdf("plots/R0_inf_vertical.pdf", width = 4, height = 12)
par(mfrow = c(3, 1))
plot.R0.poly(df = R0.5.out, modname = "Multi-species estimated")
plot.R0.poly(df = R0.6.out, modname = "An.stephensi lifetime")
plot.R0.poly(df = R0.7.out, modname = "An.stephensi estimated")
par(mfrow = c(1,1))
dev.off()
######################

# Plot all the model means together in one plot
#Plot the means and credible intervals for the main three models together in one plot
#######################
pdf("plots/R0_inf_MAIN.pdf", width = 7, height = 5)
plot.R0.poly(df = R0.5.out, add = "false", cols = "purple")
plot.R0.poly(df = R0.6.out, add = "true", cols = "black" )
plot.R0.poly(df = R0.7.out, add = "true", cols = "dodgerblue")
legend('topleft', legend = c("Multi-species estimated", expression(italic(An.~stephensi)~lifetime),expression(italic(An.~stephensi)~estimated)), lty = c(1,1,1), lwd = 2, col = c("purple","black", "dodgerblue"), bty = 'n')
dev.off()

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
    points(tmp[1,], rep((nsums+1-i), 3), pch = 20, cex = 1, col = cols)
    lines(tmp[2:3,1], rep((nsums+1-i), 2), lwd = 3, col = cols[1])
    lines(tmp[2:3,2], rep((nsums+1-i), 2), lwd = 3, col = cols[2])
    lines(tmp[2:3,3], rep((nsums+1-i), 2), lwd = 3, col = cols[3])
  }
  axis(2, at = rev(seq(1, nsums)), labels = names, tick = F)
}

# Plot the summary stats for each model
##################
pdf("plots/R0_inf_TPC_summary.pdf", width = 7 , height = 3 )
plot.summaries(R0.summaries, c("Multi-species estimated", "An.stephensi lifetime", "An.stephensi estimated"), cols = c(1,"royalblue","firebrick")) #mean, lower, upper
dev.off()
##################

#Plot the summary stats for the main three models
# apply the summary function by column to the list of R0 models
#new list of just the three I want to look at
R0s.subset = list(c("R0.5.TPdists", "R0.6.TPdists", "R0.7.TPdists"))
R0.summaries = lapply(R0s.subset, meta.summary.fun)[[1]]
names(R0.summaries) = R0s.subset[[1]]

pdf("plots2/Main_summaries_inf.pdf", width = 7 , height = 3 )
plot.summaries(R0.summaries, c("J estimated", "MSP lifetime traits", "MSP estimated traits"), cols = c(1,"royalblue","firebrick")) #mean, lower, upper
dev.off()


###########
##Modify R0 model outputs into a format that would be acceptable to give to Sadie Ryan to map over
#########
R0.5.out$temp = Temp.xs
R0.6.out$temp = Temp.xs
R0.7.out$temp = Temp.xs

#Round so that I only use up to the 0.001C spot
R0.5.out$rounded.mean <- round(R0.5.out$mean, 3)
R0.6.out$rounded.mean <- round(R0.6.out$mean, 3)
R0.7.out$rounded.mean <- round(R0.7.out$mean, 3)

R0.5.out$rounded.median <- round(R0.5.out$median, 3)
R0.6.out$rounded.median <- round(R0.6.out$median, 3)
R0.7.out$rounded.median <- round(R0.7.out$median, 3)

R0.5.out$rounded.lowerCI <- round(R0.5.out$lowerCI, 3)
R0.6.out$rounded.lowerCI <- round(R0.6.out$lowerCI, 3)
R0.7.out$rounded.lowerCI <- round(R0.7.out$lowerCI, 3)

R0.5.out$rounded.upperCI <- round(R0.5.out$upperCI, 3)
R0.6.out$rounded.upperCI <- round(R0.6.out$upperCI, 3)
R0.7.out$rounded.upperCI <- round(R0.7.out$upperCI, 3)

#write.csv(R0.6.out, "R0evaluations/04132020_MSPlifetime_R0_Sadie.csv", row.names = FALSE)
#write.csv(R0.7.out, "R0evaluations/04132020_MSPestimated_R0_Sadie.csv",row.names = FALSE)
#write.csv(R0.5.out, "R0evaluations/04132020_Jestimated_R0_Sadie.csv", row.names = FALSE)

########################Make a plot that has a histogram of the TPdists for each of the models shown for each thermal threshold
##panel A Tmin
break.n <- seq(12,28,0.4)
p1 <- hist(R0.5.TPdists$T0, breaks = break.n,freq= F)
p2 <- hist(R0.6.TPdists$T0, breaks = break.n, freq = F)
p3 <- hist(R0.7.TPdists$T0, breaks = break.n, freq = F)

#plot(p1, col =scales::alpha('seagreen3',0.5),border=F , xlim= c(10,25), ylim= c(0,1600), main = NULL)
#plot(p2, col= scales::alpha('red',0.5),border = T, add = T, density=10,angle=135)
#plot(p3, col= scales::alpha('blue',0.5),border = T, add = T, density =20, angle = 225)

#panel B Topt
break.n <- seq(22,30.8,0.4)
p4 <- hist(R0.5.TPdists$peak, breaks = break.n, freq = F)
p5 <- hist(R0.6.TPdists$peak, breaks = break.n, freq= F)
p6 <- hist(R0.7.TPdists$peak, breaks = break.n, freq = F)

#panel C Tmax
break.n <- seq(26,38,0.4)
p7 <- hist(R0.5.TPdists$Tmax, breaks = break.n, freq = F)
p8 <- hist(R0.6.TPdists$Tmax, breaks = break.n, freq= F)
p9 <- hist(R0.7.TPdists$Tmax, breaks = break.n, freq = F)
#################################################################

pdf("plots/Thresholds_density_R0s.pdf", width = 8, height = 4)
par(mfrow = c(1,3))

plot(p1, col =scales::alpha('purple',0.4),freq= F,border=scales::alpha('purple4',0.6), xlim= c(12,24),ylim = c(0,1),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, main = expression(italic(paste("T"[min]))))
plot(p2, col= scales::alpha('grey43',0.4),freq = F,border=scales::alpha('black',0.6), ylab = "Density",add = T)
plot(p3, col= scales::alpha('dodgerblue',0.4),freq= F,border=scales::alpha('royalblue',0.6), ylab = "Density",add = T)
legend('topleft', legend = c("Multi-species estimated",expression(italic(An.~stephensi)~lifetime),expression(italic(An.~stephensi)~estimated)), lty = 1, col = c("purple","black", "dodgerblue"), bty = 'n')
mtext("a", side = 3, at = 9, cex = 1.8, line = 2)

plot(p4, col =scales::alpha('purple',0.4),freq= F,border=scales::alpha('purple4',0.6), xlim= c(22,30),ylim = c(0,1.2),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, main = expression(italic(paste("T"[opt]))))
plot(p5, col= scales::alpha('grey43',0.4),freq = F,border=scales::alpha('black',0.6), ylab = "Density",add = T)
plot(p6, col= scales::alpha('dodgerblue',0.4),freq= F,border=scales::alpha('royalblue',0.6), ylab = "Density",add = T)
#legend('topleft', legend = c("Multi-species estimated",expression(italic(An.~stephensi)~lifetime),expression(italic(An.~stephensi)~estimated)), lty = 1, col = c("purple","black", "dodgerblue"), bty = 'n')
mtext("b", side = 3, at = 20, cex = 1.8, line = 2)


plot(p7, col =scales::alpha('purple',0.4),freq= F,border=scales::alpha('purple4',0.6), xlim= c(26,38),ylim = c(0,1.2),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, main = expression(italic(paste("T"[max]))))
plot(p8, col= scales::alpha('grey43',0.4),freq = F,border=scales::alpha('black',0.6), ylab = "Density",add = T)
plot(p9, col= scales::alpha('dodgerblue',0.4),freq= F,border=scales::alpha('royalblue',0.6), ylab = "Density",add = T)
#legend('topleft', legend = c("Multi-species estimated",expression(italic(An.~stephensi)~lifetime),expression(italic(An.~stephensi)~estimated)), lty = 1, col = c("purple","black", "dodgerblue"), bty = 'n')
mtext("c", side = 3, at = 22.6, cex = 1.8, line = 2)
par(mfrow = c(1,1))
dev.off()

