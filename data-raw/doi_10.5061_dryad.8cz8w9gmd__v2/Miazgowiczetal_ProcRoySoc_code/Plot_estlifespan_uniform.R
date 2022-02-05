## Erin Mordecai, Stanford University
## August 8, 2018
## Modified by KM on April 12, 2020 to compare model fits using
##            i. the full individual dataset vs. block|temp means
##            ii. normal distribution vs. truncated normal distributions

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
library(coda) #for trait plots

# Load raw trait data - for the code below to work, trait value should be called 'trait' & temperature should be called 'T'
data.mu <- read.csv("data/mu.exp.block.values.csv")

#convert to lifepan estimates
data.mu$lifespan <- 1/data.mu$Value
data.mu$temp <- data.mu$Temp
data.mu$block <- data.mu$Block

# Load temperature sequence used to calculate trajectories
load("saved posteriors/temps.Rdata") #226 seq(0,45,0.)
N.Temp.xs <-length(Temp.xs)

############ Load traits fits
# Fits from Miazgowicz and Shapiro data to be used
load("saved posteriors/est.lifespan_quad_noT_all_uniform.Rdata")
lifespan.all.noT.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
load("saved posteriors/est.lifespan_quad_withT_all_uniform.Rdata")
lifespan.all.withT.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred


#############
####### 2. Plot thermal performance curves
#############

###### Plot TPCs with the raw data - Mordecai fits
trait.plots = function(traitname, fitname, df, labs, pch = 1){
  # Specify the temperature sequences and trait to plot
  temp = df$temp
  traitplot = df[ , which(colnames(df)==traitname)]
  load(paste("saved posteriors/", fitname, sep = ""))
  model.fit = model.out
  
  # plot the data
  plot(traitplot ~ jitter(temp, 0.5), xlim = c(0, 45), ylim = c(0, max(traitplot, na.rm = T)*1.02), 
       ylab = labs, xlab = expression(paste("Temperature (",degree,"C)")), cex.lab = 1.4, cex.axis = 1.4, pch = pch)
  
  # plot the TPC mean and 95% credible interval
  lines(model.fit$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
  lines(model.fit$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
  lines(model.fit$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
}

# Plot data means and standard errors with TPCS for visual clarity
# First a function to calculate means and standarad errors
trait.fun = function(df){
  mean = mean(df$trait, na.rm=T)
  se = std.error(df$trait, na.rm=T)
  data.frame(mean, se)
}


# Plotting function to combine data and TPCs
trait.plots.meanse = function(traitname, fitname, df, labs, add = "false", cols = "black"){
  # specify which data to plot
  dat = data.frame('temp' = df$temp, 'trait' = df[ , which(colnames(df)==traitname)])
  load(paste("saved posteriors/", fitname, sep = ""))
  model.fit = model.out
  
  # calculate mean and standard error for the trait by temperature
  data.bars = ddply(dat, .(temp), trait.fun)
  
  # set the max y value for plotting
  ymax = max(data.bars$mean + data.bars$se)
  
  # if you want to make a new plot, use add = "false"
  # if you want to add to an existing plot, use add = "true"
  # plot data means
  if (add == "false"){
    plot(mean ~ temp, data = data.bars, xlim = c(0, 45), ylim = c(0, ymax), pch = 16, 
         ylab = labs, xlab = expression(paste("Temperature (",degree,"C)")), cex.lab = 1.4, cex.axis = 1.4, col = cols)
  } else {
    points(mean ~ temp, data = data.bars, pch = 16, col = cols)
  }
  # add error bars to each data point
  for (i in 1:nrow(data.bars)) {lines(rep(data.bars$temp[i],2), c(data.bars$mean[i] - data.bars$se[i], data.bars$mean[i] + data.bars$se[i]), col = cols)}
  # add TPC mean and 95% CI
  lines(model.fit$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = cols)
  lines(model.fit$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = cols)
  lines(model.fit$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = cols)
}


##########
##put all of the fits together in the same pdf

pdf("plots/est.lifespan_fits_raw.pdf", width = 10, height = 12)
par(mfrow=c(4,3))
trait.plots("lifespan", "est.lifespan_quad_noT_all_uniform.Rdata", data.mu, "Est. Lifespan")
trait.plots("lifespan", "est.lifespan_quad_withT_all_uniform.Rdata", data.mu, "Est. Lifespan")
par(mfrow = c(1,1))
dev.off()
##########


# Plot thermal performance curves with data means and standard errors
###########
pdf("plots/est.lifespan_fits_mean_se.pdf", width = 10, height = 12)
par(mfrow = c(4,3))
trait.plots.meanse("lifespan", "est.lifespan_quad_noT_all_uniform.Rdata", data.mu, "Estimated lifespan")
trait.plots.meanse("lifespan", "est.lifespan_quad_withT_all_uniform.Rdata", data.mu, "Estimated lifespan", cols = "purple", add = "true")
legend('topleft', legend = c("all_noT", "all_withT"), lty = 1, col = c("black", "purple"), bty = 'n')
mtext("A", side = 3, at = -5, cex = 1.8, line = 2)

par(mfrow = c(1,1))
dev.off()
###########
