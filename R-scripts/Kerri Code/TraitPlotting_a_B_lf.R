## Erin Mordecai, Stanford University
## August 8, 2018
##
## Purpose: Plot trait thermal response comparisons
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) Plot thermal performance curves
##           3) Plot Tmin, Topt, and Tmax mean and 95% CI

##########
###### 1. Set up workspace, load packages, get data, etc.
##########

# Set working directory
mainDir = "C:/Users/Kerri/Desktop/Fluctuation_BayesianFits"
setwd(mainDir)

# Check whether there's a folder in the directory for saving plots
# If not, create one
subDir = "plots"
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)

# Load libraties for plotting traits
library(plyr) # Slices and dices data
library(plotrix) # For standard error function
library(coda)

# Load raw trait data - for the code below to work, trait value should be called 'trait' & temperature should be called 'T'

data.constant <- read.csv("data/constant.individual.trait.csv")
names(data.constant)[names(data.constant) == 'Treatment'] <- 'temp'
data.fluc9 <- read.csv("data/fluc9.individual.trait.csv")
names(data.fluc9)[names(data.fluc9) == 'Treatment'] <- 'temp'
data.fluc12 <- read.csv("data/fluc12.individual.trait.csv")
names(data.fluc12)[names(data.fluc12) == 'Treatment'] <- 'temp'


data.bc.EIP <- read.csv("data/forErin_ShapiroData.csv")
data.bc.EIP$inverse.EIP50 = 1/data.bc.EIP$EIP50
data.pea.MDR <- read.csv("data/Krijn_Raw_Data.csv")

# Load temperature sequence used to calculate trajectories
load("saved posteriors/temps.Rdata")
N.Temp.xs <-length(Temp.xs)

############ Load traits fits
# Fits from Miazgowicz constant dataset
load("saved posteriors/constant_lifespan_quadratic_uniform.Rdata")
lifespan.constant.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
load("saved posteriors/constant_lifetime.eggs_quadratic_uniform.Rdata")
lifetime.eggs.constant.preds <-model.out$BUGSoutput$sims.list$z.trait.mu.pred
load("saved posteriors/constant_bite.rate_briere_uniform.Rdata")
bite.rate.constant.preds <-model.out$BUGSoutput$sims.list$z.trait.mu.pred

#Fits from Miazgowicz DTR9 dataset
load("saved posteriors/dtr9_lifespan_quadratic_uniform.Rdata")
lifespan.dtr9.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
load("saved posteriors/dtr9_lifetime.eggs_quadratic_uniform.Rdata")
lifetime.eggs.dtr9.preds <-model.out$BUGSoutput$sims.list$z.trait.mu.pred
load("saved posteriors/dtr9_bite.rate_briere_uniform.Rdata")
bite.rate.dtr9.preds <-model.out$BUGSoutput$sims.list$z.trait.mu.pred

#Fits from Miazgowicz DTR12 dataset
load("saved posteriors/dtr12_lifespan_quadratic_uniform.Rdata")
lifespan.dtr12.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
load("saved posteriors/dtr12_lifetime.eggs_quadratic_uniform.Rdata")
lifetime.eggs.dtr12.preds <-model.out$BUGSoutput$sims.list$z.trait.mu.pred
load("saved posteriors/dtr12_bite.rate_briere_uniform.Rdata")
bite.rate.dtr12.preds <-model.out$BUGSoutput$sims.list$z.trait.mu.pred

#Load constant temp TPCs from Shapiro 2017 - uniform
load("saved posteriors/EIP50_briere_uniform.Rdata")
EIP.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred  # actually fits for 1/EIP50
load("saved posteriors/bc_quadratic_uniform.Rdata")
bc.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred

#Load constant temp TPCs from Shapiro 2017 - informative
load("saved posteriors/EIP50_briere_inf.Rdata")
EIP.preds.inf <- model.out$BUGSoutput$sims.list$z.trait.mu.pred  # actually fits for 1/EIP50
load("saved posteriors/bc_quadratic_inf.Rdata")
bc.preds.inf <- model.out$BUGSoutput$sims.list$z.trait.mu.pred

#Load constant temp TPCs from Paaijmans 2013 Global climate change - uniform
load("saved posteriors/pea_quadratic_uniform_Krijn2013.Rdata")
pea.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
load("saved posteriors/MDR_briere_uniform_Krijn2013.Rdata")
MDR.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred

#Load constant temp TPCs from Paaijmans 2013 Global climate change - informative
load("saved posteriors/pea_quadratic_inf.Rdata")
pea.preds.inf <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
load("saved posteriors/MDR_briere_inf.Rdata")
MDR.preds.inf <- model.out$BUGSoutput$sims.list$z.trait.mu.pred


#############
####### 2. Plot thermal performance curves
#############

###### Plot TPCs with the raw data
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

##########
pdf("plots/trait_fits_raw_uniform.pdf", width = 10, height = 15)
par(mfrow=c(5,3))
trait.plots("bite.rate", "constant_bite.rate_briere_uniform.Rdata", data.constant, "Biting rate")
trait.plots("bite.rate", "dtr9_bite.rate_briere_uniform.Rdata", data.fluc9, "Biting rate")
trait.plots("bite.rate", "dtr12_bite.rate_briere_uniform.Rdata", data.fluc12, "Biting rate")
trait.plots("lifespan", "constant_lifespan_quadratic_uniform.Rdata", data.constant, "Lifespan")
trait.plots("lifespan", "dtr9_lifespan_quadratic_uniform.Rdata", data.fluc9, "Lifespan")
trait.plots("lifespan", "dtr12_lifespan_quadratic_uniform.Rdata", data.fluc12, "Lifespan")
trait.plots("lifetime.eggs", "constant_lifetime.eggs_quadratic_uniform.Rdata", data.constant, "Lifetime egg production")
trait.plots("lifetime.eggs", "dtr9_lifetime.eggs_quadratic_uniform.Rdata", data.fluc9, "Lifetime egg production")
trait.plots("lifetime.eggs", "dtr12_lifetime.eggs_quadratic_uniform.Rdata", data.fluc12, "Lifetime egg production")
trait.plots("bc", "bc_quadratic_inf.Rdata", data.bc.EIP, "Vector competence")
trait.plots("inverse.EIP50", "EIP50_briere_inf.Rdata", data.bc.EIP, "Parasite development rate (inverse EIP50)")
trait.plots("Pea", "pea_quadratic_inf.Rdata", data.pea.MDR, "Prob.E2A survival")
trait.plots("MDR", "MDR_briere_inf.Rdata", data.pea.MDR, "Mosquito development rate")
par(mfrow = c(1,1))
dev.off()
##########

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

# Plot thermal performance curves with data means and standard errors
###########
pdf("plots/trait_fits_mean_se_uniform.pdf", width = 10, height = 16)
par(mfrow = c(5,3))
trait.plots.meanse("bite.rate", "constant_bite.rate_briere_uniform.Rdata", data.constant, "Biting rate")
mtext("A", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.meanse("bite.rate", "dtr9_bite.rate_briere_uniform.Rdata", data.fluc9, "Biting rate")
mtext("B", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.meanse("bite.rate", "dtr12_bite.rate_briere_uniform.Rdata", data.fluc12, "Biting rate")
mtext("C", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.meanse("lifespan", "constant_lifespan_quadratic_uniform.Rdata", data.constant, "Lifespan")
mtext("D", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.meanse("lifespan", "dtr9_lifespan_quadratic_uniform.Rdata", data.fluc9, "Lifespan")
mtext("E", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.meanse("lifespan", "dtr12_lifespan_quadratic_uniform.Rdata", data.fluc12, "Lifespan")
mtext("F", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.meanse("lifetime.eggs", "constant_lifetime.eggs_quadratic_uniform.Rdata", data.constant, "Lifetime egg production")
mtext("G", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.meanse("lifetime.eggs", "dtr9_lifetime.eggs_quadratic_uniform.Rdata", data.fluc9, "Lifetime egg production")
mtext("H", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.meanse("lifetime.eggs", "dtr12_lifetime.eggs_quadratic_uniform.Rdata", data.fluc12, "Lifetime egg production")
mtext("I", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.meanse("bc", "bc_quadratic_inf.Rdata", data.bc.EIP, "Vector competence")
mtext("J", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.meanse("inverse.EIP50", "EIP50_briere_inf.Rdata", data.bc.EIP, "Parasite development rate (inverse EIP50)")
mtext("K", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.meanse("Pea", "pea_quadratic_inf.Rdata", data.pea.MDR, "Prob.E2A survival")
mtext("L", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.meanse("MDR", "MDR_briere_inf.Rdata", data.pea.MDR, "Mosquito development rate")
mtext("M", side = 3, at = -5, cex = 1.8, line = 2)
par(mfrow = c(1,1))
dev.off()

##############################################
##Change this plot to include polygons as well as the dashed line
##First I need to redefine the mapping function
###########################################

# Plotting function to combine data and TPCs using shaded area instead of lines for the CI -> lwd =1.2 prior to modification
trait.plots.poly.meanse = function(traitname, fitname, df, labs, add = "false", cols = "black"){
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
    plot(mean ~ temp, data = data.bars, xlim = c(0, 45), ylim = c(0, ymax), pch = 16, col= cols, 
         ylab = labs, xlab = expression(paste("Temperature (",degree,"C)")), cex.lab = 1.4, cex.axis = 1.4)
  } else {
    points(mean ~ temp, data = data.bars, pch = 16, col = cols)
  }
  # add error bars to each data point
  for (i in 1:nrow(data.bars)) {lines(rep(data.bars$temp[i],2), c(data.bars$mean[i] - data.bars$se[i], data.bars$mean[i] + data.bars$se[i]), col = cols[1])}
  # add TPC mean and 95% CI
  polygon(c(Temp.xs, rev(Temp.xs)), c(model.fit$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(model.fit$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"])),  col = adjustcolor(cols[2], alpha.f = 0.2), border = NA)
 lines(model.fit$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 2, col = cols[1])
}

# Plot thermal performance curves with data means and standard errors
pdf("plots/trait_fits_mean_poly_se_uniform.pdf", width = 10, height = 16)
par(mfrow = c(5,3))
trait.plots.poly.meanse("bite.rate", "constant_bite.rate_briere_uniform.Rdata", data.constant, "Biting rate", cols = c("black", "dodgerblue"))
mtext("A", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly.meanse("bite.rate", "dtr9_bite.rate_briere_uniform.Rdata", data.fluc9, "Biting rate", cols = c("black", "dodgerblue"))
mtext("B", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly.meanse("bite.rate", "dtr12_bite.rate_briere_uniform.Rdata", data.fluc12, "Biting rate", cols = c("black", "dodgerblue"))
mtext("C", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly.meanse("lifespan", "constant_lifespan_quadratic_uniform.Rdata", data.constant, "Lifespan", cols = c("black", "dodgerblue"))
mtext("D", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly.meanse("lifespan", "dtr9_lifespan_quadratic_uniform.Rdata", data.fluc9, "Lifespan", cols = c("black", "dodgerblue"))
mtext("E", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly.meanse("lifespan", "dtr12_lifespan_quadratic_uniform.Rdata", data.fluc12, "Lifespan", cols = c("black", "dodgerblue"))
mtext("F", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly.meanse("lifetime.eggs", "constant_lifetime.eggs_quadratic_uniform.Rdata", data.constant, "Lifetime egg production", cols = c("black", "dodgerblue"))
mtext("G", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly.meanse("lifetime.eggs", "dtr9_lifetime.eggs_quadratic_uniform.Rdata", data.fluc9, "Lifetime egg production", cols = c("black", "dodgerblue"))
mtext("H", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly.meanse("lifetime.eggs", "dtr12_lifetime.eggs_quadratic_uniform.Rdata", data.fluc12, "Lifetime egg production", cols = c("black", "dodgerblue"))
mtext("I", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly.meanse("bc", "bc_quadratic_inf.Rdata", data.bc.EIP, "Vector competence", cols = c("black", "dodgerblue"))
mtext("J", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly.meanse("inverse.EIP50", "EIP50_briere_inf.Rdata", data.bc.EIP, "Parasite development rate (inverse EIP50)", cols = c("black", "dodgerblue"))
mtext("K", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly.meanse("Pea", "pea_quadratic_inf.Rdata", data.pea.MDR, "Prob.E2A survival", cols = c("black", "dodgerblue"))
mtext("L", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly.meanse("MDR", "MDR_briere_inf.Rdata", data.pea.MDR, "Mosquito development rate", cols = c("black", "dodgerblue"))
mtext("M", side = 3, at = -5, cex = 1.8, line = 2)
par(mfrow = c(1,1))
dev.off()



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


#######
#To be able to properly compare lifetime egg production verses estimated lifetime egg production (which is EFD (1st GC)* inv(Mu.exp method))
# I need to create a function that generates this new curve and the associated CI - similar to how it is done with R0


# Alternative plotting function that uses shaded polygons instead of dashed lines --> not relative as it is in R0
plot.trait.poly = function(temp = Temp.xs, df, modname = "", add = "false", cols = "black"){
  if (add == "false") {
    par(mar=c(5,5,5,1))
    plot(temp, df$upperCI, lty = 1, xlab = expression(paste("Temperature (",degree,"C)")), ylab = "test", cex.axis = 1.2, cex.lab = 1.2, main = modname, cex.main = 1.5, yaxt = 'n', col = "white")
    polygon(c(temp, rev(temp)), c(df$upperCI, rev(df$lowerCI)), col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  } else {
    polygon(c(temp, rev(temp)), c(df$upperCI, rev(df$lowerCI)), col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  }
  
  lines(temp, df$mean, lty = 1, lwd = 2, col = cols)
}

# Alternative plotting function to display CI as dashed lines
plot.trait2 = function(temp = Temp.xs, df, modname = "", add = "false", cols = 1){
  if (add == "false") {
    par(mar=c(5,5,5,1))
    plot(temp, df$upperCI, type = "l", lty = 2, lwd = 2, xlab = expression(paste("Temperature (",degree,"C)")), ylab = modename, cex.axis = 1.2, cex.lab = 1.2, main = modname, cex.main = 1.5, yaxt = 'n', col = cols)
  } else {
    lines(temp, df$upperCI, lty = 2, lwd = 2, col = cols)
  }
  lines(temp, df$mean, lty = 1, lwd = 2, col = cols)
  lines(temp, df$lowerCI, lty = 2, lwd = 2, col = cols)
}


###########
pdf("plots/traits_temperature_overlay_uniform.pdf", width = 8, height = 12)
par(mfrow = c(3,3))
# Biting rate
trait.plots.poly.meanse("bite.rate", "constant_bite.rate_briere_uniform.Rdata", data.constant, "Biting rate")
trait.plots.poly.meanse("bite.rate", "dtr9_bite.rate_briere_uniform.Rdata", data.fluc9, "Biting rate", add = "true", cols = "royalblue")
trait.plots.poly.meanse("bite.rate", "dtr12_bite.rate_briere_uniform.Rdata", data.fluc12, "Biting rate", add = "true", cols = "purple")
legend('topleft', legend = c("dtr0", "dtr9","dtr12"), lty = 1, col = c("black", "royalblue","purple"), bty = 'n')
mtext("A", side = 3, at = -5, cex = 1.8, line = 2)
# Lifespan
trait.plots.poly.meanse("lifespan", "constant_lifespan_quadratic_uniform.Rdata", data.constant, "Lifespan")
trait.plots.poly.meanse("lifespan", "dtr9_lifespan_quadratic_uniform.Rdata", data.fluc9, "Lifespan", add = "true", cols = "royalblue")
trait.plots.poly.meanse("lifespan", "dtr12_lifespan_quadratic_uniform.Rdata", data.fluc12, "Lifespan", add = "true", cols = "purple")
legend('topleft', legend = c("dtr0", "dtr9","dtr12"), lty = 1, col = c("black", "royalblue","purple"), bty = 'n')
mtext("B", side = 3, at = -5, cex = 1.8, line = 2)
# Lifetime egg production
trait.plots.poly.meanse("lifetime.eggs", "constant_lifetime.eggs_quadratic_uniform.Rdata", data.constant, "Lifetime egg production")
trait.plots.poly.meanse("lifetime.eggs", "dtr9_lifetime.eggs_quadratic_uniform.Rdata", data.fluc9, "Lifetime egg production", add = "true", cols = "royalblue")
trait.plots.poly.meanse("lifetime.eggs", "dtr12_lifetime.eggs_quadratic_uniform.Rdata", data.fluc12, "Lifetime egg production", add = "true", cols = "purple")
legend('topleft', legend = c("dtr0", "dtr9","dtr12"), lty = 1, col = c("black", "royalblue","purple"), bty = 'n')
mtext("C", side = 3, at = -5, cex = 1.8, line = 2)
#Remaining Traits
trait.plots.poly.meanse("bc", "bc_quadratic_inf.Rdata", data.bc.EIP, "Vector competence")
mtext("D", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly.meanse("inverse.EIP50", "EIP50_briere_inf.Rdata", data.bc.EIP, "Parasite development rate")
mtext("E", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly.meanse("Pea", "pea_quadratic_inf.Rdata", data.pea.MDR, "Prob.E2A survival")
mtext("F", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly.meanse("MDR", "MDR_briere_inf.Rdata", data.pea.MDR, "Mosquito development rate")
mtext("G", side = 3, at = -5, cex = 1.8, line = 2)

par(mfrow = c(1,1))
dev.off()
############


#Compare just the three traits in this experiment bite rate, lifespan, and lifetime egg production using polys. -> No dashed lines at border of CI
pdf("plots/traits_a_lf_B_mean_poly_se_uniform.pdf", width = 6, height = 10)
par(mfrow = c(3,1))
trait.plots.poly.meanse("bite.rate", "constant_bite.rate_briere_uniform.Rdata", data.constant, "Biting rate", cols = c("black", "black"))
trait.plots.poly.meanse("bite.rate", "dtr9_bite.rate_briere_uniform.Rdata", data.fluc9, "Biting rate", cols = c("dodgerblue", "dodgerblue"), add = "true")
trait.plots.poly.meanse("bite.rate", "dtr12_bite.rate_briere_uniform.Rdata", data.fluc12, "Biting rate", cols = c("purple", "purple"), add = "true")
legend('topleft', legend = c("dtr0", "dtr9","dtr12"), lty = 1, col = c("black", "dodgerblue","purple"), bty = 'n')
mtext("A", side = 3, at = -5, cex = 1.8, line = 2)

trait.plots.poly.meanse("lifespan", "constant_lifespan_quadratic_uniform.Rdata", data.constant, "Lifespan", cols = c("black", "black"))
trait.plots.poly.meanse("lifespan", "dtr9_lifespan_quadratic_uniform.Rdata", data.fluc9, "Lifespan", cols = c("dodgerblue", "dodgerblue"), add = "true")
trait.plots.poly.meanse("lifespan", "dtr12_lifespan_quadratic_uniform.Rdata", data.fluc12, "Lifespan", cols = c("purple", "purple"), add = "true")
legend('topleft', legend = c("dtr0", "dtr9","dtr12"), lty = 1, col = c("black", "dodgerblue","purple"), bty = 'n')
mtext("B", side = 3, at = -5, cex = 1.8, line = 2)

trait.plots.poly.meanse("lifetime.eggs", "constant_lifetime.eggs_quadratic_uniform.Rdata", data.constant, "Lifetime egg production", cols = c("black", "black"))
trait.plots.poly.meanse("lifetime.eggs", "dtr9_lifetime.eggs_quadratic_uniform.Rdata", data.fluc9, "Lifetime egg production", cols = c("dodgerblue", "dodgerblue"), add = "true")
trait.plots.poly.meanse("lifetime.eggs", "dtr12_lifetime.eggs_quadratic_uniform.Rdata", data.fluc12, "Lifetime egg production", cols = c("purple", "purple"), add = "true")
legend('topleft', legend = c("dtr0", "dtr9","dtr12"), lty = 1, col = c("black", "dodgerblue","purple"), bty = 'n')
mtext("C", side = 3, at = -5, cex = 1.8, line = 2)
par(mfrow = c(1,1))
dev.off()

###Put the constant temperature plots into a single figure for the supplemental
pdf("plots/traits_constant_sources.pdf", width = 8, height = 8)
par(mfrow = c(2,2))
trait.plots.poly.meanse("bc", "bc_quadratic_inf.Rdata", data.bc.EIP, "Vector competence",cols = c("black", "black"))
mtext("A", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly.meanse("inverse.EIP50", "EIP50_briere_inf.Rdata", data.bc.EIP, "Parasite development rate",cols = c("black", "black"))
mtext("B", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly.meanse("Pea", "pea_quadratic_inf.Rdata", data.pea.MDR, "Prob.E2A survival",cols = c("black", "black"))
mtext("C", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly.meanse("MDR", "MDR_briere_inf.Rdata", data.pea.MDR, "Mosquito development rate",cols = c("black", "black"))
mtext("D", side = 3, at = -5, cex = 1.8, line = 2)
par(mfrow = c(1,1))
dev.off()

############
####### 3. Plot Tmin, Topt, and Tmax mean and 95% CI
###########

###### Function to calculate distributions of T0, Tmax, and peak for each TPC

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
# Create a list of traits you'd like to compare
traits = list("bite.rate.constant.preds", "bite.rate.dtr9.preds", "bite.rate.dtr12.preds","lifespan.constant.preds","lifespan.dtr9.preds","lifespan.dtr12.preds","lifetime.eggs.constant.preds","lifetime.eggs.dtr9.preds","lifetime.eggs.dtr12.preds", "bc.preds.inf", "EIP.preds.inf", "pea.preds.inf", "MDR.preds.inf")
trait.summaries = unlist(lapply(traits, meta.summary.TPC), recursive = F)
names(trait.summaries) = traits
##################
pdf("plots/trait_summaries_uniform.pdf", width = 8, height = 6)
plot.summaries(trait.summaries, c("biting rate (dtr0)", "biting rate (dtr9)","biting rate (dtr12)","lifespan (dtr0)", "lifespan (dtr9)", "lifespan (dtr12)", "lifetime egg production (dtr0)","lifetime egg production (dtr9)", "lifetime egg production (dtr12)", "vector competence", "pathogen development rate", "Prob. egg to adult survival", "Mosqutio development rate"), xlims = c(0, 45), cols = c(1, "royalblue", "firebrick"))
dev.off()

#Make a summary plot of the traits Miazgowicz measured
traits = list("bite.rate.constant.preds", "bite.rate.dtr9.preds", "bite.rate.dtr12.preds","lifespan.constant.preds","lifespan.dtr9.preds","lifespan.dtr12.preds","lifetime.eggs.constant.preds","lifetime.eggs.dtr9.preds","lifetime.eggs.dtr12.preds")
trait.summaries = unlist(lapply(traits, meta.summary.TPC), recursive = F)
names(trait.summaries) = traits

pdf("plots/trait_a_lf_B_summaries_uniform.pdf", width = 8, height = 6)
plot.summaries(trait.summaries, c("biting rate (dtr0)", "biting rate (dtr9)","biting rate (dtr12)","lifespan (dtr0)", "lifespan (dtr9)", "lifespan (dtr12)", "lifetime egg production (dtr0)","lifetime egg production (dtr9)", "lifetime egg production (dtr12)"), xlims = c(0, 45), cols = c(1, "royalblue", "firebrick"))
dev.off()
##################

################
################
################ Make figures which include rate summation estimates on them 
###############
###############

########Formalized approach.  -> in RateSummation.R code
# (1) Evaluate trait median curve values at 0.1C (what the temp dataframe is (Note: I can do at 0.01C intervals but I'll have to refit the Bayesian models))
# (2) Import the dataframe of temperature programs
# (3)a For each temperature (seq(0,45, by 0.1C) -> Create a function which rounds the the temperature program df to the nearest 0.1C,
# (3)b Creates a new dataframe which lists trait performance at each time interval
# (3)c Computes the average of the hourly estimates to generate a single value
# (3)d Places these values into a dataframe which consists of 'temp' & 'traitvalue'
########Bayesian fits on trait estimates -> in TraitFittingEstimates.R code

# Plotting function to add on estimated performance lines; doesn't use poly for CI of Bayesian fits or dashed lines
trait.estimates.plots.lines = function(traitname, fitname, df, labs, add = "false", cols = "black"){
  # specify which data to plot
  dat = data.frame('temp' = df$temp, 'trait' = df[ , which(colnames(df)==traitname)])
  load(paste("saved posteriors/", fitname, sep = ""))
  model.fit = model.out
  
 # add TPC mean and 95% CI
  #polygon(c(Temp.xs, rev(Temp.xs)), c(model.fit$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(model.fit$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"])),  col = adjustcolor(cols[2], alpha.f = 0.2), border = NA)
  lines(model.fit$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty= 2, lwd = 2, col = cols[1])
  #lines(model.fit$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = cols[2])
  #lines(model.fit$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = cols[2])
}


###Note the below section of code will not run properly until after RateSummation.R ,and  TraitFittingEstimates.R have been run once

#Load fluctuation temp TPCs from TraitFittingEstimates.R
load("saved posteriors/dtr9_estimated_bite.rate_briere_uniform.Rdata")
est.bite.rate.dtr9.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
load("saved posteriors/dtr12_estimated_bite.rate_briere_uniform.Rdata")
est.bite.rate.dtr12.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
load( "saved posteriors/dtr9_estimated_lifespan_quadratic_uniform.Rdata")
est.lifespan.dtr9.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
load("saved posteriors/dtr12_estimated_lifespan_quadratic_uniform.Rdata")
est.lifespan.dtr12.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
load("saved posteriors/dtr9_estimated_lifetime.eggs_quadratic_uniform.Rdata")
est.lifetime.eggs.dtr9.preds <-model.out$BUGSoutput$sims.list$z.trait.mu.pred 
load("saved posteriors/dtr12_estimated_lifetime.eggs_quadratic_uniform.Rdata")
est.lifetime.eggs.dtr12.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred

# Plot the summary stats for each model
# Create a list of traits you'd like to compare
traits = list("est.bite.rate.dtr9.preds", "est.bite.rate.dtr12.preds", "est.lifespan.dtr9.preds","est.lifespan.dtr12.preds","est.lifetime.eggs.dtr9.preds","est.lifetime.eggs.dtr12.preds")
trait.summaries = unlist(lapply(traits, meta.summary.TPC), recursive = F)
names(trait.summaries) = traits


#Plot a full figure with the observed trait values
#Compare just the three traits in this experiment bite rate, lifespan, and lifetime egg production using polys. -> No dashed lines at border of CI
pdf("plots/traits_a_lf_B_poly_withEstimates.pdf", width = 5, height = 10)
par(mfrow = c(3,1))
trait.plots.poly.meanse("bite.rate", "constant_bite.rate_briere_uniform.Rdata", data.constant, "Biting rate", cols = c("black", "black"))
trait.plots.poly.meanse("bite.rate", "dtr9_bite.rate_briere_uniform.Rdata", data.fluc9, "Biting rate", cols = c("dodgerblue", "dodgerblue"), add = "true")
trait.plots.poly.meanse("bite.rate", "dtr12_bite.rate_briere_uniform.Rdata", data.fluc12, "Biting rate", cols = c("purple", "purple"), add = "true")
trait.estimates.plots.lines("bite.rate","dtr9_estimated_bite.rate_briere_uniform.Rdata" ,data.constant, "Biting rate", cols = c("dodgerblue", "dodgerblue"), add = "true")
trait.estimates.plots.lines("bite.rate","dtr12_estimated_bite.rate_briere_uniform.Rdata" ,data.constant, "Biting rate", cols = c("purple", "purple"), add = "true")
legend('topleft', legend = c("dtr0", "dtr9","dtr12",NA, "observed", "estimated"),lwd = 2, lty = c(1,1,1,NA, 1,2), col = c("black", "dodgerblue","purple", NA, "black", "black"), bty = 'n')
mtext("A", side = 3, at = -5, cex = 1.8, line = 2)

trait.plots.poly.meanse("lifespan", "constant_lifespan_quadratic_uniform.Rdata", data.constant, "Lifespan", cols = c("black", "black"))
trait.plots.poly.meanse("lifespan", "dtr9_lifespan_quadratic_uniform.Rdata", data.fluc9, "Lifespan", cols = c("dodgerblue", "dodgerblue"), add = "true")
trait.plots.poly.meanse("lifespan", "dtr12_lifespan_quadratic_uniform.Rdata", data.fluc12, "Lifespan", cols = c("purple", "purple"), add = "true")
trait.estimates.plots.lines("lifespan", "dtr9_estimated_lifespan_quadratic_uniform.Rdata", data.constant, "Lifespan", cols = c("dodgerblue", "dodgerblue"), add = "true")
trait.estimates.plots.lines("lifespan", "dtr12_estimated_lifespan_quadratic_uniform.Rdata", data.constant, "Lifespan", cols = c("purple", "purple"), add = "true")
legend('topleft', legend = c("dtr0", "dtr9","dtr12",NA, "observed", "estimated"),lwd = 2, lty = c(1,1,1,NA, 1,2), col = c("black", "dodgerblue","purple", NA, "black", "black"), bty = 'n')
mtext("B", side = 3, at = -5, cex = 1.8, line = 2)

trait.plots.poly.meanse("lifetime.eggs", "constant_lifetime.eggs_quadratic_uniform.Rdata", data.constant, "Lifetime egg production", cols = c("black", "black"))
trait.plots.poly.meanse("lifetime.eggs", "dtr9_lifetime.eggs_quadratic_uniform.Rdata", data.fluc9, "Lifetime egg production", cols = c("dodgerblue", "dodgerblue"), add = "true")
trait.plots.poly.meanse("lifetime.eggs", "dtr12_lifetime.eggs_quadratic_uniform.Rdata", data.fluc12, "Lifetime egg production", cols = c("purple", "purple"), add = "true")
trait.estimates.plots.lines("lifetime.eggs", "dtr9_estimated_lifetime.eggs_quadratic_uniform.Rdata", data.constant, "Lifetime egg production", cols = c("dodgerblue", "dodgerblue"), add = "true")
trait.estimates.plots.lines("lifetime.eggs", "dtr12_estimated_lifetime.eggs_quadratic_uniform.Rdata", data.constant, "Lifetime egg production", cols = c("purple", "purple"), add = "true")
legend('topleft', legend = c("dtr0", "dtr9","dtr12",NA, "observed", "estimated"),lwd = 2, lty = c(1,1,1,NA, 1,2), col = c("black", "dodgerblue","purple", NA, "black", "black"), bty = 'n')
mtext("C", side = 3, at = -5, cex = 1.8, line = 2)
par(mfrow = c(1,1))
dev.off()

#############Plot for presentation to be in the same format as the source 
###Put the constant temperature plots into a single figure for the supplemental
pdf("plots/traits_presentation_model_summary.pdf", width = 8, height = 8)
par(mfrow = c(2,2))
trait.plots.poly.meanse("bite.rate", "constant_bite.rate_briere_uniform.Rdata", data.constant, "Biting rate", cols = c("black", "black"))
trait.plots.poly.meanse("bite.rate", "dtr9_bite.rate_briere_uniform.Rdata", data.fluc9, "Biting rate", cols = c("dodgerblue", "dodgerblue"), add = "true")
trait.plots.poly.meanse("bite.rate", "dtr12_bite.rate_briere_uniform.Rdata", data.fluc12, "Biting rate", cols = c("purple", "purple"), add = "true")
trait.estimates.plots.lines("bite.rate","dtr9_estimated_bite.rate_briere_uniform.Rdata" ,data.constant, "Biting rate", cols = c("dodgerblue", "dodgerblue"), add = "true")
trait.estimates.plots.lines("bite.rate","dtr12_estimated_bite.rate_briere_uniform.Rdata" ,data.constant, "Biting rate", cols = c("purple", "purple"), add = "true")
#legend('topleft', legend = c("dtr0", "dtr9","dtr12"),lwd = 2, lty = c(1,1,1), col = c("black", "dodgerblue","purple"), bty = 'n')
mtext("A", side = 3, at = -5, cex = 1.8, line = 2)

trait.plots.poly.meanse("lifespan", "constant_lifespan_quadratic_uniform.Rdata", data.constant, "Lifespan", cols = c("black", "black"))
trait.plots.poly.meanse("lifespan", "dtr9_lifespan_quadratic_uniform.Rdata", data.fluc9, "Lifespan", cols = c("dodgerblue", "dodgerblue"), add = "true")
trait.plots.poly.meanse("lifespan", "dtr12_lifespan_quadratic_uniform.Rdata", data.fluc12, "Lifespan", cols = c("purple", "purple"), add = "true")
trait.estimates.plots.lines("lifespan", "dtr9_estimated_lifespan_quadratic_uniform.Rdata", data.constant, "Lifespan", cols = c("dodgerblue", "dodgerblue"), add = "true")
trait.estimates.plots.lines("lifespan", "dtr12_estimated_lifespan_quadratic_uniform.Rdata", data.constant, "Lifespan", cols = c("purple", "purple"), add = "true")
#legend('topleft', legend = c("dtr0", "dtr9","dtr12",NA, "observed", "estimated"),lwd = 2, lty = c(1,1,1,NA, 1,2), col = c("black", "dodgerblue","purple", NA, "black", "black"), bty = 'n')
mtext("B", side = 3, at = -5, cex = 1.8, line = 2)

trait.plots.poly.meanse("lifetime.eggs", "constant_lifetime.eggs_quadratic_uniform.Rdata", data.constant, "Lifetime egg production", cols = c("black", "black"))
trait.plots.poly.meanse("lifetime.eggs", "dtr9_lifetime.eggs_quadratic_uniform.Rdata", data.fluc9, "Lifetime egg production", cols = c("dodgerblue", "dodgerblue"), add = "true")
trait.plots.poly.meanse("lifetime.eggs", "dtr12_lifetime.eggs_quadratic_uniform.Rdata", data.fluc12, "Lifetime egg production", cols = c("purple", "purple"), add = "true")
trait.estimates.plots.lines("lifetime.eggs", "dtr9_estimated_lifetime.eggs_quadratic_uniform.Rdata", data.constant, "Lifetime egg production", cols = c("dodgerblue", "dodgerblue"), add = "true")
trait.estimates.plots.lines("lifetime.eggs", "dtr12_estimated_lifetime.eggs_quadratic_uniform.Rdata", data.constant, "Lifetime egg production", cols = c("purple", "purple"), add = "true")
#legend('topleft', legend = c("dtr0", "dtr9","dtr12",NA, "observed", "estimated"),lwd = 2, lty = c(1,1,1,NA, 1,2), col = c("black", "dodgerblue","purple", NA, "black", "black"), bty = 'n')
mtext("C", side = 3, at = -5, cex = 1.8, line = 2)
par(mfrow = c(1,1))
dev.off() 

