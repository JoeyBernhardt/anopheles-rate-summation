## Erin Mordecai, Stanford University
## August 8, 2018
## Updated by Kerri Miazgowicz on April 12,2020
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
mainDir = "C:/Users/Kerri/Desktop/Chapter1 Submission"
setwd(mainDir)

# Check whether there's a folder in the directory for saving plots
# If not, create one
subDir = "plots"
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)

# Load libraties for plotting traits
library(plyr) # Slices and dices data
library(plotrix) # For standard error function

# Load raw trait data - for the code below to work, trait value should be called 'trait' & temperature should be called 'T'
data.all <- read.csv("data/UpdatedDataforErin.csv")
#NAs in the dataset. replace NAs with 0
data.all[is.na(data.all)] <- 0

#######Generate Block specific data values for fitting the curve over for comparison.
library(dplyr)
data.block.temp <- data.all %>%
  dplyr::group_by(temp,block)%>%
  dplyr::summarise(bite.rate = mean(bite.rate),
                   lifespan = mean(lifespan),
                   lifetime.eggs = mean(lifetime.eggs),
                   EFD.1stGC = mean(EFD.1stGC)) %>%
  ungroup()

data.no0 <- data.all[!(data.all$est.bite == 0),]
est.bite.means <- data.no0 %>%
  dplyr::group_by(temp,block) %>%
  dplyr::summarise(est.bite = mean(est.bite)) %>%
  ungroup()

data.mu <- read.csv("data/mu.exp.block.values.csv")
data.mu$inverse.mu = 1/data.mu$Value
data.mu$temp <- data.mu$Temp
data.bc.EIP <- read.csv("data/forErin_ShapiroData.csv")
data.bc.EIP$inverse.EIP50 = 1/data.bc.EIP$EIP50
data.pea.MDR <- read.csv("data/Krijn_Raw_Data.csv")

# Load temperature sequence used to calculate trajectories
load("saved posteriors/temps.Rdata")
N.Temp.xs <-length(Temp.xs)


############ Load traits fits Uniform prios
# Fits from Miazgowicz and Shapiro data
load("saved posteriors/lifespan_quad_withT_means_uniform.Rdata")
uni.lifespan.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
load("saved posteriors/bite.rate_briere_withT_means_uniform.Rdata")
uni.bite.rate.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
load("saved posteriors/lifeeggs_briere_withT_means_uniform.Rdata")
uni.lifetime.eggs.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
load("saved posteriors/est.EFD_briere_withT_means_uniform.Rdata")
uni.EFD.1stGC.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
load("saved posteriors/est.bite_briere_withT_means_uniform.Rdata")
uni.est.bite.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
load("saved posteriors/est.lifespan_quad_withT_all_uniform.Rdata")
uni.mu.exp.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred  # actually fits for 1/mu
load("saved posteriors/shapiro_PDR_briere_withT_uniform.Rdata")
uni.EIP.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred  # actually fits for 1/EIP50
load("saved posteriors/shapiro_bc_quad_withT_uniform.Rdata")
uni.bc.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
#Note: to prevent overlap in preds variable names I add a "K" after pea and MDR
load("saved posteriors/krijn_pea_quad_withT_uniform.Rdata")
uni.peaK.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
load("saved posteriors/krijn_MDR_briere_withT_uniform.Rdata")
uni.MDRK.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred

############ Load traits fits from informative priors
# Fits from Miazgowicz and Shapiro data informative
load("saved posteriors inf/lifespan_quadratic_inf.Rdata")
inf.lifespan.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
load("saved posteriors inf/bite.rate_briere_inf.Rdata")
inf.bite.rate.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
load("saved posteriors inf/lifetime.eggs_briere_inf.Rdata")
inf.lifetime.eggs.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
load("saved posteriors inf/EFD1stGC_briere_inf.Rdata")
inf.EFD.1stGC.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred

load("saved posteriors inf/est.bite_briere_inf.Rdata")
inf.est.bite.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
load("saved posteriors inf/mu_exponential_quadratic_inf.Rdata")
inf.mu.exp.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred  # actually fits for 1/mu
load("saved posteriors inf/shapiro_PDR_briere_inf.Rdata")
inf.EIP.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred  # actually fits for 1/EIP50
load("saved posteriors inf/shapiro_bc_quadratic_inf.Rdata")
inf.bc.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
#Note: to prevent overlap in preds variable names I add a "K" after pea and MDR
load("saved posteriors inf/peaK_quadratic_inf.Rdata")
inf.peaK.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
load("saved posteriors inf/MDRK_briere_inf.Rdata")
inf.MDRK.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred

#############
####### 2. Plot thermal performance curves for each trait with both uniform and informative priors
#############
#put the inf fits second
# Plotting function to combine data and TPCs using shaded area instead of lines for the CI

trait.plots.meanse = function(traitname, fitname, df, labs, add = "false", cols = "black"){
  # specify which data to plot
  dat = data.frame('temp' = df$temp, 'trait' = df[ , which(colnames(df)==traitname)])
  colnames(dat)[2] <- "trait"
  traitplot = df[ , which(colnames(df)==traitname)]
  ifelse(add == "true",load(paste("saved posteriors inf/", fitname, sep = "")),load(paste("saved posteriors/", fitname, sep = "")))
  model.fit = model.out
  
  # set the max y value for plotting to the max CI value
  ymax = max(model.fit$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"])
  
  # if you want to make a new plot, use add = "false"
  # if you want to add to an existing plot, use add = "true"
  # plot data means
  if (add == "false"){
    plot(trait ~ temp, data = dat, xlim = c(0, 45), ylim = c(0, ymax), pch = 16, 
         ylab = labs, xlab = expression(paste("Temperature (",degree,"C)")), cex.lab = 1.4, cex.axis = 1.4, col = "black")
  } else {
    points(trait ~ temp, data = dat, pch = 16, col = "black")
  }
  # add TPC mean and 95% CI
  polygon(c(Temp.xs, rev(Temp.xs)), c(model.fit$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(model.fit$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"])),  col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  lines(model.fit$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = cols)
}




# Plot thermal performance curves with only mean and CI ; using shaded area instead of lines to display the CI
###########
pdf("plots/04122020_traitsUNIvsINFPriors.pdf", width = 8, height = 10)
par(mfrow = c(4,3))
trait.plots.meanse("bite.rate", "bite.rate_briere_withT_means_uniform.Rdata " ,data.block.temp, "Bite rate (a)" )
trait.plots.meanse("bite.rate",  "bite.rate_briere_inf.Rdata", data.block.temp, "Bite rate (a)", add = "true", cols = "#C75000")
legend('none', legend = c("uniform", "informative"), lty = 1, col = c("black", "#C75000"), bty = 'n')
mtext("A", side = 3, at = -5, cex = 1.8, line = 2)

trait.plots.meanse("lifespan", "lifespan_quad_withT_means_uniform.Rdata", data.block.temp, "Lifespan (lf)")
trait.plots.meanse("lifespan", "lifespan_quadratic_inf.Rdata", data.block.temp, "Lifespan (lf)", add = "true", cols = "#C75000")
legend('none', legend = c("uniform", "informative"), lty = 1, col = c("black", "#C75000"), bty = 'n')
mtext("B", side = 3, at = -5, cex = 1.8, line = 2)

trait.plots.meanse("lifetime.eggs", "lifeeggs_briere_withT_means_uniform.Rdata", data.block.temp, "Lifetime egg production (B)")
#trait.plots.meanse("lifetime.eggs", "lifetime.eggs_briere_inf.Rdata", data.block.temp, "Lifetime egg production (B)", add = "true", cols = "#C75000")
legend('none', legend = c("uniform", "informative"), lty = 1, col = c("black", "#C75000"), bty = 'n')
mtext("C", side = 3, at = -5, cex = 1.8, line = 2)

trait.plots.meanse("est.bite", "est.bite_briere_withT_means_uniform.Rdata", est.bite.means, "Estimated biting rate (a*)")
trait.plots.meanse("est.bite","est.bite_briere_inf.Rdata" , est.bite.means, "Estimated biting rate (a*)", add = "true", cols = "#C75000")
legend('none', legend = c("uniform", "informative"), lty = 1, col = c("black", "#C75000"), bty = 'n')
mtext("D", side = 3, at = -5, cex = 1.8, line = 2)

trait.plots.meanse("inverse.mu", "est.lifespan_quad_withT_all_uniform.Rdata", data.mu, "Estimated lifespan (lf*)")
trait.plots.meanse("inverse.mu", "mu_exponential_quadratic_inf.Rdata", data.mu, "Estimated lifespan (lf*) ", add = "true", cols = "#C75000")
legend('none', legend = c("uniform", "informative"), lty = 1, col = c("black", "#C75000"), bty = 'n')
mtext("E", side = 3, at = -5, cex = 1.8, line = 2)

trait.plots.meanse("EFD.1stGC", "est.EFD_briere_withT_means_uniform.Rdata", data.block.temp, "Estimated daily eggs (EFD*)")
trait.plots.meanse("EFD.1stGC","EFD1stGC_briere_inf.Rdata" , data.block.temp, "Estimated daily eggs (EFD*')", add = "true", cols = "#C75000")
legend('none', legend = c("uniform", "informative"), lty = 1, col = c("black", "#C75000"), bty = 'n')
mtext("F", side = 3, at = -5, cex = 1.8, line = 2)

trait.plots.meanse("bc", "shapiro_bc_quad_withT_uniform.Rdata", data.bc.EIP, "Vector competence (bc)")
trait.plots.meanse("bc", "shapiro_bc_quadratic_inf.Rdata", data.bc.EIP, "Vector competence (bc)", add = "true", cols = "#C75000")
legend('none', legend = c("uniform", "informative"), lty = 1, col = c("black", "#C75000"), bty = 'n')
mtext("G", side = 3, at = -5, cex = 1.8, line = 2)

trait.plots.meanse("inverse.EIP50","shapiro_PDR_briere_withT_uniform.Rdata", data.bc.EIP, "Parasite development rate (PDR)")
trait.plots.meanse("inverse.EIP50",  "shapiro_PDR_briere_inf.Rdata", data.bc.EIP, "Parasite development rate (PDR)", add ="true", cols = "#C75000")
legend('none', legend = c("uniform", "informative"), lty = 1, col = c("black", "#C75000"), bty = 'n')
mtext("H", side = 3, at = -5, cex = 1.8, line = 2)

trait.plots.meanse("Pea","krijn_pea_quad_withT_uniform.Rdata", data.pea.MDR, "Prob egg2adult survival (pEA)")
trait.plots.meanse("Pea",  "peaK_quadratic_inf.Rdata", data.pea.MDR, "Prob egg2adult survival (pEA)", add ="true", cols = "#C75000")
legend('none', legend = c("uniform", "informative"), lty = 1, col = c("black", "#C75000"), bty = 'n')
mtext("I", side = 3, at = -5, cex = 1.8, line = 2)

trait.plots.meanse("MDR","krijn_MDR_briere_withT_uniform.Rdata", data.pea.MDR, "Mosquito development rate (MDR)")
trait.plots.meanse("MDR",  "MDRK_briere_inf.Rdata", data.pea.MDR, "Mosquito development rate (MDR)", add ="true", cols = "#C75000")
legend('none', legend = c("uniform", "informative"), lty = 1, col = c("black", "#C75000"), bty = 'n')
mtext("J", side = 3, at = -5, cex = 1.8, line = 2)

#add in panels for An. stephensi estimated
#For this to work; need to run the Calc_R0_uniform code
#For this to work; need to run the Calc_R0_informative code
plot.R0.poly(df = R0.6.uni.out, add = "false", cols = "black" )
plot.R0.poly(df = R0.6.out, add = "true", cols = "#C75000")
legend('none', legend = c(expression(italic(An.~stephensi)~estimated),"uniform", "informative"), lty = 1, col = c(NULL,"black", "#C75000"), bty = 'n')
mtext("K", side = 3, at = -5, cex = 1.8, line = 2)


#add in panels for An. stephensi lifetime
#For this to work; need to run the Calc_R0_uniform code
#For this to work; need to run the Calc_R0_informative code
plot.R0.poly(df = R0.7.uni.out, add = "false", cols = "black" )
plot.R0.poly(df = R0.7.out, add = "true", cols = "#C75000")
legend('none', legend = c(expression(italic(An.~stephensi)~lifetime),"uniform", "informative"), lty = 1, col = c(NULL,"black", "#C75000"), bty = 'n')
mtext("L", side = 3, at = -5, cex = 1.8, line = 2)

par(mfrow = c(1,1))
dev.off()

