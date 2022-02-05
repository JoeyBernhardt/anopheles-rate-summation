## Erin Mordecai, Stanford University
## August 8, 2018
## Updated by Kerri Miazgowicz on November 19th, 2019 to include the truncated model fits.
## Updated by KM on April 13, 2020 to revise for resubmission
## Purpose: Plot trait thermal response comparisons
##
## Contents: 1) Set-up,load packages, get data, etc.
##           1b) Use the model mean function parameters to calculate the max trait value
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

# Load raw trait data - for the code below to work, trait value should be called 'trait' & temperature should be called 'T'
data.all <- read.csv("data/UpdatedDataforErin.csv")

#Remove 0's from est.bite rate
data.no0 <- data.all[!(data.all$est.bite == 0),]
#######Generate Block specific data values for fitting the curve over for comparison - removed 0s
data.block.temp.est.bite <- data.no0 %>%
  dplyr::group_by(temp,block)%>%
  dplyr::summarise(est.bite = mean(est.bite)) %>%
  ungroup()

data.block.temp.est.bite <- as.data.frame(data.block.temp.est.bite)

#######Generate Block specific data values for fitting the curve over for comparison.
data.block.temp <- data.all %>%
  dplyr::group_by(temp,block)%>%
  dplyr::summarise(bite.rate = mean(bite.rate),
                   lifespan = mean(lifespan),
                   lifetime.eggs = mean(lifetime.eggs),
                   EFD.1stGC = mean(EFD.1stGC)) %>%
  ungroup()

data.block.temp <- as.data.frame(data.block.temp)

# Load raw trait data - for the code below to work, trait value should be called 'trait' & temperature should be called 'T'
data.mu <- read.csv("data/mu.exp.block.values.csv")
#convert to lifepan estimates
data.mu$lifespan <- 1/data.mu$Value
data.mu$temp <- data.mu$Temp
data.mu$block <- data.mu$Block

data.gamma <- read.csv("data/Miaz_gamma.values.csv")

# Load raw trait data - for the code below to work, trait value should be called 'trait' & temperature should be called 'T'
data.pea.MDR <- read.csv("data/Krijn_Raw_Data.csv")
#replace Na with 0s
data.pea.MDR[is.na(data.pea.MDR)] = 0

# Load raw trait data - for the code below to work, trait value should be called 'trait' & temperature should be called 'T'
data.bc.EIP <- read.csv("data/forErin_ShapiroData.csv")
data.bc.EIP$inverse.EIP50 = 1/data.bc.EIP$EIP50

# Load temperature sequence used to calculate trajectories
load("saved posteriors/temps.Rdata")
N.Temp.xs <-length(Temp.xs)

############ Load traits fits
#########NOTE: this time as I load I'm going to extract out the model parameters

#Lifetime associated traits
load("saved posteriors inf/lifespan_quadratic_inf.Rdata")
lifespan.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
lifespan.params <- model.out$BUGSoutput$summary[1:5,1:9]
#Lifespan mean function -Quad; 1st row-T0, 2nd row-Tm, 3rd row-q
lifespan.values.fun <- function(x){  -1 *lifespan.params[3,1] * (x - lifespan.params[1,1]) * (x - lifespan.params[2,1]) * (lifespan.params[2,1] > x) * (lifespan.params[1,1] < x)}
#     -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
#Calc values from 0,50 and find the max trait value
max.lifespan <- max(lifespan.values.fun(seq(0,50,0.2))) #38.39

load("saved posteriors inf/bite.rate_briere_inf.Rdata")
bite.rate.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
bite.rate.params <- model.out$BUGSoutput$summary[1:5,1:9]
#bite.rate mean functio - Briere
## 1st row T0, 2nd row Tm, 3rd row q
bite.rate.values.fun <- function(x){bite.rate.params[3,1] * x * (x - bite.rate.params[1,1]) * sqrt((bite.rate.params[2,1] - x) * (bite.rate.params[2,1] > x)) * (bite.rate.params[1,1] < x)}
max.bite.rate <- max(bite.rate.values.fun(seq(0,50,0.2)), na.rm = T)

#Note: due to note having an appropiate trait to fit informative priors to we use the uniformative fits
load("saved posteriors/lifeeggs_briere_withT_means_uniform.Rdata")
lifetime.eggs.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
lifetime.eggs.params <- model.out$BUGSoutput$summary[1:5,1:9]
lifetime.eggs.values.fun <- function(x){lifetime.eggs.params[3,1] * x * (x - lifetime.eggs.params[1,1]) * sqrt((lifetime.eggs.params[2,1] - x) * (lifetime.eggs.params[2,1] > x)) * (lifetime.eggs.params[1,1] < x)}
max.lifetime.eggs <- max(lifetime.eggs.values.fun(seq(0,50,0.2)), na.rm = T)



#Estimated associated traits
#Note: due to the high amount of uncertainty on the inf fit; we are using the uni for this trait
load("saved posteriors/est.EFD_briere_withT_means_uniform.Rdata")
EFD.1stGC.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
EFD.1stGC.params <- model.out$BUGSoutput$summary[1:5,1:9]
EFD.1stGC.values.fun <- function(x){EFD.1stGC.params[3,1] * x * (x - EFD.1stGC.params[1,1]) * sqrt((EFD.1stGC.params[2,1] - x) * (EFD.1stGC.params[2,1] > x)) * (EFD.1stGC.params[1,1] < x)}


load("saved posteriors inf/est.bite_briere_inf.Rdata")
est.bite.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
est.bite.params <- model.out$BUGSoutput$summary[1:5,1:9]
est.bite.values.fun <- function(x){est.bite.params[3,1] * x * (x - est.bite.params[1,1]) * sqrt((est.bite.params[2,1] - x) * (est.bite.params[2,1] > x)) * (est.bite.params[1,1] < x)}
max.est.rate <- max(est.bite.values.fun(seq(0,50,0.2)), na.rm = T)


load("saved posteriors inf/mu_exponential_quadratic_inf.Rdata")
mu.exp.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred  # actually fits for 1/mu
mu.exp.params <- model.out$BUGSoutput$summary[1:5,1:9]
#est.lifespan mean function -Quad; 1st row-T0, 2nd row-Tm, 3rd row-q
est.lifespan.values.fun <- function(x){  -1 *mu.exp.params[3,1] * (x - mu.exp.params[1,1]) * (x - mu.exp.params[2,1]) * (mu.exp.params[2,1] > x) * (mu.exp.params[1,1] < x)}
#     -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
#Calc values from 0,50 and find the max trait value
max.est.lifespan <- max(est.lifespan.values.fun(seq(0,50,0.2))) #38.39

#Calculate the max of the estimated lifetime egg production curve
max.est.lifetime.eggs <- max(est.lifespan.values.fun(seq(0,50,0.2))*EFD.1stGC.values.fun(seq(0,50,0.2)))


#Shapiro traits -> for determining the TPC thresholds
load("saved posteriors inf/shapiro_PDR_briere_inf.Rdata")
EIP.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred  # actually fits for 1/EIP50
EIP.params <- model.out$BUGSoutput$summary[1:5,1:9]

load("saved posteriors inf/shapiro_bc_quadratic_inf.Rdata")
bc.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
bc.params <- model.out$BUGSoutput$summary[1:5,1:9]

#Krijn traits -> for determining the TPC thresholds
#Note: to prevent overlap in preds variable names I add a "K" after pea and MDR
load("saved posteriors inf/peaK_quadratic_inf.Rdata")
peaK.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
peaK.params <- model.out$BUGSoutput$summary[1:5,1:9]

load("saved posteriors inf/MDRK_briere_inf.Rdata")
MDRK.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
MDRK.params <- model.out$BUGSoutput$summary[1:5,1:9]

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
all.param.list = list(c("lifespan.params", "bite.rate.params","lifetime.eggs.params","EFD.1stGC.params","est.bite.params","mu.exp.params","EIP.params","bc.params","peaK.params","MDRK.params"))
trait.params = lapply(all.param.list, param.compile.fun)[[1]]
names(trait.params) <- all.param.list[[1]]
trait.params

#Export this for quick reference
#sink("output_final_TPC_params.txt")
#print(trait.params)
#sink()

#############
####### 2. Plot thermal performance curves
#############

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


##############################################
##First I need to redefine the mapping function
###########################################

# Plot thermal performance curves with data they were fit over
pdf("plots/trait_fits_final_means.pdf", width = 8, height = 10)
par(mfrow = c(4,3))
trait.plots("bite.rate", "bite.rate_briere_inf.Rdata", data.block.temp, "Biting rate", cols = "black", pch = 1)
mtext("A", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots("lifespan", "lifespan_quadratic_inf.Rdata", data.block.temp, "Lifespan", cols = c("black"), pch = 1)
mtext("B", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots("lifetime.eggs", "lifeeggs_briere_withT_means_uniform.Rdata", data.block.temp, "Lifetime egg production", cols = c("black"),pch=1, source = "uni")
mtext("C", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots("est.bite", "est.bite_briere_inf.Rdata", data.block.temp.est.bite, "Estimated biting rate", cols = c("dodgerblue"), pch = 1)
mtext("D", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots("lifespan", "mu_exponential_quadratic_inf.Rdata", data.mu, "Estimated lifespan", cols = c("dodgerblue"),pch =1)
mtext("E", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots("EFD.1stGC", "est.EFD_briere_withT_means_uniform.Rdata", data.block.temp, "Estimated daily eggs", cols = c("dodgerblue"), pch =1, source = "uni")
mtext("F", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots("bc", "shapiro_bc_quadratic_inf.Rdata", data.bc.EIP, "Vector competence", cols = c("seagreen3"), pch=1)
mtext("G", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots("inverse.EIP50", "shapiro_PDR_briere_inf.Rdata", data.bc.EIP, "Parasite development rate", cols = c("seagreen3"), pch =1)
mtext("H", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots("Pea", "peaK_quadratic_inf.Rdata", data.pea.MDR, "Prob. egg to adult survival", cols = c("seagreen3"), pch=1)
mtext("I", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots("MDR", "MDRK_briere_inf.Rdata", data.pea.MDR, "Mosquito development rate", cols = c("seagreen3"), pch=1)
mtext("J", side = 3, at = -5, cex = 1.8, line = 2)
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
#To be able to properly compare lifetime egg production verses estimated lifetime egg production (which is EFD (1st GC)* inv(Mu.exp method))
# I need to create a function that generates this new curve and the associated CI - similar to how it is done with R0

### function to define estimated lifetime egg production:
# constant to keep lifespan from being numerically zero
# assume minimum survival time is half an hour
ec.lf = 1/48

est.lifetime.eggs = function(EFD, inverse.mu){
  (EFD*(inverse.mu+ec.lf))
}

######  Calculate estimated lifetime egg production posteriors and quantiles

# Calculate full posteriors and quantiles

estimated.lifeeggs = est.lifetime.eggs(EFD.1stGC.preds, mu.exp.preds)
estimated.lifeeggs.out = calcPostQuants(estimated.lifeeggs, Temp.xs)
estimated.lifeeggs.out$temp = Temp.xs
estimated.lifeeggs.TPdists = calcThreshPeakDists(estimated.lifeeggs, Temp.xs)


##############
#Calculate post Quants for the other traits so that I can find the corresponding Bayesian fit value for all the other traits. Give in preds and temp vector
bite.rate.out = calcPostQuants(bite.rate.preds, Temp.xs)
bite.rate.out$temp = Temp.xs
estimated.bite.rate.out = calcPostQuants(est.bite.preds, Temp.xs)
estimated.bite.rate.out$temp = Temp.xs

lifespan.out = calcPostQuants(lifespan.preds, Temp.xs)
lifespan.out$temp = Temp.xs
estimated.lifespan.out = calcPostQuants(mu.exp.preds, Temp.xs)
estimated.lifespan.out$temp = Temp.xs

lifetime.eggs.out = calcPostQuants(lifetime.eggs.preds,Temp.xs)
lifetime.eggs.out$temp = Temp.xs


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


###########
pdf("plots/trait_lifetime_estimated_inf_versionB.pdf", width = 4.5, height = 11)
par(mfrow = c(3,1))
# Biting rate
trait.plots("bite.rate", "bite.rate_briere_inf.Rdata", data.block.temp, "Biting rate (a)", cols = "black", pch =16)
trait.plots("est.bite", "est.bite_briere_inf.Rdata", data.block.temp.est.bite, "Estimated biting rate", add = "true", cols = "dodgerblue", pch =16)
legend('topleft', legend = c("observed (a)", "estimated (a*)"), lty = 1, col = c("black", "dodgerblue"), bty = 'n')
mtext("a", side = 3, at = -5, cex = 1.8, line = 2)

# Lifespan
trait.plots("lifespan", "lifespan_quadratic_inf.Rdata", data.block.temp, "Lifespan (lf)", cols = "black",pch =16)
trait.plots("lifespan", "mu_exponential_quadratic_inf.Rdata", data.mu, "Estimated lifespan", add = "true", cols = "dodgerblue", pch =16)
legend('topleft', legend = c("observed (lf)", "estimated (lf*)"), lty = 1, col = c("black", "dodgerblue"), bty = 'n')
mtext("b", side = 3, at = -5, cex = 1.8, line = 2)

# Egg production
# We'll have to do it differently without showing data points, 
# since we have to combine fits from multiple traits
trait.plots("lifetime.eggs", "lifeeggs_briere_withT_means_uniform.Rdata", data.block.temp, "Lifetime egg production (B)",cols="black", source = "uni",pch =16)
plot.trait2(df = estimated.lifeeggs.out, add = "true", cols = "dodgerblue")
legend('topleft', legend = c("observed (B)", "estimated (B*)"), lty = 1, col = c("black", "dodgerblue"), bty = 'n')
mtext("c", side = 3, at = -5, cex = 1.8, line = 2)

par(mfrow = c(1,1))
dev.off()
############

##############################################Need to run the R0 calc R script before running the below code
###############################################


# Alternative function that uses shaded polygons instead of dashed lines for R0
plot.R0.poly = function(temp = Temp.xs, df, modname = "", add = "false", cols = "black", ltys = 1){
  if (add == "false") {
    par(mar=c(5,6,5,1))
    plot(temp, df$upperCI/max(df$upperCI, na.rm = T), lty = 1, xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, ylab = expression("Relative R"[0]), yaxt = 'n', col = "white")
    polygon(c(temp, rev(temp)), c(df$upperCI/max(df$upperCI, na.rm = T), rev(df$lowerCI/max(df$upperCI, na.rm = T))), col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  } else {
    polygon(c(temp, rev(temp)), c(df$upperCI/max(df$upperCI, na.rm = T), rev(df$lowerCI/max(df$upperCI, na.rm = T))), col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  }
  
  lines(temp, df$mean/max(df$upperCI, na.rm = T), lty = ltys, lwd = 2, col = cols)
}


##################
###############################3
##########Main text Figure 2- the three trait plot overlays along with R0#############3
#Note for this to run properly the according code from R0_calc_inf.R needs to be run prior
#Required R0.5.out, R0.6.out, R0.7.out from the final fits with mostly informative priors

pdf("plots/Main_Figure2_lifetime_estimated_final_.pdf", width = 13.3333, height = 16)
par(mfrow = c(4,3))

# Biting rate
trait.plots("bite.rate", "bite.rate_briere_inf.Rdata", data.block.temp, "Biting rate (a)", cols = "black", pch =16)
trait.plots("est.bite", "est.bite_briere_inf.Rdata", data.block.temp.est.bite, "Estimated biting rate", add = "true", cols = "dodgerblue", pch =16)
legend('topleft', legend = c("observed (a)", "estimated (a*)"), lty = 1, col = c("black", "dodgerblue"), bty = 'n')
mtext("a", side = 3, at = -5, cex = 1.8, line = 2)

# Lifespan
trait.plots("lifespan", "lifespan_quadratic_inf.Rdata", data.block.temp, "Lifespan (lf)", cols = "black",pch =16)
trait.plots("lifespan", "mu_exponential_quadratic_inf.Rdata", data.mu, "Estimated lifespan", add = "true", cols = "dodgerblue", pch =16)
legend('topleft', legend = c("observed (lf)", "estimated (lf*)"), lty = 1, col = c("black", "dodgerblue"), bty = 'n')
mtext("b", side = 3, at = -5, cex = 1.8, line = 2)

plot(NULL)

# Egg production
# We'll have to do it differently without showing data points, 
# since we have to combine fits from multiple traits
trait.plots("lifetime.eggs", "lifeeggs_briere_withT_means_uniform.Rdata", data.block.temp, "Lifetime egg production (B)",cols="black", source = "uni",pch =16)
plot.trait2(df = estimated.lifeeggs.out, add = "true", cols = "dodgerblue")
legend('topleft', legend = c("observed (B)", "estimated (B*)"), lty = 1, col = c("black", "dodgerblue"), bty = 'n')
mtext("c", side = 3, at = -5, cex = 1.8, line = 2)


plot.R0.poly(df = R0.5.out, add = "false", cols = "purple")
plot.R0.poly(df = R0.6.out, add = "true", cols = "black" )
plot.R0.poly(df = R0.7.out, add = "true", cols = "dodgerblue")
legend('topleft', legend = c("Multi-species estimated", expression(italic(An.~stephensi)~lifetime),expression(italic(An.~stephensi)~estimated)), lty = c(1,1,1), lwd = 2, col = c("purple","black", "dodgerblue"), bty = 'n')
mtext("d", side = 3, at = -5, cex = 1.8, line = 2)
par(mfrow = c(1,1))
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
traits = list("bite.rate.preds", "lifespan.preds", "lifetime.eggs.preds", "est.bite.preds", "mu.exp.preds", "EFD.1stGC.preds", "bc.preds", "EIP.preds", "peaK.preds", "MDRK.preds")
trait.summaries = unlist(lapply(traits, meta.summary.TPC), recursive = F)
names(trait.summaries) = traits
trait.summaries
#Export this list for quick reference later on

#sink("output_final_trait_thresholds.txt")
#print(trait.summaries)
#sink()

#Deterimine T-thresholds for estimated lifetime egg production (was product of two other traits, evaluate like R0)
meta.summary.fun = function(txt) {
  out = list()
  for (i in 1:length(txt)){
    df = eval(parse(text = txt[i]))
    out[[i]] = apply(df, 2, summary.fn)
  }
  out
}

#looking at thermal limits for estimated.lifeeggs
trait.subset = list(c("estimated.lifeeggs.TPdists"))
trait.subset.summaries = lapply(trait.subset, meta.summary.fun)[[1]]
names(trait.subset.summaries) = trait.subset[[1]]
trait.subset.summaries
#export this as well for quick reference later
#sink("output_final_est.lifeeggs_thresholds.txt")
#print(trait.subset.summaries)
#sink()

##################
pdf("plots/trait_summaries_inf.pdf", width = 7, height = 5)
plot.summaries(trait.summaries, c("biting rate", "lifespan", "lifetime eggs", "est. biting rate", "est. lifespan", "EFD 1st GC", "vector competence", "inverse of EIP50", "Prob. egg to adult survival", "Mosqutio development rate"), xlims = c(0, 45), cols = c(1, "royalblue", "firebrick"))
dev.off()
##################

#########################################Supplmental Figuure
####Posterior distributions of the peak, Tmax, and Tmin for each trait in Figure 2.
##panel A Tmin

#bite rate
br.TPdists <- calcThreshPeakDists(bite.rate.preds, Temp.xs)
est.br.TPdists <- calcThreshPeakDists(est.bite.preds,Temp.xs)
#lifespan
lf.TPdists <- calcThreshPeakDists(lifespan.preds, Temp.xs)
est.lf.TPdists <- calcThreshPeakDists(mu.exp.preds, Temp.xs)
#lifetime egg production
life.eggs.TPdists <- calcThreshPeakDists(lifetime.eggs.preds, Temp.xs)
#estimated.lifeeggs.TPdists <- calculated in the above code

break.n <- seq(0,50,0.4)
#Tmin
p1 <- hist(br.TPdists$T0, breaks = break.n,freq= F)
p2 <- hist(est.br.TPdists$T0, breaks = break.n, freq = F)
p3 <- hist(lf.TPdists$T0, breaks = break.n,freq=F)
p4 <- hist(est.lf.TPdists$T0, breaks = break.n, freq=F)
p5 <- hist(life.eggs.TPdists$T0, breaks = break.n, freq = F)
p6 <- hist(estimated.lifeeggs.TPdists$T0, breaks = break.n, freq = F)
#Peak
p7 <- hist(br.TPdists$peak, breaks = break.n,freq= F)
p8 <- hist(est.br.TPdists$peak, breaks = break.n, freq = F)
p9 <- hist(lf.TPdists$peak, breaks = break.n,freq=F)
p10 <- hist(est.lf.TPdists$peak, breaks = break.n, freq=F)
p11 <- hist(life.eggs.TPdists$peak, breaks = break.n, freq = F)
p12 <- hist(estimated.lifeeggs.TPdists$peak, breaks = break.n, freq = F)
#Tmax
p13 <- hist(br.TPdists$Tmax, breaks = break.n,freq= F)
p14 <- hist(est.br.TPdists$Tmax, breaks = break.n, freq = F)
p15 <- hist(lf.TPdists$Tmax, breaks = break.n,freq=F)
p16 <- hist(est.lf.TPdists$Tmax, breaks = break.n, freq=F)
p17 <- hist(life.eggs.TPdists$Tmax, breaks = break.n, freq = F)
p18 <- hist(estimated.lifeeggs.TPdists$Tmax, breaks = break.n, freq = F)

###############################Ok make a seperate plot for each trait
pdf("plots/Thresholds_density_a_B_lf.pdf", width = 8, height = 12)
par(mar=c(5,6,5,1))
par(mfrow = c(3,3))

plot(p1, col =scales::alpha('grey43',0.4),freq= F,border=scales::alpha('black',0.6), xlim= c(0,20),ylim = c(0,0.3),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, main = expression(italic(paste("T"[min])- a)))
plot(p2, col= scales::alpha('dodgerblue',0.4),freq= F,border=scales::alpha('royalblue',0.6), ylab = "Density",add = T)
legend('topleft', legend = c(expression(italic(An.~stephensi)~lifetime),expression(italic(An.~stephensi)~estimated)), lty = 1, col = c("black", "dodgerblue"), bty = 'n')
mtext("a", side = 3, at = -7, cex = 1.8, line = 2)


plot(p7, col =scales::alpha('grey43',0.4),freq= F,border=scales::alpha('black',0.6), xlim= c(32,40),ylim = c(0,0.6),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, main = expression(italic(paste("T"[opt])-a)))
plot(p8, col= scales::alpha('dodgerblue',0.4),freq= F,border=scales::alpha('royalblue',0.6), ylab = "Density",add = T)
legend('topleft', legend = c(expression(italic(An.~stephensi)~lifetime),expression(italic(An.~stephensi)~estimated)), lty = 1, col = c("black", "dodgerblue"), bty = 'n')

plot(p13, col =scales::alpha('grey43',0.4),freq= F,border=scales::alpha('black',0.6), xlim= c(38,46),ylim = c(0,0.425),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, main = expression(italic(paste("T"[max])-a)))
plot(p14, col= scales::alpha('dodgerblue',0.4),freq= F,border=scales::alpha('royalblue',0.6), ylab = "Density",add = T)
legend('topleft', legend = c(expression(italic(An.~stephensi)~lifetime),expression(italic(An.~stephensi)~estimated)), lty = 1, col = c("black", "dodgerblue"), bty = 'n')

plot(p3, col =scales::alpha('grey43',0.4),freq= F,border=scales::alpha('black',0.6), xlim= c(0,6),ylim = c(0,.6),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, main = expression(italic(paste("T"[min])-lf)))
plot(p4, col= scales::alpha('dodgerblue',0.4),freq= F,border=scales::alpha('royalblue',0.6), ylab = "Density",add = T)
legend('topleft', legend = c(expression(italic(An.~stephensi)~lifetime),expression(italic(An.~stephensi)~estimated)), lty = 1, col = c("black", "dodgerblue"), bty = 'n')
mtext("b", side = 3, at = -2, cex = 1.8, line = 2)

plot(p9, col =scales::alpha('grey43',0.4),freq= F,border=scales::alpha('black',0.6), xlim= c(17,22),ylim = c(0,0.7),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, main = expression(italic(paste("T"[opt])-lf)))
plot(p10, col= scales::alpha('dodgerblue',0.4),freq= F,border=scales::alpha('royalblue',0.6), ylab = "Density",add = T)
legend('topleft', legend = c(expression(italic(An.~stephensi)~lifetime),expression(italic(An.~stephensi)~estimated)), lty = 1, col = c("black", "dodgerblue"), bty = 'n')

plot(p15, col =scales::alpha('grey43',0.4),freq= F,border=scales::alpha('black',0.6), xlim= c(32,42),ylim = c(0,0.6),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, main = expression(italic(paste("T"[max])-lf)))
plot(p16, col= scales::alpha('dodgerblue',0.4),freq= F,border=scales::alpha('royalblue',0.6), ylab = "Density",add = T)
legend('topleft', legend = c(expression(italic(An.~stephensi)~lifetime),expression(italic(An.~stephensi)~estimated)), lty = 1, col = c("black", "dodgerblue"), bty = 'n')

plot(p5, col =scales::alpha('grey43',0.4),freq= F,border=scales::alpha('black',0.3), xlim= c(0,20),ylim = c(0,.1),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, main = expression(italic(paste("T"[min])-B)))
plot(p6, col= scales::alpha('dodgerblue',0.4),freq= F,border=scales::alpha('royalblue',0.6), ylab = "Density",add = T)
legend('topleft', legend = c(expression(italic(An.~stephensi)~lifetime),expression(italic(An.~stephensi)~estimated)), lty = 1, col = c("black", "dodgerblue"), bty = 'n')
mtext("c", side = 3, at = -7, cex = 1.8, line = 2)


plot(p11, col =scales::alpha('grey43',0.4),freq= F,border=scales::alpha('black',0.6), xlim= c(22,32),ylim = c(0,0.6),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, main = expression(italic(paste("T"[opt])-B)))
plot(p12, col= scales::alpha('dodgerblue',0.4),freq= F,border=scales::alpha('royalblue',0.6), ylab = "Density",add = T)
legend('topleft', legend = c(expression(italic(An.~stephensi)~lifetime),expression(italic(An.~stephensi)~estimated)), lty = 1, col = c("black", "dodgerblue"), bty = 'n')

plot(p17, col =scales::alpha('grey43',0.4),freq= F,border=scales::alpha('black',0.6), xlim= c(26,45),ylim = c(0,.4),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, main = expression(italic(paste("T"[max])-B)))
plot(p18, col= scales::alpha('dodgerblue',0.4),freq= F,border=scales::alpha('royalblue',0.6), ylab = "Density",add = T)
legend('topleft', legend = c(expression(italic(An.~stephensi)~lifetime),expression(italic(An.~stephensi)~estimated)), lty = 1, col = c("black", "dodgerblue"), bty = 'n')
par(mfrow = c(1,1))
dev.off()
