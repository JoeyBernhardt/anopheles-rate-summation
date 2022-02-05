## Marta Shocket, Stanford University
## Updated by Erin Mordecai, August 8, 2018
## Updated by Kerri Miazgowicz, November 19,2019 
## Updated by KM April 12,2020
##
## Purpose: Use Bayesian Inference (JAGS) to fit temperature-dependent functions traits using informative priors
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) JAGS models
##           3) Shared settings for all models
##           4) Code for fitting prior distributions from the Johnson et al TPCs
##           5) Code for fitting a trait - informative priors fit from data

##########
###### 1. Set up workspace, load packages, get data, etc.
##########

mainDir = "C:/Users/Kerri/Desktop/Chapter1 Submission"
setwd(mainDir)

# Check whether there's a folder in the directory for saving posteriors
# If not, create one
subDir = "saved posteriors inf"
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)

# Load libraties for fitting traits
library('R2jags') # Fits Bayesian models
library('mcmcplots') # Diagnostic plots for fits
library('MASS') # Fits distributions for informative priors
library(plyr) # Slices and dices data
library(plotrix) # For standard error function
library(coda)

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

data.block.temp <- as.data.frame(data.block.temp)
data.mu <- read.csv("data/mu.exp.block.values.csv")
data.mu$inverse.mu = 1/data.mu$Value
data.mu$temp <- data.mu$Temp
data.bc.EIP <- read.csv("data/forErin_ShapiroData.csv")
data.bc.EIP$inverse.EIP50 = 1/data.bc.EIP$EIP50
data.pea.MDR <- read.csv("data/Krijn_Raw_Data.csv")


# Load prior fits
load("Johnson et al posteriors/PriorFitsList.Rdata")

##########
###### 2. JAGS Models
##########

# NOTES:
# - Running the code below writes a .txt file for each model to your current working directory.
# - The models include a section for 'derived quanities' that calculates the trait across a
#     temperature gradient for each saved sample in the MCMC chain. This output is what we
#     use later on to calculate R0. 

############## Quadratic Model with gamma priors

sink("quad_T_inf.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dgamma(hypers[1,1], hypers[2,1])
    cf.T0 ~ dgamma(hypers[1,2], hypers[2,2])
    cf.Tm ~ dgamma(hypers[1,3], hypers[2,3])
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)T(0,)
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])
    }
    
    } # close model
    ",fill=T)
sink()

############## Briere Model with gamma priors

sink("briere_T_inf.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dgamma(hypers[1,1], hypers[2,1])
    cf.T0 ~ dgamma(hypers[1,2], hypers[2,2])
    cf.Tm ~ dgamma(hypers[1,3], hypers[2,3])
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- cf.q * temp[i] * (temp[i] - cf.T0) * sqrt((cf.Tm - temp[i]) * (cf.Tm > temp[i])) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)T(0,)
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i])
    }
    
    } # close model
    ",fill=T)
sink()

############## Briere Model with gamma priors for T0 and Tm only

sink("briere_T_inf_T0_Tm.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dgamma(hypers[1,1], hypers[2,1])
    cf.Tm ~ dgamma(hypers[1,2], hypers[2,2])
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- cf.q * temp[i] * (temp[i] - cf.T0) * sqrt((cf.Tm - temp[i]) * (cf.Tm > temp[i])) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)T(0,)
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i])
    }
    
    } # close model
    ",fill=T)
sink()

##########
###### 3. Shared settings for all models
##########

##### inits Function
inits <- function(){list(
  cf.q = 0.01,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1))}

##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")

##### MCMC Settings
# Number of posterior distribution elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 25000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains


# load saved temperature sequence used for uniform prior analyses
load("saved posteriors/temps.Rdata")
N.Temp.xs = length(Temp.xs)

##########
###### 4. Fitting prior distributions from the Johnson et al TPCs
##########

# Function to fit gamma distributions for q, T0, Tm, and sigma from fitted TPCs
prior.fitting = function(model.out){
  # Fit gamma distributions for each parameter posterior dists
  prior.fits = suppressWarnings(apply(model.out, 2, function(df) fitdistr(df, "gamma")$estimate))
  prior.fits
} 

# Load Prior Fits
posts = c(
  "a_posterior1",
  "PDR_post_prior1",
  "MDR_posterior1",
  "EFD_posterior1",
  "e2a_posterior1",
  "bc_posterior1_quad",
  "mu_posterior1")
fits = paste("Johnson et al posteriors/", posts, ".Rsave", sep = "")
for (i in 1:length(fits)){
  load(fits[i])
  assign(posts[i], samps)
}

prior.fits.list = lapply(list(a_posterior1, PDR_post_prior1, MDR_posterior1, EFD_posterior1[,-3], e2a_posterior1[,-3], bc_posterior1_quad[,-3], mu_posterior1), prior.fitting)
names(prior.fits.list) = posts

# Since I fit a TPC for lifespan from the Johnson et al data separately, load that and derive priors
load("saved posteriors/johnson_lifespan_quadratic_uniform.Rdata")
lf.prior.out = model.out

# Get the posterior dists for TPC parameters into a data frame
lf.prior.cf.dists <- data.frame(q = as.vector(lf.prior.out$BUGSoutput$sims.list$cf.q),
                                T0 = as.vector(lf.prior.out$BUGSoutput$sims.list$cf.T0), 
                                Tm = as.vector(lf.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
lf.prior.gamma.fits = apply(lf.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)
prior.fits.list$lf_posterior1 = lf.prior.gamma.fits

save(prior.fits.list, file = "Johnson et al posteriors/PriorFitsList.Rdata")



# ##########
# ###### 5. Code for fitting a trait - informative priors fit from data
# ##########
# 
# ############## Trait name and parameter symbol - model type (e.g. quadratic, Briere)

########## Biting rate (a)
# Get data
data.specific <- with(data.block.temp, data.frame('T' = temp, 'trait' = bite.rate)) # subset specific trait data from complete list
data <- data.specific # assign trait data to variable 'data'
trait.prior.fits = prior.fits.list$a_posterior1[,c('c', 'T0', 'Tm')]
hypers <- trait.prior.fits * 0.5 # assign priors to variable 'priors' and weight importance if desired

# Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS - **select correct model file**
model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_T_inf.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model.out$BUGSoutput$summary[1:10,]
mcmcplot(model.out)

# Save model output
save(model.out, file = "saved posteriors inf/bite.rate_briere_inf.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, xaxt = "n",
     ylab = "Biting rate", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

############ Lifetime eggs
# Get data
data.specific <- with(data.block.temp, data.frame('T' = temp, 'trait' = lifetime.eggs)) # subset specific trait data from complete list
data <- data.specific # assign trait data to variable 'data'
trait.prior.fits = prior.fits.list$EFD_posterior1[,c('n.qd', 'T0', 'Tm')]
hypers <- trait.prior.fits * 0.1 # assign priors to variable 'priors' and weight importance if desired


# Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS - **select correct model file**
model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_T_inf.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model.out$BUGSoutput$summary[1:10,]
mcmcplot(model.out)

# Save model output 
save(model.out, file = "saved posteriors inf/lifetime.eggs_briere_inf.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
     ylab = "Lifetime eggs", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

############ EFD 1st gonotrophic cycle
# Get data
data.specific <- with(data.block.temp, data.frame('T' = temp, 'trait' = EFD.1stGC)) # subset specific trait data from complete list
data <- data.specific # assign trait data to variable 'data'
trait.prior.fits = prior.fits.list$EFD_posterior1[,c('n.qd', 'T0', 'Tm')]
hypers <- trait.prior.fits * 0.2 # assign priors to variable 'priors' and weight importance if desired

# Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS - **select correct model file**
model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_T_inf.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model.out$BUGSoutput$summary[1:10,]
mcmcplot(model.out)

# Save model output 
save(model.out, file = "saved posteriors inf/EFD1stGC_briere_inf.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
     ylab = "EFD first GC", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

############ Estimated biting rate from gonotrophic cycle
# Get data
#remove 0's from dataset
data.no0 <- data.all[!(data.all$est.bite == 0),]
est.bite.means <- data.no0 %>%
               dplyr::group_by(temp,block) %>%
               dplyr::summarise(est.bite = mean(est.bite)) %>%
               ungroup()

est.bite.means <- as.data.frame(est.bite.means)
data.specific <- with(est.bite.means, data.frame('T' = temp, 'trait' = est.bite)) # subset specific trait data from complete list
data <- data.specific # assign trait data to variable 'data'
trait.prior.fits = prior.fits.list$a_posterior1[,c('c', 'T0', 'Tm')]
hypers <- trait.prior.fits * 0.5 # assign priors to variable 'priors' and weight importance if desired

# Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS - **select correct model file**
model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_T_inf.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model.out$BUGSoutput$summary[1:10,]
mcmcplot(model.out)

# Save model output 
save(model.out, file = "saved posteriors inf/est.bite_briere_inf.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
     ylab = "Estimated biting rate", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

############ EIP50 from Shapiro et al. 2017 Plos Biology
# Get data
data.specific <- with(data.bc.EIP, data.frame('T' = temp, 'trait' = 1/EIP50)) # subset specific trait data from complete list
data <- data.specific # assign trait data to variable 'data'
trait.prior.fits = prior.fits.list$PDR_post_prior1[,c('c', 'T0', 'Tm')]
hypers <- trait.prior.fits * 0.008 # assign priors to variable 'priors' and weight importance if desired
#changed to 0.03 as priors were having too strong of an influence on the current data

# Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS - **select correct model file**
model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_T_inf.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model.out$BUGSoutput$summary[1:10,]
mcmcplot(model.out)

# Save model output 
save(model.out, file = "saved posteriors inf/shapiro_PDR_briere_inf.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), ylim = c(0, 0.2), data = data.specific, 
     ylab = "PDR", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

############ Vector competence from Shapiro et al. 2017 Plos Biology
# Get data
data.specific <- with(data.bc.EIP, data.frame('T' = temp, 'trait' = bc)) # subset specific trait data from complete list
data <- data.specific # assign trait data to variable 'data'
trait.prior.fits = prior.fits.list$bc_posterior1_quad[,c('n.qd', 'T0', 'Tm')]
hypers <- trait.prior.fits * 0.1 # assign priors to variable 'priors' and weight importance if desired

# Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS - **select correct model file**
model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_T_inf.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model.out$BUGSoutput$summary[1:10,]
mcmcplot(model.out)

# Save model output 
save(model.out, file = "saved posteriors inf/shapiro_bc_quadratic_inf.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), ylim = c(0,0.7), data = data.specific, 
     ylab = "Vector competence", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


########### Lifespan
# Get data
data.specific <- with(data.block.temp, data.frame('T' = temp, 'trait' = lifespan)) # subset specific trait data from complete list
data <- data.specific # assign trait data to variable 'data'
trait.prior.fits = prior.fits.list$lf_posterior1[,c('q', 'T0', 'Tm')]
hypers <- trait.prior.fits * 0.5 # assign priors to variable 'priors' and weight importance if desired

# Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS - **select correct model file**
model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_T_inf.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model.out$BUGSoutput$summary[1:10,]
mcmcplot(model.out)

# Save model output 
save(model.out, file = "saved posteriors inf/lifespan_quadratic_inf.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
     ylab = "Lifespan", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

############ Estimated lifespan
# Get data
data.specific <- with(data.mu, data.frame('T' = temp, 'trait' = inverse.mu)) # subset specific trait data from complete list
data <- data.specific # assign trait data to variable 'data'
trait.prior.fits = prior.fits.list$lf_posterior1[,c('q', 'T0', 'Tm')]
hypers <- trait.prior.fits * 0.5 # assign priors to variable 'priors' and weight importance if desired

# Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS - **select correct model file**
model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_T_inf.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model.out$BUGSoutput$summary[1:10,]
mcmcplot(model.out)

# Save model output 
save(model.out, file = "saved posteriors inf/mu_exponential_quadratic_inf.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
     ylab = "1/mu", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
############ Mosqutio developement rate from Paaijmans et al. 2013 Global Climate - Briere
# Get data
data.specific <- with(data.pea.MDR, data.frame('T' = temp, 'trait' = MDR)) # subset specific trait data from complete list
data <- data.specific # assign trait data to variable 'data'
trait.prior.fits = prior.fits.list$MDR_posterior1[,c('c', 'T0', 'Tm')]
hypers <- trait.prior.fits * 0.5 # assign priors to variable 'priors' and weight importance if desired

# Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS - **select correct model file**
model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_T_inf.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model.out$BUGSoutput$summary[1:10,]
mcmcplot(model.out)

# Save model output
save(model.out, file = "saved posteriors inf/MDRK_briere_inf.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, xaxt = "n",
     ylab = "Mosquito development rate", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

########## Prob. egg to adult survival from Paaijmans et al. 2013 Global Climate - Quadratic
# Get data
data.specific <- with(data.pea.MDR, data.frame('T' = temp, 'trait' = Pea)) # subset specific trait data from complete list
data <- data.specific # assign trait data to variable 'data'
trait.prior.fits = prior.fits.list$e2a_posterior1[,c('n.qd', 'T0', 'Tm')]
hypers <- trait.prior.fits * 0.5 # assign priors to variable 'priors' and weight importance if desired

# Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS - **select correct model file**
model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_T_inf.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model.out$BUGSoutput$summary[1:10,]
mcmcplot(model.out)

# Save model output 
save(model.out, file = "saved posteriors inf/peaK_quadratic_inf.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), ylim = c(0,1), data = data.specific, 
     ylab = "Prob. e2a", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


