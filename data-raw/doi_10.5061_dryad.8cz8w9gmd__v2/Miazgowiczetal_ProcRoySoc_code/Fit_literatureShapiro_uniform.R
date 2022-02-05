## Marta Shocket, Stanford University
## Updated by Erin Mordecai, August 6, 2018
## Updated by Kerri Miazgowicz between August 6, 2018 and August 5, 2019
### Modified by KM on April 12, 2020 to compare model fits using
##            i. the full individual dataset vs. block|temp means
##            ii. normal distribution vs. truncated normal distributions
## Purpose: Use Bayesian Inference (JAGS) to fit temperature-dependent functions traits
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) JAGS models
##           3) Shared settings for all models
##           4) Code for fitting a trait - uniform priors
##           5) Example code for fitting a trait - informative priors fit from data
##           6) Example code for fitting prior distribution



##########
###### 1. Set up workspace, load packages, get data, etc.
##########

mainDir = "C:/Users/Kerri/Desktop/Chapter1 Submission"
setwd(mainDir)

# Check whether there's a folder in the directory for saving plots
# If not, create one
subDir = "saved posteriors"
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)

# Load libraties for fitting traits
library('R2jags') # Fits Bayesian models
library('mcmcplots') # Diagnostic plots for fits
library('MASS') # Fits distributions for informative priors
library(plyr) # Slices and dices data
library(plotrix) # For standard error function

#change default plot specifications
par(mar=c(1,1,1,1))

# Load raw trait data - for the code below to work, trait value should be called 'trait' & temperature should be called 'T'
data.bc.EIP <- read.csv("data/forErin_ShapiroData.csv")
data.bc.EIP$inverse.EIP50 = 1/data.bc.EIP$EIP50
data.pea.MDR <- read.csv("data/Krijn_Raw_Data.csv")
#replace Na with 0s
data.pea.MDR[is.na(data.pea.MDR)] = 0

##########
###### 2. JAGS Models
##########

# NOTES:
# - Running the code below writes a .txt file for each model to your current working directory.
# - The models include a section for 'derived quanities' that calculates the trait across a
#     temperature gradient for each saved sample in the MCMC chain. This output is what we
#     use later on to calculate R0. 

############## Quadratic Model with uniform priors

sink("quad_T.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 24)
    cf.Tm ~ dunif(25, 45)
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

############## Quadratic Model with uniform priors - derived quantities always =< 1 (for traits that are proportions)

sink("quadprob_T.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 24)
    cf.Tm ~ dunif(25, 45)
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)T(0,) 
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])) * (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) < 1) + (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) > 1)
    }
    
    } # close model
    ",fill=T)
sink()

############## Briere Model with uniform priors

sink("briere_T.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 20)
    cf.Tm ~ dunif(28, 45)
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

###########Briere Model with extended contraints on the Tmax
sink("Mod.briere_T.txt") #"Mod.briere_T.txt"
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 20)
    cf.Tm ~ dunif(28, 55)
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


#################3
############## Linear Model with uniform priors and no truncation

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

##### Derived Quantity Settings
Temp.xs <- seq(0, 45, 0.2) # temperature gradient to calculate derived quantities over
N.Temp.xs <-length(Temp.xs)

# save the temperature sequence for future analyses
save(Temp.xs, file = "saved posteriors/temps.Rdata")

##########
###### 4. Fitting traits - uniform priors
##########

############## First plot each trait to check functional forms
#create a new column for individual EFD.1stGC*lifespan to get an individual-based estimated lifetime egg production
boxplot(bc ~ temp, data = data.bc.EIP) # quad
boxplot(inverse.EIP50 ~ temp, data = data.bc.EIP) # briere


############## Trait name and parameter symbol - model type (e.g. quadratic, Briere)

############ Vector competence bc - no truncation
# Get data
data.specific <- with(data.bc.EIP, data.frame('T' = temp, 'trait' = bc)) # subset specific trait data from complete list
data <- data.specific # assign trait data to variable 'data'

# Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                      n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model.out$BUGSoutput$summary[1:10,]
mcmcplot(model.out)
DIC.bc.briere.noT <- model.out$BUGSoutput$DIC

# Save model output 
save(model.out, file = "saved posteriors/shapiro_bc_briere_noT_uniform.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
     ylab = "bc", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#################################

#############same data but now the truncated briere
# Run JAGS - **select correct model file**
model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_T.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model.out$BUGSoutput$summary[1:10,]
mcmcplot(model.out)
DIC.bc.briere.withT <- model.out$BUGSoutput$DIC

# Save model output 
save(model.out, file = "saved posteriors/shapiro_bc_briere_withT_uniform.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
     ylab = "bc", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

#############same data for shapiro bc but now the quadratic no T - only save the DIC value for comparison
# Run JAGS - **select correct model file**
model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model.out$BUGSoutput$summary[1:10,]
DIC.bc.quad.noT <- model.out$BUGSoutput$DIC

# Save model output 
save(model.out, file = "saved posteriors/shapiro_bc_quad_noT_uniform.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
     ylab = "bc", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

#############same data for shapiro bc but now the quadratic with Truncation - only save the DIC value for comparison
# Run JAGS - **select correct model file**
model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_T.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model.out$BUGSoutput$summary[1:10,]
DIC.bc.quad.withT <- model.out$BUGSoutput$DIC

# Save model output 
save(model.out, file = "saved posteriors/shapiro_bc_quad_withT_uniform.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
     ylab = "bc", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


#############################################################################
#############
#############################################################################
######### Pathogen developement rate
###############
#####################################################################
############ PDR - no truncation
# Get data
data.specific <- with(data.bc.EIP, data.frame('T' = temp, 'trait' = inverse.EIP50)) # subset specific trait data from complete list
data <- data.specific # assign trait data to variable 'data'

# Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model.out$BUGSoutput$summary[1:10,]
mcmcplot(model.out)
DIC.PDR.briere.noT <- model.out$BUGSoutput$DIC

# Save model output 
save(model.out, file = "saved posteriors/shapiro_PDR_briere_noT_uniform.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
     ylab = "PDR", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#################################

#############same data but now the truncated briere
# Run JAGS - **select correct model file**
model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_T.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model.out$BUGSoutput$summary[1:10,]
mcmcplot(model.out)
DIC.PDR.briere.withT <- model.out$BUGSoutput$DIC

# Save model output 
save(model.out, file = "saved posteriors/shapiro_PDR_briere_withT_uniform.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
     ylab = "PDR", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

#############same data for shapiro bc but now the quadratic no T - only save the DIC value for comparison
# Run JAGS - **select correct model file**
model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model.out$BUGSoutput$summary[1:10,]
DIC.PDR.quad.noT <- model.out$BUGSoutput$DIC

# Save model output 
save(model.out, file = "saved posteriors/shapiro_PDR_quad_noT_uniform.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
     ylab = "PDR", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

#############same data for shapiro PDR but now the quadratic with Truncation - only save the DIC value for comparison
# Run JAGS - **select correct model file**
model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_T.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model.out$BUGSoutput$summary[1:10,]
DIC.PDR.quad.withT <- model.out$BUGSoutput$DIC

# Save model output 
save(model.out, file = "saved posteriors/shapiro_PDR_quad_withT_uniform.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
     ylab = "PDR", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
