## Marta Shocket, Stanford University
## Updated Oct 2017
##
## Purpose: Use Bayesian Inference (JAGS) to fit temperature-dependent functions traits
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) JAGS models
##           3) Shared settings for all models
##           4) Example code for fitting a trait - uniform priors
##           5) Example code for fitting a trait - informative priors fit from data
##           6) Example code for fitting prior distribution



##########
###### 1. Set up workspace, load packages, get data, etc.
##########

# Set working directory
setwd("~/Dropbox/Articles/Vector-borne Disease Temp Project/Fitting Traits")

# Load libraties for fitting traits
library('R2jags') # Fits Bayesian models
library('mcmcplots') # Diagnostic plots for fits
library('MASS') # Fits distributions for informative priors

# Load raw trait data - for the code below to work, trait value should be called 'trait' & temperature should be called 'T'
data.all <- read.csv("file.csv")

# Load prior fits if needed
load("PriorFitsList.Rdata")


##########
###### 2. JAGS Models
##########

# NOTES:
# - Running the code below writes a .txt file for each model to your current working directory.
# - The models include a section for 'derived quanities' that calculates the trait across a
#     temperature gradient for each saved sample in the MCMC chain. This output is what we
#     use later on to calculate R0. 

############## Quadratic Model with uniform priors

sink("quad.txt")
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
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])
    }
    
    } # close model
    ",fill=T)
sink()

############## Quadratic Model with uniform priors - derived quantities always =< 1 (for traits that are proportions)

sink("quadprob.txt")
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
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])) * (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) < 1) + (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) > 1)
    }
    
    } # close model
    ",fill=T)
sink()

############## Briere Model with uniform priors

sink("briere.txt")
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
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i])
    }
    
    } # close model
    ",fill=T)
sink()

############## Quadratic Model with gamma priors

sink("quad_inf.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dgamma(hypers[1,1], hypers[2,1])
    cf.T0 ~ dgamma(hypers[1,2], hypers[2,2])
    cf.Tm ~ dgamma(hypers[1,3], hypers[2,3])
    cf.sigma ~ dgamma(hypers[1,4], hypers[2,4])
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])
    }
    
    } # close model
    ",fill=T)
sink()

############## Briere Model with gamma priors

sink("briere_inf.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dgamma(hypers[1,1], hypers[2,1])
    cf.T0 ~ dgamma(hypers[1,2], hypers[2,2])
    cf.Tm ~ dgamma(hypers[1,3], hypers[2,3])
    cf.sigma ~ dgamma(hypers[1,4], hypers[2,4])
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- cf.q * temp[i] * (temp[i] - cf.T0) * sqrt((cf.Tm - temp[i]) * (cf.Tm > temp[i])) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
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

##### Derived Quantity Settings
Temp.xs <- seq(5, 45, 0.2) # temperature gradient to calculate derived quantities over
N.Temp.xs <-length(Temp.xs)


##########
###### 4. Example code for fitting a trait - uniform priors
##########

############## Trait name and parameter symbol - model type (e.g. quadratic, Briere)

# Get data
data.specific <- subset(data.all, trait == "traitname") # subset specific trait data from complete list
data <- data.specific # assign trait data to variable 'data'

# Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                      n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model.out$BUGSoutput$summary[1:10,]
mcmcplot(model.out)

# Save model output 
save(model.out, file = "model.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,10), data = data.specific, xaxt = "n",
     ylab = "Trait name", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##########
###### 5. Example code for fitting a trait - informative priors fit from data
##########

############## Trait name and parameter symbol - model type (e.g. quadratic, Briere)

# Get data
data.specific <- subset(data.all, trait == "traitname") # subset specific trait data from complete list
data <- data.specific # assign trait data to variable 'data'
hypers <- trait.prior.fits * 0.5 # assign priors to variable 'priors' and weight importance if desired

# Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS - **select correct model file**
model.inf.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model.inf.out$BUGSoutput$summary[1:10,]
mcmcplot(model.inf.out)

# Save model output 
save(model.inf.out, file = "model_inf.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,10), data = data.specific, xaxt = "n",
     ylab = "Trait name", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.inf.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.inf.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.inf.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##########
###### 6. Example code for fitting prior distributions
##########

# Load Prior Fits
load("model.Rdata")

# Get the posterior dists for all 4 parameters into a data frame - **replace 'model.out' with correct name**
prior.dists <- data.frame(q = as.vector(model.out$BUGSoutput$sims.list$cf.q),
                                T0 = as.vector(model.out$BUGSoutput$sims.list$cf.T0), 
                                Tm = as.vector(model.out$BUGSoutput$sims.list$cf.Tm), 
                                sigma = as.vector(model.out$BUGSoutput$sims.list$cf.sigma))

# Fit gamma distributions for each parameter posterior dists  - **replace 'prior.dists' with correct name**
prior.fits = apply(prior.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

# Save prior dist fits either separately or multiple in a list
prior.fits.list <- list(prior.fit.1, prior.fit.2)
save(prior.fits.list, file = "PriorFitsList.Rdata")