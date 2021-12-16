## Kerri Miazgowicz, University of Georgia
## Modified from code provided by Erin Mordecai on August 6, 2018
## October 1, 2018
## 
## Purpose: Use Bayesian Inference (JAGS) to fit temperature-dependent functions traits over constant and fluctuating temperatures
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

mainDir = "C:/Users/Kerri/Desktop/Fluctuation_BayesianFits"
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
#Change these to the appropiate files
#Data collected from Miazgowicz for lf, a, and B are provided for each individual and generated in code file DataFormatting.R

data.constant <- read.csv("data-raw/constant.individual.trait.csv")
names(data.constant)[names(data.constant) == 'Treatment'] <- 'temp'
data.fluc9 <- read.csv("data-raw/fluc9.individual.trait.csv")
names(data.fluc9)[names(data.fluc9) == 'Treatment'] <- 'temp'
data.fluc12 <- read.csv("data-raw/fluc12.individual.trait.csv")
names(data.fluc12)[names(data.fluc12) == 'Treatment'] <- 'temp'

# data.bc.EIP <- read.csv("data/forErin_ShapiroData.csv")
# data.bc.EIP$inverse.EIP50 = 1/data.bc.EIP$EIP50
# data.pea.MDR <- read.csv("data/Krijn_Raw_Data.csv")

#Note: We kept the data at 36C at constant temp for the TPC fits 


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
Temp.xs <- seq(0, 45, 0.1) # temperature gradient to calculate derived quantities over
N.Temp.xs <-length(Temp.xs)

# save the temperature sequence for future analyses
save(Temp.xs, file = "saved-posteriors/temps.Rdata")

##########
###### 4. Fitting traits - uniform priors
##########

############## First plot each trait to check functional forms
boxplot(bite.rate ~ temp, data = data.constant) # briere
boxplot(bite.rate ~ temp, data = data.fluc9) # briere
boxplot(bite.rate ~ temp, data = data.fluc12) # briere /linear
boxplot(lifespan ~ temp, data = data.constant) # quad
boxplot(lifespan ~ temp, data = data.fluc9) # quad
boxplot(lifespan ~ temp, data = data.fluc12) #quad
boxplot(lifetime.eggs ~ temp, data = data.constant) # quad
boxplot(lifetime.eggs ~ temp, data = data.fluc9) #quad
boxplot(lifetime.eggs ~ temp, data = data.fluc12) #quad

plot(bc ~ temp, data = data.bc.EIP, pch = 16) # quad?
plot(1/EIP50 ~ temp, data = data.bc.EIP, pch = 16) # briere / linear
plot(Pea~ temp, data = data.pea.MDR) #quadratic 
plot(MDR ~temp, data = data.pea.MDR) #briere

############## Trait name and parameter symbol - model type (e.g. quadratic, Briere)

############ Bite rate at constant temperature
# Get data
data.specific <- with(data.constant, data.frame('T' = temp, 'trait' = bite.rate)) # subset specific trait data from complete list
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

# Save model output 
save(model.out, file = "saved-posteriors/constant_bite.rate_briere_uniform.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
     ylab = "Bite rate", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


############ Bite rate at DTR 9C
# Get data
data.specific <- with(data.fluc9, data.frame('T' = temp, 'trait' = bite.rate)) # subset specific trait data from complete list
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

# Save model output 
save(model.out, file = "saved-posteriors/dtr9_bite.rate_briere_uniform.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
     ylab = "Bite rate", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

############ Bite rate at DTR 12C
# Get data
data.specific <- with(data.fluc12, data.frame('T' = temp, 'trait' = bite.rate)) # subset specific trait data from complete list
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

# Save model output 
save(model.out, file = "saved-posteriors/dtr12_bite.rate_briere_uniform.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
     ylab = "Bite rate", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



############ Lifespan at constant temperature
# Get data
data.specific <- with(data.constant, data.frame('T' = temp, 'trait' = lifespan)) # subset specific trait data from complete list
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
save(model.out, file = "saved-posteriors/constant_lifespan_quadratic_uniform.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
     ylab = "Lifespan", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

############ Lifespan at DTR 9C
# Get data
data.specific <- with(data.fluc9, data.frame('T' = temp, 'trait' = lifespan)) # subset specific trait data from complete list
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
save(model.out, file = "saved posteriors/dtr9_lifespan_quadratic_uniform.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
     ylab = "Lifespan", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

############ Lifespan at DTR 12C
# Get data
data.specific <- with(data.fluc12, data.frame('T' = temp, 'trait' = lifespan)) # subset specific trait data from complete list
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
save(model.out, file = "saved posteriors/dtr12_lifespan_quadratic_uniform.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
     ylab = "Lifespan", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


############ Lifetime eggs at constant temperature
# Get data
data.specific <- with(data.constant, data.frame('T' = temp, 'trait' = lifetime.eggs)) # subset specific trait data from complete list
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
save(model.out, file = "saved posteriors/constant_lifetime.eggs_quadratic_uniform.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
     ylab = "Lifetime eggs", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

############ Lifetime eggs at DTR 9C
# Get data
data.specific <- with(data.fluc9, data.frame('T' = temp, 'trait' = lifetime.eggs)) # subset specific trait data from complete list
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
save(model.out, file = "saved posteriors/dtr9_lifetime.eggs_quadratic_uniform.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
     ylab = "Lifetime eggs", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

############ Lifetime eggs at dtr12
# Get data
data.specific <- with(data.fluc12, data.frame('T' = temp, 'trait' = lifetime.eggs)) # subset specific trait data from complete list
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
save(model.out, file = "saved posteriors/dtr12_lifetime.eggs_quadratic_uniform.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
     ylab = "Lifetime eggs", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


############ EIP50 from Shapiro et al. 2017 Plos Biology
# Get data
data.specific <- with(data.bc.EIP, data.frame('T' = temp, 'trait' = 1/EIP50)) # subset specific trait data from complete list
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

# Save model output 
save(model.out, file = "saved posteriors/EIP50_briere_uniform.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
     ylab = "EIP-50", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

############ Vector competence from Shapiro et al. 2017 Plos Biology
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
model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model.out$BUGSoutput$summary[1:10,]
mcmcplot(model.out)

# Save model output 
save(model.out, file = "saved posteriors/bc_quadratic_uniform.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
     ylab = "Vector competence", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


############ Pea - Quadratic from Paaijmans 2013 Global Climate Change
# Get data
data.specific <- with(data.pea.MDR, data.frame('T' = temp, 'trait' = Pea)) # subset specific trait data from complete list
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
save(model.out, file = "saved posteriors/pea_quadratic_uniform_Krijn2013.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
     ylab = "Pea", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

############ Mosqutio development rate from Paaijmans 2013 Global Climate Change
# Get data
data.specific <- with(data.pea.MDR, data.frame('T' = temp, 'trait' = MDR)) # subset specific trait data from complete list
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

# Save model output 
save(model.out, file = "saved posteriors/MDR_briere_uniform_Krijn2013.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
     ylab = "MDR", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


#####################
#####Load Johnson fits, use Johnson fits for EIP, bc, MDR, and pEA as informative priors
####################

# Load prior fits
load("Johnson et al posteriors/PriorFitsList.Rdata")



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



############ EIP50 from Shapiro et al. 2017 Plos Biology
# Get data
data.specific <- with(data.bc.EIP, data.frame('T' = temp, 'trait' = 1/EIP50)) # subset specific trait data from complete list
data <- data.specific # assign trait data to variable 'data'
trait.prior.fits = prior.fits.list$PDR_post_prior1[,c('c', 'T0', 'Tm')]
hypers <- trait.prior.fits * 0.1 # assign priors to variable 'priors' and weight importance if desired

# Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS - **select correct model file**
model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model.out$BUGSoutput$summary[1:10,]
mcmcplot(model.out)

# Save model output 
save(model.out, file = "saved posteriors/EIP50_briere_inf.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), ylim = c(0, 0.2), data = data.specific, 
     ylab = "EIP-50", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


############ Vector competence from Shapiro et al. 2017 Plos Biology
# Get data
data.specific <- with(data.bc.EIP, data.frame('T' = temp, 'trait' = bc)) # subset specific trait data from complete list
data <- data.specific # assign trait data to variable 'data'
trait.prior.fits = prior.fits.list$bc_posterior1_quad[,c('n.qd', 'T0', 'Tm')]
hypers <- trait.prior.fits * 0.5 # assign priors to variable 'priors' and weight importance if desired

# Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS - **select correct model file**
model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model.out$BUGSoutput$summary[1:10,]
mcmcplot(model.out)

# Save model output 
save(model.out, file = "saved posteriors/bc_quadratic_inf.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), ylim = c(0,0.7), data = data.specific, 
     ylab = "Vector competence", xlab = expression(paste("Temperature (",degree,"C)")))
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
model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model.out$BUGSoutput$summary[1:10,]
mcmcplot(model.out)

# Save model output
save(model.out, file = "saved posteriors/MDR_briere_inf.Rdata")

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
model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model.out$BUGSoutput$summary[1:10,]
mcmcplot(model.out)

# Save model output 
save(model.out, file = "saved posteriors/pea_quadratic_inf.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), ylim = c(0,0.7), data = data.specific, 
     ylab = "Prob. e2a", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


