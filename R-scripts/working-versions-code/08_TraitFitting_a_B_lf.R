## Rate Summation Project
## Script for fitting TPCs using JAGS
## Written by Joey Bernhardt & Marta Shocket
## Modified from code provided by Kerri Miazgowicz in 2020, which was modified from code provided by Marta Shocket in 2018


##	1. Prepare workspace
##	2. JAGS models
##	3. Shared settings for all models
##	4. Fitting traits - individual data from fluctuating temperature experiment: biting rate (a), lifespan (lf), & lifetime eggs (B)
##	5. Fitting traits - group mean data from previous experiments: vector competence (bc), extrinsic incubation period (EIP), 
##		juvenile survival (pEA), and mosquito development rate (MDR)
##	6. Fitting traits - gamma - block mean data - calculated by combining our data with data from previous experiments
##	7. Process and save model output for plotting




##########
###### 1. Prepare work space, get data, etc. 
##########

##### Load libraties for fitting traits
library(tidyverse)
library('R2jags') # Fits Bayesian models
library('mcmcplots') # Diagnostic plots for fits
#library('MASS') # Fits distributions for informative priors
#library(plotrix) # For standard error function

##### Load raw trait data from fluctuating temperature experiment
data.constant <- read.csv("data-raw/constant.individual.trait.csv")
names(data.constant)[names(data.constant) == 'Treatment'] <- 'temp'
data.fluc9 <- read.csv("data-raw/fluc9.individual.trait.csv")
names(data.fluc9)[names(data.fluc9) == 'Treatment'] <- 'temp'
data.fluc12 <- read.csv("data-raw/fluc12.individual.trait.csv")
names(data.fluc12)[names(data.fluc12) == 'Treatment'] <- 'temp'
# Note: We kept the data at 36C at constant temp for the TPC fits 

##### Load data from previosu studies for other traits
data.pea.MDR <- read.csv("data-raw/Krijn_Raw_Data.csv")
data.bc.EIP <- read.csv("data-raw/forErin_ShapiroData.csv")
data.bc.EIP$inverse.EIP50 = 1/data.bc.EIP$EIP50

##### Load gamma calculations
data.gamma.constant <- read.csv("data-processed/gamma_values_constant.csv")
data.gamma.dtr9 <- read.csv("data-processed/gamma_values_dtr9.csv")
data.gamma.dtr12 <- read.csv("data-processed/gamma_values_dtr12.csv")

# const_data <- read.csv("data-raw/kerri-data-from-drive/constant_master.csv")
# fluc9_data <- read.csv("data-raw/kerri-data-from-drive/fluc_dtr9_master.csv")
# fluc12_data <- read.csv("data-raw/kerri-data-from-drive/fluc_dtr12_master.csv")
# 
# const_surv_KM_values <- read.csv("data-raw/kerri-data-from-drive/KM_estimates_constant.csv")
# fluc9_surv_KM_values <- read.csv("data-raw/kerri-data-from-drive/KM_estimates_dtr9.csv")
# fluc12_surv_KM_values <- read.csv("data-raw/kerri-data-from-drive/KM_estimates_dtr12.csv")




##########
###### 2. JAGS Models
##########

# NOTES:
# - Running the code below writes a .txt file for each model to your current working directory.
# - The models include a section for 'derived quantities' that calculates the trait across a
#     temperature gradient for each saved sample in the MCMC chain. This output is what we
#     use later on to calculate R0. 


################### Truncated normal distributed Briere - for individual biting rate data

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


################### Gamma distributed Quadratic - for individual lifespan data

sink("quad_gamma.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 20)
    cf.Tm ~ dunif(28, 45)
    cf.ra ~ dunif(0, 100)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait[i] ~ dgamma(sh[i], cf.ra)
    sh[i] <- cf.ra * trait.mu[i]
    trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
    
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])
    }
    
    } # close model
    ",fill=T)
sink()


################### Negative Binomial distributed Briere - for individual lifetime eggs data

sink("briere_negbin.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 20)
    cf.Tm ~ dunif(28, 45)
    cf.r ~ dunif(0, 100)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait[i] ~ dnegbin(p[i], cf.r)
    p[i] <- cf.r / (cf.r + trait.mu[i])
    trait.mu[i] <- cf.q * temp[i] * (temp[i] - cf.T0) * sqrt((cf.Tm - temp[i]) * (cf.Tm > temp[i])) * (cf.T0 < temp[i])
    
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i])
    }
    
    } # close model
    ",fill=T)
sink()


# ############## Quadratic Model with uniform priors
# 
# sink("quad.txt")
# cat("
#     model{
#     
#     ## Priors
#     cf.q ~ dunif(0, 1)
#     cf.T0 ~ dunif(0, 24)
#     cf.Tm ~ dunif(25, 45)
#     cf.sigma ~ dunif(0, 1000)
#     cf.tau <- 1 / (cf.sigma * cf.sigma)
#     
#     ## Likelihood
#     for(i in 1:N.obs){
#     trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
#     trait[i] ~ dnorm(trait.mu[i], cf.tau)
#     }
#     
#     ## Derived Quantities and Predictions
#     for(i in 1:N.Temp.xs){
#     z.trait.mu.pred[i] <- -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])
#     }
#     
#     } # close model
#     ",fill=T)
# sink()
# 
# ############## Quadratic Model with uniform priors - derived quantities always =< 1 (for traits that are proportions)
# 
# sink("quadprob.txt")
# cat("
#     model{
#     
#     ## Priors
#     cf.q ~ dunif(0, 1)
#     cf.T0 ~ dunif(0, 24)
#     cf.Tm ~ dunif(25, 45)
#     cf.sigma ~ dunif(0, 1000)
#     cf.tau <- 1 / (cf.sigma * cf.sigma)
#     
#     ## Likelihood
#     for(i in 1:N.obs){
#     trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
#     trait[i] ~ dnorm(trait.mu[i], cf.tau)
#     }
#     
#     ## Derived Quantities and Predictions
#     for(i in 1:N.Temp.xs){
#     z.trait.mu.pred[i] <- (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])) * (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) < 1) + (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) > 1)
#     }
#     
#     } # close model
#     ",fill=T)
# sink()
# 
# ############## Briere Model with uniform priors
# 
# sink("briere.txt")
# cat("
#     model{
#     
#     ## Priors
#     cf.q ~ dunif(0, 1)
#     cf.T0 ~ dunif(0, 20)
#     cf.Tm ~ dunif(28, 45)
#     cf.sigma ~ dunif(0, 1000)
#     cf.tau <- 1 / (cf.sigma * cf.sigma)
#     
#     ## Likelihood
#     for(i in 1:N.obs){
#     trait.mu[i] <- cf.q * temp[i] * (temp[i] - cf.T0) * sqrt((cf.Tm - temp[i]) * (cf.Tm > temp[i])) * (cf.T0 < temp[i])
#     trait[i] ~ dnorm(trait.mu[i], cf.tau)
#     }
#     
#     ## Derived Quantities and Predictions
#     for(i in 1:N.Temp.xs){
#     z.trait.mu.pred[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i])
#     }
#     
#     } # close model
#     ",fill=T)
# sink()




##########
###### 3. Shared settings for all models
##########

##### MCMC Settings
# Number of posterior distribution elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 25000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains

##### Derived Quantity Settings
Temp.xs <- seq(0, 45, 0.1) # temperature gradient to calculate derived quantities over
N.Temp.xs <-length(Temp.xs)

# Save the temperature sequence for future analyses
# save(Temp.xs, file = "saved-posteriors/temps.Rdata")



##########
###### 4. Fitting traits - individual data from fluctuating temperature experiment: biting rate (a), lifespan (lf), & lifetime eggs (B)
##########


################################################ Model fitting: bite rate (a)


###################################### Settings for fitting normal / truncated normal distribution

#####  inits Function
inits<-function(){list(
	cf.q = 0.01,
	cf.Tm = 35,
	cf.T0 = 5,
	cf.sigma = rlnorm(1))}

#####  Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")


##### Bite rate at constant temperature ---------------------------------------

##### Get data
data_bite_rate_constant <- with(data.constant, data.frame('T' = temp, 'trait' = bite.rate)) # subset specific trait data from complete list

##### Organize Data for JAGS
trait <- data_bite_rate_constant$trait
N.obs <- length(trait)
temp <- data_bite_rate_constant$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS - **select correct model file**
model_bite_rate_constant <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_T.txt",
                      n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine model output & run diagnostics
model_bite_rate_constant$BUGSoutput$summary[1:5,]
# mcmcplot(model_bite_rate_constant)

# ##### Save model output 
# save(model_bite_rate_constant, file = "saved-posteriors/constant_biterate.Rdata")

##### Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_bite_rate_constant, 
     ylab = "Bite rate", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_bite_rate_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_bite_rate_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_bite_rate_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### Bite rate at DTR 9 ---------------------------------------

##### Get data
data_bite_rate_dtr9 <- with(data.fluc9, data.frame('T' = temp, 'trait' = bite.rate)) # subset specific trait data from complete list

##### Organize Data for JAGS
trait <- data_bite_rate_dtr9$trait
N.obs <- length(trait)
temp <- data_bite_rate_dtr9$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS - **select correct model file**
model_bite_rate_dtr9 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_T.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine model output & run diagnostics
model_bite_rate_dtr9$BUGSoutput$summary[1:5,]
# mcmcplot(model_bite_rate_dtr9)

# ##### Save model output 
# save(model_bite_rate_dtr9, file = "saved-posteriors/dtr9_biterate.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_bite_rate_dtr9,
     ylab = "Bite rate", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_bite_rate_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_bite_rate_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_bite_rate_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### Bite rate at DTR 12 ---------------------------------------

##### Get data
data_bite_rate_dtr12 <- with(data.fluc12, data.frame('T' = temp, 'trait' = bite.rate)) # subset specific trait data from complete list

##### Organize Data for JAGS
trait <- data_bite_rate_dtr12$trait
N.obs <- length(trait)
temp <- data_bite_rate_dtr12$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS - **select correct model file**
model_bite_rate_dtr12 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_T.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine model output & run diagnostics
model_bite_rate_dtr12$BUGSoutput$summary[1:5,]
# mcmcplot(model_bite_rate_dtr12)
 
# ##### Save model output
# save(model_bite_rate_dtr12, file = "saved-posteriors/dtr12_biterate.Rdata")

##### Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_bite_rate_dtr12,
     ylab = "Bite rate", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_bite_rate_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_bite_rate_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_bite_rate_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



################################################ Model fitting: lifespan (lf)


###################################### Settings for fitting gamma distribution

##### inits Function
inits <- function(){list(
	cf.q = 0.01,
	cf.Tm = 40,
	cf.T0 = 5,
	cf.ra = 0.003)}

##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm", "cf.ra", "z.trait.mu.pred")


##### Lifespan at constant temperature ---------------------------------------

##### Get data
data_lifespan_constant <- with(data.constant, data.frame('T' = temp, 'trait' = lifespan)) # subset specific trait data from complete list

##### Organize Data for JAGS
trait <- data_lifespan_constant$trait
N.obs <- length(trait)
temp <- data_lifespan_constant$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS - **select correct model file**
model_lifespan_constant <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_gamma.txt",
				  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine model output & run diagnostics
model_lifespan_constant$BUGSoutput$summary[1:5,]
# mcmcplot(model_lifespan_constant)

##### Save model output
# save(model_lifespan_constant, file = "saved-posteriors/constant_lifespan.Rdata")

##### Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_lifespan_constant,
     ylab = "Lifespan", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_lifespan_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_lifespan_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_lifespan_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### Lifespan at DTR 9 ---------------------------------------

##### Get data
data_lifespan_dtr9 <- with(data.fluc9, data.frame('T' = temp, 'trait' = lifespan)) # subset specific trait data from complete list

##### Organize Data for JAGS
trait <- data_lifespan_dtr9$trait
N.obs <- length(trait)
temp <- data_lifespan_dtr9$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS - **select correct model file**
model_lifespan_dtr9 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_gamma.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine model output & run diagnostics
model_lifespan_dtr9$BUGSoutput$summary[1:5,]
# mcmcplot(model_lifespan_dtr9)

##### Save model output
# save(model_lifespan_dtr9, file = "saved-posteriors/dtr9_lifespan.Rdata")

##### Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_lifespan_dtr9,
     ylab = "Lifespan", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_lifespan_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_lifespan_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_lifespan_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### Lifespan at DTR 12 ---------------------------------------

##### Get data
data_lifespan_dtr12 <- with(data.fluc12, data.frame('T' = temp, 'trait' = lifespan)) # subset specific trait data from complete list

##### Organize Data for JAGS
trait <- data_lifespan_dtr12$trait
N.obs <- length(trait)
temp <- data_lifespan_dtr12$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS - **select correct model file**
model_lifespan_dtr12 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_gamma.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine model output & run diagnostics
model_lifespan_dtr12$BUGSoutput$summary[1:5,]
# mcmcplot(model_lifespan_dtr12)

##### Save model output
# save(model_lifespan_dtr12, file = "saved-posteriors/dtr12_lifespan.Rdata")

##### Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_lifespan_dtr12,
     ylab = "Lifespan", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_lifespan_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_lifespan_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_lifespan_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


################################################ Model fitting: lifetime eggs (B)


###################################### Settings for fitting negative binomial distribution

##### inits Function
inits <- function(){list(
	cf.q = 0.01,
	cf.Tm = 40,
	cf.T0 = 5,
	cf.r = 1)}

##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm", "cf.r", "z.trait.mu.pred")

##### Lifetime eggs at constant -----------------------------------

##### Get data
data_eggs_constant <- with(data.constant, data.frame('T' = temp, 'trait' = lifetime.eggs)) # subset specific trait data from complete list

##### Organize Data for JAGS
trait <- data_eggs_constant$trait
N.obs <- length(trait)
temp <- data_eggs_constant$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS - **select correct model file**
model_eggs_constant <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_negbin.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine model output & run diagnostics
model_eggs_constant$BUGSoutput$summary[1:5,]
# mcmcplot(model_eggs_constant)

##### Save model output
# save(model_eggs_constant, file = "saved-posteriors/constant_eggs.Rdata")

##### Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_eggs_constant,
     ylab = "Lifetime eggs", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_eggs_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_eggs_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_eggs_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### Lifetime eggs at DTR 9 ------------------------------------

##### Get data
data_eggs_dtr9 <- with(data.fluc9, data.frame('T' = temp, 'trait' = lifetime.eggs)) # subset specific trait data from complete list

##### Organize Data for JAGS
trait <- data_eggs_dtr9$trait
N.obs <- length(trait)
temp <- data_eggs_dtr9$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS - **select correct model file**
model_eggs_dtr9 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_negbin.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model_eggs_dtr9$BUGSoutput$summary[1:5,]
# mcmcplot(model_eggs_dtr9)

##### Save model output
# save(model_eggs_dtr9, file = "saved-posteriors/dtr9_eggs.Rdata")

##### Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_eggs_dtr9,
     ylab = "Lifetime eggs", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_eggs_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_eggs_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_eggs_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### Lifetime eggs at DTR 12 ------------------------------------

##### Get data
data_eggs_dtr12 <- with(data.fluc12, data.frame('T' = temp, 'trait' = lifetime.eggs)) # subset specific trait data from complete list

##### Organize Data for JAGS
trait <- data_eggs_dtr12$trait
N.obs <- length(trait)
temp <- data_eggs_dtr12$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS - **select correct model file**
model_eggs_dtr12 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_negbin.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine model output & run diagnostics
model_eggs_dtr12$BUGSoutput$summary[1:5,]
# mcmcplot(model_eggs_dtr12)

##### Save model output
# save(model_eggs_dtr12, file = "saved-posteriors/dtr12_eggs_briere_uniform.Rdata")

##### Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_eggs_dtr12,
     ylab = "Lifetime eggs", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_eggs_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_eggs_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_eggs_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)




##########
###### 5. Fitting traits - group mean data from previous experiments: vector competence (bc), extrinic incubation period (EIP), juvenile survival (pEA), and mosquito development rate (MDR)
##########

#### The fluctating temperature experiment measured bite rate (a), adult mosquito lifespan (lf), lifetime egg production (B). 
# To calculate S(T) we need addtional parameters: vector competence (bc; the proportion of infectious mosquitoes),
#   parasite development rate or extrinsic incubation period (PDR = 1/EIP), probability of egg to adult survival (pEA), 
#   and mosquito development rate (MDR)


###################################### Settings for fitting normal / truncated normal distribution

#####  inits Function
inits<-function(){list(
	cf.q = 0.01,
	cf.Tm = 35,
	cf.T0 = 5,
	cf.sigma = rlnorm(1))}

#####  Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")



# Model fitting: EIP50 ----------------------------------------------------

############ EIP50 from Shapiro et al. 2017 Plos Biology
# Get data
data_eip50 <- with(data.bc.EIP, data.frame('T' = temp, 'trait' = 1/EIP50)) # subset specific trait data from complete list
# data <- data.specific # assign trait data to variable 'data'

# Organize Data for JAGS
trait <- data_eip50$trait
N.obs <- length(trait)
temp <- data_eip50$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
model_eip50 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model_eip50$BUGSoutput$summary[1:5,]
mcmcplot(model_eip50)

# Save model output
save(model_eip50, file = "saved-posteriors/EIP50.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), ylim = c(0,0.2), data = data_eip50,
     ylab = "EIP-50", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_eip50$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_eip50$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_eip50$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


# Model fitting: Vector competence ----------------------------------------

############ Vector competence from Shapiro et al. 2017 Plos Biology
# Get data
data_bc <- with(data.bc.EIP, data.frame('T' = temp, 'trait' = bc)) # subset specific trait data from complete list

# Organize Data for JAGS
trait <- data_bc$trait
N.obs <- length(trait)
temp <- data_bc$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
model_bc <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model_bc$BUGSoutput$summary[1:5,]
mcmcplot(model_bc)

# Save model output
save(model_bc, file = "saved-posteriors/bc.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), ylim = c(0,1), data = data_bc,
     ylab = "Vector competence", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_bc$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_bc$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_bc$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



# Model fitting: Pea ------------------------------------------------------

############ Pea - Quadratic from Paaijmans 2013 Global Climate Change
# Get data
data_pea <- with(data.pea.MDR, data.frame('T' = temp, 'trait' = Pea)) # subset specific trait data from complete list

# Organize Data for JAGS
trait <- data_pea$trait
N.obs <- length(trait)
temp <- data_pea$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
model_pea <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model_pea$BUGSoutput$summary[1:5,]
mcmcplot(model_pea)

# Save model output
save(model_pea, file = "saved-posteriors/pea.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_pea,
     ylab = "Pea", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_pea$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_pea$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_pea$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



# Model fitting: MDR ------------------------------------------------------

############ Mosqutio development rate from Paaijmans 2013 Global Climate Change
# Get data
data_mdr <- with(data.pea.MDR, data.frame('T' = temp, 'trait' = MDR)) # subset specific trait data from complete list
# data <- data.specific # assign trait data to variable 'data'

# Organize Data for JAGS
trait <- data_mdr$trait
N.obs <- length(trait)
temp <- data_mdr$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
model_mdr <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model_mdr$BUGSoutput$summary[1:5,]
mcmcplot(model_mdr)

# Save model output
save(model_mdr, file = "saved-posteriors/MDR.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_mdr,
     ylab = "MDR", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_mdr$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_mdr$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_mdr$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)




##########
###### 6. Fitting traits - gamma - block mean data - calculated by combining our data with data from previous experiments
##########


#### Gamma - the proportion of mosquitoes surviving the EIP50 - is calculated based on Kaplan-Meier survival 
# estimates from our experiment (with constant and fluctuating temperatures) and the EIP50 values estimated  
# for our temperature treatments based on the TPC fit to data from previous studies (above).


################################################ Model fitting: Gamma


###################################### Settings for fitting normal / truncated normal distribution

#####  inits Function
inits<-function(){list(
	cf.q = 0.01,
	cf.Tm = 35,
	cf.T0 = 5,
	cf.sigma = rlnorm(1))}

#####  Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")


##### Gamma at constant -----------------------------------

# Get data
data_gamma_constant <- with(data.gamma.constant, data.frame('T' = temp, 'trait' = gamma)) # subset specific trait data from complete list

# Organize Data for JAGS
trait <- data_gamma_constant$trait
N.obs <- length(trait)
temp <- data_gamma_constant$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
model_gamma_constant <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
					n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model_gamma_constant$BUGSoutput$summary[1:5,]
mcmcplot(model_gamma_constant)

# Save model output
save(model_gamma_constant, file = "saved-posteriors/constant_gamma.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), ylim = c(0,1), data = data_gamma_constant,
	 ylab = "Gamma", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_gamma_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_gamma_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_gamma_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### Gamma at dtr9 -----------------------------------

# Get data
data_gamma_dtr9 <- with(data.gamma.dtr9, data.frame('T' = temp, 'trait' = gamma)) # subset specific trait data from complete list

# Organize Data for JAGS
trait <- data_gamma_dtr9$trait
N.obs <- length(trait)
temp <- data_gamma_dtr9$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
model_gamma_dtr9 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
							 n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model_gamma_dtr9$BUGSoutput$summary[1:5,]
mcmcplot(model_gamma_dtr9)

# Save model output
save(model_gamma_dtr9, file = "saved-posteriors/dtr9_gamma.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), ylim = c(0,1), data = data_gamma_dtr9,
	 ylab = "Gamma", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_gamma_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_gamma_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_gamma_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### Gamma at dtr12 -----------------------------------

# Get data
data_gamma_dtr12 <- with(data.gamma.dtr12, data.frame('T' = temp, 'trait' = gamma)) # subset specific trait data from complete list

# Organize Data for JAGS
trait <- data_gamma_dtr12$trait
N.obs <- length(trait)
temp <- data_gamma_dtr12$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
model_gamma_dtr12 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
							 n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model_gamma_dtr12$BUGSoutput$summary[1:5,]
mcmcplot(model_gamma_dtr12)

# Save model output
save(model_gamma_dtr12, file = "saved-posteriors/dtr12_gamma.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), ylim = c(0,1), data = data_gamma_dtr12,
	 ylab = "Gamma", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_gamma_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_gamma_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_gamma_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)










##########
###### 7. Process and save model output for plotting
##########



############## 1A. Bite rate: constant ---------------------------------------------------------------

##### Extract matrix of predicted values from fitted model (451 col = temperature gradient, 7500 rows = MCMC steps)
predictions_bite_rate_constant <- as.data.frame(model_bite_rate_constant$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
colnames(predictions_bite_rate_constant) <- Temp.xs

##### Create & save summary df with 2.5% and 97.5% CIs, mean, and median trait values for each temperature
predictions_bite_rate_constant_summary <- predictions_bite_rate_constant %>%
	mutate(iteration = rownames(.)) %>% 
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(temperature) %>%  
	summarise(lowerCI = quantile(trait_value, probs = 0.025),
			  upperCI = quantile(trait_value, probs = 0.975),
			  mean = mean(trait_value),
			  median = median(trait_value)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "bite_rate_constant")

# write_csv(predictions_bite_rate_constant_summary, "data-processed/predictions_bite_rate_constant_summary.csv")

##### Create summary df with 2.5% and 97.5% CIs, mean, and median Topt values for the TPC
topt_bite_rate_constant <- predictions_bite_rate_constant %>% 
	mutate(iteration = rownames(.)) %>%
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(iteration) %>% 
	slice_max(order_by = trait_value, n = 1) %>%
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(lowerCI = quantile(temperature, probs = 0.025),
			  upperCI = quantile(temperature, probs = 0.975),
			  mean = mean(temperature),
			  median = median(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "bite_rate_constant")

##### Create & save df with fitted parameters and Topt values from above
params_bite_rate_constant <- as.data.frame(model_bite_rate_constant$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>%
	rename(lowerCI = `2.5%`, upperCI = `97.5%`, median = `50%`) %>%
	mutate(treatment = "bite_rate_constant")

params_bite_rate_constant_all <- bind_rows(params_bite_rate_constant, topt_bite_rate_constant)

# write_csv(params_bite_rate_constant_all, "data-processed/params_bite_rate_constant_all.csv")

##### Create & save df with summarized data for plotting
data_bite_rate_constant_sum <- data_bite_rate_constant %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = sd(trait)/sqrt(n())) %>% 
	mutate(trait = "bite rate") %>% 
	mutate(treatment = "bite_rate_constant")

# write_csv(data_bite_rate_constant_sum, "data-processed/data_bite_rate_constant_sum.csv")



############## 1B. Bite rate: DTR 9 ---------------------------------------------------------------

##### Extract matrix of predicted values from fitted model (451 col = temperature gradient, 7500 rows = MCMC steps)
predictions_bite_rate_dtr9 <- as.data.frame(model_bite_rate_dtr9$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
colnames(predictions_bite_rate_dtr9) <- Temp.xs

##### Create & save summary df with 2.5% and 97.5% CIs, mean, and median trait values for each temperature
predictions_bite_rate_dtr9_summary <- predictions_bite_rate_dtr9 %>%
	mutate(iteration = rownames(.)) %>% 
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(temperature) %>%  
	summarise(lowerCI = quantile(trait_value, probs = 0.025),
			  upperCI = quantile(trait_value, probs = 0.975),
			  mean = mean(trait_value),
			  median = median(trait_value)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "bite_rate_dtr9")

# write_csv(predictions_bite_rate_dtr9_summary, "data-processed/predictions_bite_rate_dtr9_summary.csv")

##### Create summary df with 2.5% and 97.5% CIs, mean, and median Topt values for the TPC
topt_bite_rate_dtr9 <- predictions_bite_rate_dtr9 %>% 
	mutate(iteration = rownames(.)) %>%
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(iteration) %>% 
	slice_max(order_by = trait_value, n = 1) %>%
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(lowerCI = quantile(temperature, probs = 0.025),
			  upperCI = quantile(temperature, probs = 0.975),
			  mean = mean(temperature),
			  median = median(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "bite_rate_dtr9")

##### Create & save df with fitted parameters and Topt values from above
params_bite_rate_dtr9 <- as.data.frame(model_bite_rate_dtr9$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>%
	rename(lowerCI = `2.5%`, upperCI = `97.5%`, median = `50%`) %>%
	mutate(treatment = "bite_rate_dtr9")

params_bite_rate_dtr9_all <- bind_rows(params_bite_rate_dtr9, topt_bite_rate_dtr9)

# write_csv(params_bite_rate_dtr9_all, "data-processed/params_bite_rate_dtr9_all.csv")

##### Create & save df with summarized data for plotting
data_bite_rate_dtr9_sum <- data_bite_rate_dtr9 %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = sd(trait)/sqrt(n())) %>% 
	mutate(trait = "bite rate") %>% 
	mutate(treatment = "bite_rate_dtr9")

# write_csv(data_bite_rate_dtr9_sum, "data-processed/data_bite_rate_dtr9_sum.csv")



############## 1C. Bite rate: DTR 12 ---------------------------------------------------------------

##### Extract matrix of predicted values from fitted model (451 col = temperature gradient, 7500 rows = MCMC steps)
predictions_bite_rate_dtr12 <- as.data.frame(model_bite_rate_dtr12$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
colnames(predictions_bite_rate_dtr12) <- Temp.xs

##### Create & save summary df with 2.5% and 97.5% CIs, mean, and median trait values for each temperature
predictions_bite_rate_dtr12_summary <- predictions_bite_rate_dtr12 %>%
	mutate(iteration = rownames(.)) %>% 
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(temperature) %>%  
	summarise(lowerCI = quantile(trait_value, probs = 0.025),
			  upperCI = quantile(trait_value, probs = 0.975),
			  mean = mean(trait_value),
			  median = median(trait_value)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "bite_rate_dtr12")

# write_csv(predictions_bite_rate_dtr12_summary, "data-processed/predictions_bite_rate_dtr12_summary.csv")

##### Create summary df with 2.5% and 97.5% CIs, mean, and median Topt values for the TPC
topt_bite_rate_dtr12 <- predictions_bite_rate_dtr12 %>% 
	mutate(iteration = rownames(.)) %>%
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(iteration) %>% 
	slice_max(order_by = trait_value, n = 1) %>%
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(lowerCI = quantile(temperature, probs = 0.025),
			  upperCI = quantile(temperature, probs = 0.975),
			  mean = mean(temperature),
			  median = median(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "bite_rate_dtr12")

##### Create & save df with fitted parameters and Topt values from above
params_bite_rate_dtr12 <- as.data.frame(model_bite_rate_dtr12$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>%
	rename(lowerCI = `2.5%`, upperCI = `97.5%`, median = `50%`) %>%
	mutate(treatment = "bite_rate_dtr12")

params_bite_rate_dtr12_all <- bind_rows(params_bite_rate_dtr12, topt_bite_rate_dtr12)

# write_csv(params_bite_rate_dtr12_all, "data-processed/params_bite_rate_dtr12_all.csv")

##### Create & save df with summarized data for plotting
data_bite_rate_dtr12_sum <- data_bite_rate_dtr12 %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = sd(trait)/sqrt(n())) %>% 
	mutate(trait = "bite rate") %>% 
	mutate(treatment = "bite_rate_dtr12")

# write_csv(data_bite_rate_dtr12_sum, "data-processed/data_bite_rate_dtr12_sum.csv")



############## 2A. Lifespan: constant ---------------------------------------------------------------

##### Extract matrix of predicted values from fitted model (451 col = temperature gradient, 7500 rows = MCMC steps)
predictions_lifespan_constant <- as.data.frame(model_lifespan_constant$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
colnames(predictions_lifespan_constant) <- Temp.xs

##### Create & save summary df with 2.5% and 97.5% CIs, mean, and median trait values for each temperature
predictions_lifespan_constant_summary <- predictions_lifespan_constant %>%
	mutate(iteration = rownames(.)) %>% 
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(temperature) %>%  
	summarise(lowerCI = quantile(trait_value, probs = 0.025),
			  upperCI = quantile(trait_value, probs = 0.975),
			  mean = mean(trait_value),
			  median = median(trait_value)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "lifespan_constant")

# write_csv(predictions_lifespan_constant_summary, "data-processed/predictions_lifespan_constant_summary.csv")

##### Create summary df with 2.5% and 97.5% CIs, mean, and median Topt values for the TPC
topt_lifespan_constant <- predictions_lifespan_constant %>% 
	mutate(iteration = rownames(.)) %>%
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(iteration) %>% 
	slice_max(order_by = trait_value, n = 1) %>%
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(lowerCI = quantile(temperature, probs = 0.025),
			  upperCI = quantile(temperature, probs = 0.975),
			  mean = mean(temperature),
			  median = median(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "lifespan_constant")

##### Create & save df with fitted parameters and Topt values from above
params_lifespan_constant <- as.data.frame(model_lifespan_constant$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>%
	rename(lowerCI = `2.5%`, upperCI = `97.5%`, median = `50%`) %>%
	mutate(treatment = "lifespan_constant")

params_lifespan_constant_all <- bind_rows(params_lifespan_constant, topt_lifespan_constant)

# write_csv(params_lifespan_constant_all, "data-processed/params_lifespan_constant_all.csv")

##### Create & save df with summarized data for plotting
data_lifespan_constant_sum <- data_lifespan_constant %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = sd(trait)/sqrt(n())) %>% 
	mutate(trait = "bite rate") %>% 
	mutate(treatment = "lifespan_constant")

# write_csv(data_lifespan_constant_sum, "data-processed/data_lifespan_constant_sum.csv")



############## 2B. Lifespan: DTR 9 ---------------------------------------------------------------

##### Extract matrix of predicted values from fitted model (451 col = temperature gradient, 7500 rows = MCMC steps)
predictions_lifespan_dtr9 <- as.data.frame(model_lifespan_dtr9$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
colnames(predictions_lifespan_dtr9) <- Temp.xs

##### Create & save summary df with 2.5% and 97.5% CIs, mean, and median trait values for each temperature
predictions_lifespan_dtr9_summary <- predictions_lifespan_dtr9 %>%
	mutate(iteration = rownames(.)) %>% 
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(temperature) %>%  
	summarise(lowerCI = quantile(trait_value, probs = 0.025),
			  upperCI = quantile(trait_value, probs = 0.975),
			  mean = mean(trait_value),
			  median = median(trait_value)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "lifespan_dtr9")

# write_csv(predictions_lifespan_dtr9_summary, "data-processed/predictions_lifespan_dtr9_summary.csv")

##### Create summary df with 2.5% and 97.5% CIs, mean, and median Topt values for the TPC
topt_lifespan_dtr9 <- predictions_lifespan_dtr9 %>% 
	mutate(iteration = rownames(.)) %>%
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(iteration) %>% 
	slice_max(order_by = trait_value, n = 1) %>%
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(lowerCI = quantile(temperature, probs = 0.025),
			  upperCI = quantile(temperature, probs = 0.975),
			  mean = mean(temperature),
			  median = median(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "lifespan_dtr9")

##### Create & save df with fitted parameters and Topt values from above
params_lifespan_dtr9 <- as.data.frame(model_lifespan_dtr9$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>%
	rename(lowerCI = `2.5%`, upperCI = `97.5%`, median = `50%`) %>%
	mutate(treatment = "lifespan_dtr9")

params_lifespan_dtr9_all <- bind_rows(params_lifespan_dtr9, topt_lifespan_dtr9)

# write_csv(params_lifespan_dtr9_all, "data-processed/params_lifespan_dtr9_all.csv")

##### Create & save df with summarized data for plotting
data_lifespan_dtr9_sum <- data_lifespan_dtr9 %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = sd(trait)/sqrt(n())) %>% 
	mutate(trait = "bite rate") %>% 
	mutate(treatment = "lifespan_dtr9")

# write_csv(data_lifespan_dtr9_sum, "data-processed/data_lifespan_dtr9_sum.csv")



############## 2C. Lifespan: DTR 12 ---------------------------------------------------------------

##### Extract matrix of predicted values from fitted model (451 col = temperature gradient, 7500 rows = MCMC steps)
predictions_lifespan_dtr12 <- as.data.frame(model_lifespan_dtr12$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
colnames(predictions_lifespan_dtr12) <- Temp.xs

##### Create & save summary df with 2.5% and 97.5% CIs, mean, and median trait values for each temperature
predictions_lifespan_dtr12_summary <- predictions_lifespan_dtr12 %>%
	mutate(iteration = rownames(.)) %>% 
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(temperature) %>%  
	summarise(lowerCI = quantile(trait_value, probs = 0.025),
			  upperCI = quantile(trait_value, probs = 0.975),
			  mean = mean(trait_value),
			  median = median(trait_value)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "lifespan_dtr12")

# write_csv(predictions_lifespan_dtr12_summary, "data-processed/predictions_lifespan_dtr12_summary.csv")

##### Create summary df with 2.5% and 97.5% CIs, mean, and median Topt values for the TPC
topt_lifespan_dtr12 <- predictions_lifespan_dtr12 %>% 
	mutate(iteration = rownames(.)) %>%
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(iteration) %>% 
	slice_max(order_by = trait_value, n = 1) %>%
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(lowerCI = quantile(temperature, probs = 0.025),
			  upperCI = quantile(temperature, probs = 0.975),
			  mean = mean(temperature),
			  median = median(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "lifespan_dtr12")

##### Create & save df with fitted parameters and Topt values from above
params_lifespan_dtr12 <- as.data.frame(model_lifespan_dtr12$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>% 
	rename(lowerCI = `2.5%`, upperCI = `97.5%`, median = `50%`) %>%
	mutate(treatment = "lifespan_dtr12")

params_lifespan_dtr12_all <- bind_rows(params_lifespan_dtr12, topt_lifespan_dtr12)

# write_csv(params_lifespan_dtr12_all, "data-processed/params_lifespan_dtr12_all.csv")

##### Create & save df with summarized data for plotting
data_lifespan_dtr12_sum <- data_lifespan_dtr12 %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = sd(trait)/sqrt(n())) %>% 
	mutate(trait = "bite rate") %>% 
	mutate(treatment = "lifespan_dtr12")

# write_csv(data_lifespan_dtr12_sum, "data-processed/data_lifespan_dtr12_sum.csv")



############## 3A. Lifetime Eggs: constant ---------------------------------------------------------------

##### Extract matrix of predicted values from fitted model (451 col = temperature gradient, 7500 rows = MCMC steps)
predictions_eggs_constant <- as.data.frame(model_eggs_constant$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
colnames(predictions_eggs_constant) <- Temp.xs

##### Create & save summary df with 2.5% and 97.5% CIs, mean, and median trait values for each temperature
predictions_eggs_constant_summary <- predictions_eggs_constant %>%
	mutate(iteration = rownames(.)) %>% 
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(temperature) %>%  
	summarise(lowerCI = quantile(trait_value, probs = 0.025),
			  upperCI = quantile(trait_value, probs = 0.975),
			  mean = mean(trait_value),
			  median = median(trait_value)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "eggs_constant")

# write_csv(predictions_eggs_constant_summary, "data-processed/predictions_eggs_constant_summary.csv")

##### Create summary df with 2.5% and 97.5% CIs, mean, and median Topt values for the TPC
topt_eggs_constant <- predictions_eggs_constant %>% 
	mutate(iteration = rownames(.)) %>%
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(iteration) %>% 
	slice_max(order_by = trait_value, n = 1) %>%
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(lowerCI = quantile(temperature, probs = 0.025),
			  upperCI = quantile(temperature, probs = 0.975),
			  mean = mean(temperature),
			  median = median(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "eggs_constant")

##### Create & save df with fitted parameters and Topt values from above
params_eggs_constant <- as.data.frame(model_eggs_constant$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>%
	rename(lowerCI = `2.5%`, upperCI = `97.5%`, median = `50%`) %>%
	mutate(treatment = "eggs_constant")

params_eggs_constant_all <- bind_rows(params_eggs_constant, topt_eggs_constant)

# write_csv(params_eggs_constant_all, "data-processed/params_eggs_constant_all.csv")

##### Create & save df with summarized data for plotting
data_eggs_constant_sum <- data_eggs_constant %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = sd(trait)/sqrt(n())) %>% 
	mutate(trait = "bite rate") %>% 
	mutate(treatment = "eggs_constant")

# write_csv(data_eggs_constant_sum, "data-processed/data_eggs_constant_sum.csv")



############## 3B. Lifetime Eggs: DTR 9 ---------------------------------------------------------------

##### Extract matrix of predicted values from fitted model (451 col = temperature gradient, 7500 rows = MCMC steps)
predictions_eggs_dtr9 <- as.data.frame(model_eggs_dtr9$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
colnames(predictions_eggs_dtr9) <- Temp.xs

##### Create & save summary df with 2.5% and 97.5% CIs, mean, and median trait values for each temperature
predictions_eggs_dtr9_summary <- predictions_eggs_dtr9 %>%
	mutate(iteration = rownames(.)) %>% 
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(temperature) %>%  
	summarise(lowerCI = quantile(trait_value, probs = 0.025),
			  upperCI = quantile(trait_value, probs = 0.975),
			  mean = mean(trait_value),
			  median = median(trait_value)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "eggs_dtr9")

# write_csv(predictions_eggs_dtr9_summary, "data-processed/predictions_eggs_dtr9_summary.csv")

##### Create summary df with 2.5% and 97.5% CIs, mean, and median Topt values for the TPC
topt_eggs_dtr9 <- predictions_eggs_dtr9 %>% 
	mutate(iteration = rownames(.)) %>%
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(iteration) %>% 
	slice_max(order_by = trait_value, n = 1) %>%
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(lowerCI = quantile(temperature, probs = 0.025),
			  upperCI = quantile(temperature, probs = 0.975),
			  mean = mean(temperature),
			  median = median(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "eggs_dtr9")

##### Create & save df with fitted parameters and Topt values from above
params_eggs_dtr9 <- as.data.frame(model_eggs_dtr9$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>%
	rename(lowerCI = `2.5%`, upperCI = `97.5%`, median = `50%`) %>%
	mutate(treatment = "eggs_dtr9")

params_eggs_dtr9_all <- bind_rows(params_eggs_dtr9, topt_eggs_dtr9)

# write_csv(params_eggs_dtr9_all, "data-processed/params_eggs_dtr9_all.csv")

##### Create & save df with summarized data for plotting
data_eggs_dtr9_sum <- data_eggs_dtr9 %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = sd(trait)/sqrt(n())) %>% 
	mutate(trait = "bite rate") %>% 
	mutate(treatment = "eggs_dtr9")

# write_csv(data_eggs_dtr9_sum, "data-processed/data_eggs_dtr9_sum.csv")



############## 3C. Lifetime Eggs: DTR 12 ---------------------------------------------------------------

##### Extract matrix of predicted values from fitted model (451 col = temperature gradient, 7500 rows = MCMC steps)
predictions_eggs_dtr12 <- as.data.frame(model_eggs_dtr12$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
colnames(predictions_eggs_dtr12) <- Temp.xs

##### Create & save summary df with 2.5% and 97.5% CIs, mean, and median trait values for each temperature
predictions_eggs_dtr12_summary <- predictions_eggs_dtr12 %>%
	mutate(iteration = rownames(.)) %>% 
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(temperature) %>%  
	summarise(lowerCI = quantile(trait_value, probs = 0.025),
			  upperCI = quantile(trait_value, probs = 0.975),
			  mean = mean(trait_value),
			  median = median(trait_value)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "eggs_dtr12")

# write_csv(predictions_eggs_dtr12_summary, "data-processed/predictions_eggs_dtr12_summary.csv")

##### Create summary df with 2.5% and 97.5% CIs, mean, and median Topt values for the TPC
topt_eggs_dtr12 <- predictions_eggs_dtr12 %>% 
	mutate(iteration = rownames(.)) %>%
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(iteration) %>% 
	slice_max(order_by = trait_value, n = 1) %>%
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(lowerCI = quantile(temperature, probs = 0.025),
			  upperCI = quantile(temperature, probs = 0.975),
			  mean = mean(temperature),
			  median = median(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "eggs_dtr12")

##### Create & save df with fitted parameters and Topt values from above
params_eggs_dtr12 <- as.data.frame(model_eggs_dtr12$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>%
	rename(lowerCI = `2.5%`, upperCI = `97.5%`, median = `50%`) %>%
	mutate(treatment = "eggs_dtr12")

params_eggs_dtr12_all <- bind_rows(params_eggs_dtr12, topt_eggs_dtr12)

# write_csv(params_eggs_dtr12_all, "data-processed/params_eggs_dtr12_all.csv")

##### Create & save df with summarized data for plotting
data_eggs_dtr12_sum <- data_eggs_dtr12 %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = sd(trait)/sqrt(n())) %>% 
	mutate(trait = "bite rate") %>% 
	mutate(treatment = "eggs_dtr12")

# write_csv(data_eggs_dtr12_sum, "data-processed/data_eggs_dtr12_sum.csv")



############## 4. vector competence (bc) -------------------------------------------------------

##### Extract matrix of predicted values from fitted model (451 col = temperature gradient, 7500 rows = MCMC steps)
predictions_bc <- as.data.frame(model_bc$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
colnames(predictions_bc) <- Temp.xs

##### Create & save summary df with 2.5% and 97.5% CIs, mean, and median trait values for each temperature
predictions_bc_summary <- predictions_bc %>%
	mutate(iteration = rownames(.)) %>% 
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(temperature) %>%  
	summarise(lowerCI = quantile(trait_value, probs = 0.025),
			  upperCI = quantile(trait_value, probs = 0.975),
			  mean = mean(trait_value),
			  median = median(trait_value)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "bc")

# write_csv(predictions_bc_summary, "data-processed/predictions_bc_summary.csv")

##### Create summary df with 2.5% and 97.5% CIs, mean, and median Topt values for the TPC
topt_bc <- predictions_bc %>% 
	mutate(iteration = rownames(.)) %>%
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(iteration) %>% 
	slice_max(order_by = trait_value, n = 1) %>%
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(lowerCI = quantile(temperature, probs = 0.025),
			  upperCI = quantile(temperature, probs = 0.975),
			  mean = mean(temperature),
			  median = median(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "bc")

##### Create & save df with fitted parameters and Topt values from above
params_bc <- as.data.frame(model_bc$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>%
	rename(lowerCI = `2.5%`, upperCI = `97.5%`, median = `50%`) %>%
	mutate(treatment = "bc")

params_bc_all <- bind_rows(params_bc, topt_bc)

# write_csv(params_bc_all, "data-processed/params_bc_all.csv")

##### Create & save df with summarized data for plotting
data_bc_sum <- data_bc %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = sd(trait)/sqrt(n())) %>%
	mutate(trait = "vector competence") %>% 
	mutate(treatment = "bc")

# write_csv(data_bc_sum, "data-processed/data_bc_sum.csv")



############## 5. EIP50 ---------------------------------------------------------------

##### Extract matrix of predicted values from fitted model (451 col = temperature gradient, 7500 rows = MCMC steps)
predictions_eip50 <- as.data.frame(model_eip50$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
colnames(predictions_eip50) <- Temp.xs

##### Create & save summary df with 2.5% and 97.5% CIs, mean, and median trait values for each temperature
predictions_eip50_summary <- predictions_eip50 %>%
	mutate(iteration = rownames(.)) %>% 
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(temperature) %>%  
	summarise(lowerCI = quantile(trait_value, probs = 0.025),
			  upperCI = quantile(trait_value, probs = 0.975),
			  mean = mean(trait_value),
			  median = median(trait_value)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "eip50")

# write_csv(predictions_eip50_summary, "data-processed/predictions_eip50_summary.csv")

##### Create summary df with 2.5% and 97.5% CIs, mean, and median Topt values for the TPC
topt_eip50 <- predictions_eip50 %>%
	mutate(iteration = rownames(.)) %>%
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(iteration) %>% 
	slice_max(order_by = trait_value, n = 1) %>%
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(lowerCI = quantile(temperature, probs = 0.025),
			  upperCI = quantile(temperature, probs = 0.975),
			  mean = mean(temperature),
			  median = median(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "eip50")

##### Create & save df with fitted parameters and Topt values from above
params_eip50 <- as.data.frame(model_eip50$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>%
	rename(lowerCI = `2.5%`, upperCI = `97.5%`, median = `50%`) %>%
	mutate(treatment = "eip50")

params_eip50_all <- bind_rows(params_eip50, topt_eip50)

# write_csv(params_eip50_all, "data-processed/params_eip50_all.csv")

##### Create & save df with summarized data for plotting
data_eip50_sum <- data_eip50 %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = sd(trait)/sqrt(n())) %>%
	mutate(trait = "eip50") %>% 
	mutate(treatment = "eip50")

# write_csv(data_eip50_sum, "data-processed/data_eip50_sum.csv")



############## 6. pEA -----------------------------------------------------------------

##### Extract matrix of predicted values from fitted model (451 col = temperature gradient, 7500 rows = MCMC steps)
predictions_pea <- as.data.frame(model_pea$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
colnames(predictions_pea) <- Temp.xs

##### Create & save summary df with 2.5% and 97.5% CIs, mean, and median trait values for each temperature
predictions_pea_summary <- predictions_pea %>%
	mutate(iteration = rownames(.)) %>% 
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(temperature) %>%  
	summarise(lowerCI = quantile(trait_value, probs = 0.025),
			  upperCI = quantile(trait_value, probs = 0.975),
			  mean = mean(trait_value),
			  median = median(trait_value)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "pea")

# write_csv(predictions_pea_summary, "data-processed/predictions_pea_summary.csv")

##### Create summary df with 2.5% and 97.5% CIs, mean, and median Topt values for the TPC
topt_pea <- predictions_pea %>% 
	mutate(iteration = rownames(.)) %>%
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(iteration) %>% 
	slice_max(order_by = trait_value, n = 1) %>%
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(lowerCI = quantile(temperature, probs = 0.025),
			  upperCI = quantile(temperature, probs = 0.975),
			  mean = mean(temperature),
			  median = median(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "pea")

##### Create & save df with fitted parameters and Topt values from above
params_pea <- as.data.frame(model_pea$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>%
	rename(lowerCI = `2.5%`, upperCI = `97.5%`, median = `50%`) %>%
	mutate(treatment = "pea")

params_pea_all <- bind_rows(params_pea, topt_pea)

# write_csv(params_pea_all, "data-processed/params_pea_all.csv")

##### Create & save df with summarized data for plotting
data_pea_sum <- data_pea %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = sd(trait)/sqrt(n())) %>%
	mutate(trait = "pea") %>% 
	mutate(treatment = "pea")

# write_csv(data_pea_sum, "data-processed/data_pea_sum.csv")



############## 7. MDR -----------------------------------------------------------------

##### Extract matrix of predicted values from fitted model (451 col = temperature gradient, 7500 rows = MCMC steps)
predictions_mdr <- as.data.frame(model_mdr$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
colnames(predictions_mdr) <- Temp.xs

##### Create & save summary df with 2.5% and 97.5% CIs, mean, and median trait values for each temperature
predictions_mdr_summary <- predictions_mdr %>%
	mutate(iteration = rownames(.)) %>% 
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(temperature) %>%  
	summarise(lowerCI = quantile(trait_value, probs = 0.025),
			  upperCI = quantile(trait_value, probs = 0.975),
			  mean = mean(trait_value),
			  median = median(trait_value)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "mdr")

# write_csv(predictions_mdr_summary, "data-processed/predictions_mdr_summary.csv")

##### Create summary df with 2.5% and 97.5% CIs, mean, and median Topt values for the TPC
topt_mdr <- predictions_mdr %>% 
	mutate(iteration = rownames(.)) %>%
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(iteration) %>% 
	slice_max(order_by = trait_value, n = 1) %>%
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(lowerCI = quantile(temperature, probs = 0.025),
			  upperCI = quantile(temperature, probs = 0.975),
			  mean = mean(temperature),
			  median = median(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "mdr")

##### Create & save df with fitted parameters and Topt values from above
params_mdr <- as.data.frame(model_mdr$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>%
	rename(lowerCI = `2.5%`, upperCI = `97.5%`, median = `50%`) %>%
	mutate(treatment = "mdr")

params_mdr_all <- bind_rows(params_mdr, topt_mdr)

# write_csv(params_mdr_all, "data-processed/params_mdr_all.csv")

##### Create & save df with summarized data for plotting
data_mdr_sum <- data_mdr %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = sd(trait)/sqrt(n())) %>%
	mutate(trait = "mdr") %>% 
	mutate(treatment = "mdr")

# write_csv(data_mdr_sum, "data-processed/data_mdr_sum.csv")



############## 8A. Gamma: constant ---------------------------------------------------------------

##### Extract matrix of predicted values from fitted model (451 col = temperature gradient, 7500 rows = MCMC steps)
predictions_gamma_constant <- as.data.frame(model_gamma_constant$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
colnames(predictions_gamma_constant) <- Temp.xs

##### Create & save summary df with 2.5% and 97.5% CIs, mean, and median trait values for each temperature
predictions_gamma_constant_summary <- predictions_gamma_constant %>%
	mutate(iteration = rownames(.)) %>% 
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(temperature) %>%  
	summarise(lowerCI = quantile(trait_value, probs = 0.025),
			  upperCI = quantile(trait_value, probs = 0.975),
			  mean = mean(trait_value),
			  median = median(trait_value)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "gamma_constant")

# write_csv(predictions_gamma_constant_summary, "data-processed/predictions_gamma_constant_summary.csv")

##### Create summary df with 2.5% and 97.5% CIs, mean, and median Topt values for the TPC
topt_gamma_constant <- predictions_gamma_constant %>% 
	mutate(iteration = rownames(.)) %>%
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(iteration) %>% 
	slice_max(order_by = trait_value, n = 1) %>%
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(lowerCI = quantile(temperature, probs = 0.025),
			  upperCI = quantile(temperature, probs = 0.975),
			  mean = mean(temperature),
			  median = median(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "gamma_constant")

##### Create & save df with fitted parameters and Topt values from above
params_gamma_constant <- as.data.frame(model_gamma_constant$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>%
	rename(lowerCI = `2.5%`, upperCI = `97.5%`, median = `50%`) %>%
	mutate(treatment = "gamma_constant")

params_gamma_constant_all <- bind_rows(params_gamma_constant, topt_gamma_constant)

# write_csv(params_gamma_constant_all, "data-processed/params_gamma_constant_all.csv")

##### Create & save df with summarized data for plotting
data_gamma_constant_sum <- data_gamma_constant %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = sd(trait)/sqrt(n())) %>% 
	mutate(trait = "bite rate") %>% 
	mutate(treatment = "gamma_constant")

# write_csv(data_gamma_constant_sum, "data-processed/data_gamma_constant_sum.csv")


############## 8B. Gamma: DTR 9 ---------------------------------------------------------------

##### Extract matrix of predicted values from fitted model (451 col = temperature gradient, 7500 rows = MCMC steps)
predictions_gamma_dtr9 <- as.data.frame(model_gamma_dtr9$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
colnames(predictions_gamma_dtr9) <- Temp.xs

##### Create & save summary df with 2.5% and 97.5% CIs, mean, and median trait values for each temperature
predictions_gamma_dtr9_summary <- predictions_gamma_dtr9 %>%
	mutate(iteration = rownames(.)) %>% 
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(temperature) %>%  
	summarise(lowerCI = quantile(trait_value, probs = 0.025),
			  upperCI = quantile(trait_value, probs = 0.975),
			  mean = mean(trait_value),
			  median = median(trait_value)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "gamma_dtr9")

# write_csv(predictions_gamma_dtr9_summary, "data-processed/predictions_gamma_dtr9_summary.csv")

##### Create summary df with 2.5% and 97.5% CIs, mean, and median Topt values for the TPC
topt_gamma_dtr9 <- predictions_gamma_dtr9 %>% 
	mutate(iteration = rownames(.)) %>%
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(iteration) %>% 
	slice_max(order_by = trait_value, n = 1) %>%
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(lowerCI = quantile(temperature, probs = 0.025),
			  upperCI = quantile(temperature, probs = 0.975),
			  mean = mean(temperature),
			  median = median(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "gamma_dtr9")

##### Create & save df with fitted parameters and Topt values from above
params_gamma_dtr9 <- as.data.frame(model_gamma_dtr9$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>%
	rename(lowerCI = `2.5%`, upperCI = `97.5%`, median = `50%`) %>%
	mutate(treatment = "gamma_dtr9")

params_gamma_dtr9_all <- bind_rows(params_gamma_dtr9, topt_gamma_dtr9)

# write_csv(params_gamma_dtr9_all, "data-processed/params_gamma_dtr9_all.csv")

##### Create & save df with summarized data for plotting
data_gamma_dtr9_sum <- data_gamma_dtr9 %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = sd(trait)/sqrt(n())) %>% 
	mutate(trait = "bite rate") %>% 
	mutate(treatment = "gamma_dtr9")

# write_csv(data_gamma_dtr9_sum, "data-processed/data_gamma_dtr9_sum.csv")


############## 8C. Gamma: DTR 12 ---------------------------------------------------------------

##### Extract matrix of predicted values from fitted model (451 col = temperature gradient, 7500 rows = MCMC steps)
predictions_gamma_dtr12 <- as.data.frame(model_gamma_dtr12$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
colnames(predictions_gamma_dtr12) <- Temp.xs

##### Create & save summary df with 2.5% and 97.5% CIs, mean, and median trait values for each temperature
predictions_gamma_dtr12_summary <- predictions_gamma_dtr12 %>%
	mutate(iteration = rownames(.)) %>% 
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(temperature) %>%  
	summarise(lowerCI = quantile(trait_value, probs = 0.025),
			  upperCI = quantile(trait_value, probs = 0.975),
			  mean = mean(trait_value),
			  median = median(trait_value)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "gamma_dtr12")

# write_csv(predictions_gamma_dtr12_summary, "data-processed/predictions_gamma_dtr12_summary.csv")

##### Create summary df with 2.5% and 97.5% CIs, mean, and median Topt values for the TPC
topt_gamma_dtr12 <- predictions_gamma_dtr12 %>% 
	mutate(iteration = rownames(.)) %>%
	pivot_longer(!iteration, names_to = "temperature", values_to = "trait_value") %>%
	group_by(iteration) %>% 
	slice_max(order_by = trait_value, n = 1) %>%
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(lowerCI = quantile(temperature, probs = 0.025),
			  upperCI = quantile(temperature, probs = 0.975),
			  mean = mean(temperature),
			  median = median(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "gamma_dtr12")

##### Create & save df with fitted parameters and Topt values from above
params_gamma_dtr12 <- as.data.frame(model_gamma_dtr12$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>%
	rename(lowerCI = `2.5%`, upperCI = `97.5%`, median = `50%`) %>%
	mutate(treatment = "gamma_dtr12")

params_gamma_dtr12_all <- bind_rows(params_gamma_dtr12, topt_gamma_dtr12)

# write_csv(params_gamma_dtr12_all, "data-processed/params_gamma_dtr12_all.csv")

##### Create & save df with summarized data for plotting
data_gamma_dtr12_sum <- data_gamma_dtr12 %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = sd(trait)/sqrt(n())) %>% 
	mutate(trait = "bite rate") %>% 
	mutate(treatment = "gamma_dtr12")

# write_csv(data_gamma_dtr12_sum, "data-processed/data_gamma_dtr12_sum.csv")







# #####################
# #####Load Johnson fits, use Johnson fits for EIP, bc, MDR, and pEA as informative priors
# ####################
# 
# # Load prior fits
# load("Johnson et al posteriors/PriorFitsList.Rdata")
# 
# 
# 
# ##########
# ###### 4. Fitting prior distributions from the Johnson et al TPCs
# ##########
# 
# # Function to fit gamma distributions for q, T0, Tm, and sigma from fitted TPCs
# prior.fitting = function(model.out){
#   # Fit gamma distributions for each parameter posterior dists
#   prior.fits = suppressWarnings(apply(model.out, 2, function(df) fitdistr(df, "gamma")$estimate))
#   prior.fits
# } 
# 
# # Load Prior Fits
# posts = c(
#   "a_posterior1",
#   "PDR_post_prior1",
#   "MDR_posterior1",
#   "EFD_posterior1",
#   "e2a_posterior1",
#   "bc_posterior1_quad",
#   "mu_posterior1")
# fits = paste("Johnson et al posteriors/", posts, ".Rsave", sep = "")
# for (i in 1:length(fits)){
#   load(fits[i])
#   assign(posts[i], samps)
# }
# 
# prior.fits.list = lapply(list(a_posterior1, PDR_post_prior1, MDR_posterior1, EFD_posterior1[,-3], e2a_posterior1[,-3], bc_posterior1_quad[,-3], mu_posterior1), prior.fitting)
# names(prior.fits.list) = posts
# 
# 
# 
# ############ EIP50 from Shapiro et al. 2017 Plos Biology
# # Get data
# data.specific <- with(data.bc.EIP, data.frame('T' = temp, 'trait' = 1/EIP50)) # subset specific trait data from complete list
# data <- data.specific # assign trait data to variable 'data'
# trait.prior.fits = prior.fits.list$PDR_post_prior1[,c('c', 'T0', 'Tm')]
# hypers <- trait.prior.fits * 0.1 # assign priors to variable 'priors' and weight importance if desired
# 
# # Organize Data for JAGS
# trait <- data$trait
# N.obs <- length(trait)
# temp <- data$T
# 
# # Bundle Data
# jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)
# 
# # Run JAGS - **select correct model file**
# model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
#                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())
# 
# # Examine model output & run diagnostics
# model.out$BUGSoutput$summary[1:10,]
# mcmcplot(model.out)
# 
# # Save model output 
# save(model.out, file = "saved posteriors/EIP50_briere_inf.Rdata")
# 
# # Plot trait data, model mean and CIs
# plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), ylim = c(0, 0.2), data = data.specific, 
#      ylab = "EIP-50", xlab = expression(paste("Temperature (",degree,"C)")))
# lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
# lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
# lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
# 
# 
# ############ Vector competence from Shapiro et al. 2017 Plos Biology
# # Get data
# data.specific <- with(data.bc.EIP, data.frame('T' = temp, 'trait' = bc)) # subset specific trait data from complete list
# data <- data.specific # assign trait data to variable 'data'
# trait.prior.fits = prior.fits.list$bc_posterior1_quad[,c('n.qd', 'T0', 'Tm')]
# hypers <- trait.prior.fits * 0.5 # assign priors to variable 'priors' and weight importance if desired
# 
# # Organize Data for JAGS
# trait <- data$trait
# N.obs <- length(trait)
# temp <- data$T
# 
# # Bundle Data
# jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)
# 
# # Run JAGS - **select correct model file**
# model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
#                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())
# 
# # Examine model output & run diagnostics
# model.out$BUGSoutput$summary[1:10,]
# mcmcplot(model.out)
# 
# # Save model output 
# save(model.out, file = "saved posteriors/bc_quadratic_inf.Rdata")
# 
# # Plot trait data, model mean and CIs
# plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), ylim = c(0,0.7), data = data.specific, 
#      ylab = "Vector competence", xlab = expression(paste("Temperature (",degree,"C)")))
# lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
# lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
# lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
# 
# ############ Mosqutio developement rate from Paaijmans et al. 2013 Global Climate - Briere
# # Get data
# data.specific <- with(data.pea.MDR, data.frame('T' = temp, 'trait' = MDR)) # subset specific trait data from complete list
# data <- data.specific # assign trait data to variable 'data'
# trait.prior.fits = prior.fits.list$MDR_posterior1[,c('c', 'T0', 'Tm')]
# hypers <- trait.prior.fits * 0.5 # assign priors to variable 'priors' and weight importance if desired
# 
# # Organize Data for JAGS
# trait <- data$trait
# N.obs <- length(trait)
# temp <- data$T
# 
# # Bundle Data
# jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)
# 
# # Run JAGS - **select correct model file**
# model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
#                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())
# 
# # Examine model output & run diagnostics
# model.out$BUGSoutput$summary[1:10,]
# mcmcplot(model.out)
# 
# # Save model output
# save(model.out, file = "saved posteriors/MDR_briere_inf.Rdata")
# 
# # Plot trait data, model mean and CIs
# plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, xaxt = "n",
#      ylab = "Mosquito development rate", xlab = expression(paste("Temperature (",degree,"C)")))
# lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
# lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
# lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
# 
# ########## Prob. egg to adult survival from Paaijmans et al. 2013 Global Climate - Quadratic
# # Get data
# data.specific <- with(data.pea.MDR, data.frame('T' = temp, 'trait' = Pea)) # subset specific trait data from complete list
# data <- data.specific # assign trait data to variable 'data'
# trait.prior.fits = prior.fits.list$e2a_posterior1[,c('n.qd', 'T0', 'Tm')]
# hypers <- trait.prior.fits * 0.5 # assign priors to variable 'priors' and weight importance if desired
# 
# # Organize Data for JAGS
# trait <- data$trait
# N.obs <- length(trait)
# temp <- data$T
# 
# # Bundle Data
# jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)
# 
# # Run JAGS - **select correct model file**
# model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
#                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())
# 
# # Examine model output & run diagnostics
# model.out$BUGSoutput$summary[1:10,]
# mcmcplot(model.out)
# 
# # Save model output 
# save(model.out, file = "saved posteriors/pea_quadratic_inf.Rdata")
# 
# # Plot trait data, model mean and CIs
# plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), ylim = c(0,0.7), data = data.specific, 
#      ylab = "Prob. e2a", xlab = expression(paste("Temperature (",degree,"C)")))
# lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
# lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
# lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
