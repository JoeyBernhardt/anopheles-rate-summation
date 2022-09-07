## Joey Bernhardt
## Modified from code provided by Kerri Miazgowicz, which was modified from code by Erin Mordecai on August 6, 2018


# Load libraties for fitting traits
library('R2jags') # Fits Bayesian models
library('mcmcplots') # Diagnostic plots for fits
library('MASS') # Fits distributions for informative priors
library(plyr) # Slices and dices data
library(plotrix) # For standard error function
library(tidyverse)


data.constant <- read.csv("data-raw/constant.individual.trait.csv")
names(data.constant)[names(data.constant) == 'Treatment'] <- 'temp'
data.fluc9 <- read.csv("data-raw/fluc9.individual.trait.csv")
names(data.fluc9)[names(data.fluc9) == 'Treatment'] <- 'temp'
data.fluc12 <- read.csv("data-raw/fluc12.individual.trait.csv")
names(data.fluc12)[names(data.fluc12) == 'Treatment'] <- 'temp'

data.bc.EIP <- read.csv("data-raw/forErin_ShapiroData.csv")
data.bc.EIP$inverse.EIP50 = 1/data.bc.EIP$EIP50
data.pea.MDR <- read.csv("data-raw/Krijn_Raw_Data.csv")

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
	cf.Tm = 40,
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
# save(Temp.xs, file = "saved-posteriors/temps.Rdata")

##########
###### 4. Fitting traits - uniform priors
##########

############## First plot each trait to check functional forms
# boxplot(bite.rate ~ temp, data = data.constant) # briere
# boxplot(bite.rate ~ temp, data = data.fluc9) # briere
# boxplot(bite.rate ~ temp, data = data.fluc12) # briere /linear
# boxplot(lifespan ~ temp, data = data.constant) # quad
# boxplot(lifespan ~ temp, data = data.fluc9) # quad
# boxplot(lifespan ~ temp, data = data.fluc12) #quad
# boxplot(lifetime.eggs ~ temp, data = data.constant) # quad
# boxplot(lifetime.eggs ~ temp, data = data.fluc9) #quad
# boxplot(lifetime.eggs ~ temp, data = data.fluc12) #quad
# 
# plot(bc ~ temp, data = data.bc.EIP, pch = 16) # quad?
# plot(1/EIP50 ~ temp, data = data.bc.EIP, pch = 16) # briere / linear
# plot(Pea~ temp, data = data.pea.MDR) #quadratic 
# plot(MDR ~temp, data = data.pea.MDR) #briere

############## Trait name and parameter symbol - model type (e.g. quadratic, Briere)



# Model fitting: Bite rate constant temperature ------------------------------------------


############ Bite rate at constant temperature
# Get data
# data.specific <- with(data.constant, data.frame('T' = temp, 'trait' = bite.rate)) # subset specific trait data from complete list
# data_constant_bite_rate <- data.specific # assign trait data to variable 'data'

data_constant_bite_rate <- with(data.constant, data.frame('T' = temp, 'trait' = bite.rate)) # subset specific trait data from complete list
data_constant_bite_rate_2 <- data_constant_bite_rate %>% 
	mutate(treatment = "bite_rate_constant")
# Organize Data for JAGS
trait <- data_constant_bite_rate_2$trait
N.obs <- length(trait)
temp <- data_constant_bite_rate_2$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
model_out_bite_rate_constant_briere <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
											n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model_out_bite_rate_constant_briere$BUGSoutput$summary[1:10,]
mcmcplot(model_out_bite_rate_constant_briere)

# Save model output 
# save(model_out_bite_rate_constant_briere, file = "saved-posteriors/constant_bite.rate_briere_uniform.Rdata")
# constant_briere_bite_rate <- load("saved-posteriors/constant_bite.rate_briere_uniform.Rdata")

#### now plot in ggplot
b_params_bite_rate <- as.data.frame(model_out_bite_rate_constant_briere$BUGSoutput$summary[1:5,]) %>% 
	rownames_to_column(var = "term")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_constant_bite_rate_2, 
	 ylab = "Bite rate", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_out_bite_rate_constant_briere$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_out_bite_rate_constant_briere$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_out_bite_rate_constant_briere$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(trait ~ T, xlim = c(0, 45), ylim = c(0,1), data = data, 
	 ylab = "Bite rate", xlab = expression(paste("Temperature (",degree,"C)")))



# Model fitting: Bite rate at DTR 9 ---------------------------------------


############ Bite rate at DTR 9C
# Get data
data.specific_bite_rate_dtr9 <- with(data.fluc9, data.frame('T' = temp, 'trait' = bite.rate)) # subset specific trait data from complete list
# data <- data.specific # assign trait data to variable 'data'

# Organize Data for JAGS
trait <- data.specific_bite_rate_dtr9$trait
N.obs <- length(trait)
temp <- data.specific_bite_rate_dtr9$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
model.out_bite_rate_dtr9_briere <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
										n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# # Examine model output & run diagnostics
# model.out$BUGSoutput$summary[1:10,]
# mcmcplot(model.out)
# 
# # Save model output 
# save(model.out_bite_rate_dtr9_briere, file = "saved-posteriors/dtr9_bite.rate_briere_uniform.Rdata")
# 
# # Plot trait data, model mean and CIs
# plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
#      ylab = "Bite rate", xlab = expression(paste("Temperature (",degree,"C)")))
# lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
# lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
# lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



# Model fitting: Bite rate at DTR 12 --------------------------------------


############ Bite rate at DTR 12C
# Get data
data_bite_rate_dtr12 <- with(data.fluc12, data.frame('T' = temp, 'trait' = bite.rate)) # subset specific trait data from complete list
# data_bite_rate_dtr12 <- data.specific # assign trait data to variable 'data'

# Organize Data for JAGS
trait <- data_bite_rate_dtr12$trait
N.obs <- length(trait)
temp <- data_bite_rate_dtr12$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
model.out_bite_rate_dtr12_briere <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
										 n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# # Examine model output & run diagnostics
# model.out$BUGSoutput$summary[1:10,]
# mcmcplot(model.out)
# 
# # Save model output 
# save(model.out_bite_rate_dtr12_briere, file = "saved-posteriors/dtr12_bite.rate_briere_uniform.Rdata")
# 
# # Plot trait data, model mean and CIs
# plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data.specific, 
#      ylab = "Bite rate", xlab = expression(paste("Temperature (",degree,"C)")))
# lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
# lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
# lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
# 


# Model fitting: Lifespan constant ---------------------------------------------------------



############ Lifespan at constant temperature
# Get data
data_constant_lifespan <- with(data.constant, data.frame('T' = temp, 'trait' = lifespan)) # subset specific trait data from complete list
# data <- data.specific # assign trait data to variable 'data'

# Organize Data for JAGS
trait <- data_constant_lifespan$trait
N.obs <- length(trait)
temp <- data_constant_lifespan$T

# trait <- data$trait
# N.obs <- length(trait)
# temp <- data$T

jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
model_out_constant_lifespan_briere <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
										   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# # Bundle Data
# # jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
# 
# # Run JAGS - **select correct model file**
# # data_constant_lifespan
# # model_out_constant_lifespan_briere <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
# #                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())
# 
# # Examine model output & run diagnostics
# # model.out$BUGSoutput$summary[1:10,]
# # mcmcplot(model.out)
# 
# # Save model output 
# save(model_out_constant_lifespan_briere, file = "saved-posteriors/constant_lifespan_briere_uniform.Rdata")
# 
# # Plot trait data, model mean and CIs
# plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_constant_lifespan,
#      ylab = "Lifespan", xlab = expression(paste("Temperature (",degree,"C)")))
# lines(model_out_constant_lifespan_briere$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
# lines(model_out_constant_lifespan_briere$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
# lines(model_out_constant_lifespan_briere$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
# 

# Model fitting: Lifespan at DTR 9 ----------------------------------------


############ Lifespan at DTR 9C
# Get data
data_dtr9_lifespan <- with(data.fluc9, data.frame('T' = temp, 'trait' = lifespan)) # subset specific trait data from complete list
# data <- data.specific # assign trait data to variable 'data'

# Organize Data for JAGS
trait <- data_dtr9_lifespan$trait
N.obs <- length(trait)
temp <- data_dtr9_lifespan$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
data_dtr9_lifespan
model_out_lifespan_dtr9_briere <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
									   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
# model_out_lifespan_dtr9_briere$BUGSoutput$summary[1:10,]
# mcmcplot(model_out_lifespan_dtr9_briere)
# 
# # Save model output 
# save(model_out_lifespan_dtr9_briere, file = "saved-posteriors/dtr9_lifespan_briere_uniform.Rdata")

# # Plot trait data, model mean and CIs
# plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_dtr9_lifespan, 
#      ylab = "Lifespan", xlab = expression(paste("Temperature (",degree,"C)")))
# lines(model_out_lifespan_dtr9_briere$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
# lines(model_out_lifespan_dtr9_briere$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
# lines(model_out_lifespan_dtr9_briere$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



# Model fitting: Lifespan DTR12 ----------------------------------------------------------

############ Lifespan at DTR 12C
# Get data
data_dtr12_lifespan <- with(data.fluc12, data.frame('T' = temp, 'trait' = lifespan)) # subset specific trait data from complete list
# data <- data.specific # assign trait data to variable 'data'

# Organize Data for JAGS
trait <- data_dtr12_lifespan$trait
N.obs <- length(trait)
temp <- data_dtr12_lifespan$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
# data_dtr12_lifespan
model_out_lifespan_dtr12 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
								 n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# # Examine model output & run diagnostics
# model_out_lifespan_dtr12$BUGSoutput$summary[1:10,]
# mcmcplot(model_out_lifespan_dtr12)
# 
# # Save model output 
# save(model_out_lifespan_dtr12, file = "saved-posteriors/dtr12_lifespan_briere_uniform.Rdata")
# 
# # Plot trait data, model mean and CIs
# plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_dtr12_lifespan, 
#      ylab = "Lifespan", xlab = expression(paste("Temperature (",degree,"C)")))
# lines(model_out_lifespan_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
# lines(model_out_lifespan_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
# lines(model_out_lifespan_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



# Model fitting: Lifetime eggs constant -----------------------------------


############ Lifetime eggs at constant temperature
# Get data
# data_constant_eggs <- with(data.constant, data.frame('T' = temp, 'trait' = lifetime.eggs)) # subset specific trait data from complete list
# # data <- data.specific # assign trait data to variable 'data'
# 
# # Organize Data for JAGS
# trait <- data_constant_eggs$trait
# N.obs <- length(trait)
# temp <- data_constant_eggs$T
# 
# # Bundle Data
# jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
# 
# # Run JAGS - **select correct model file**
# data_constant_eggs
# model_out_eggs_constant <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
#                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# # Examine model output & run diagnostics
# model_out_eggs_constant$BUGSoutput$summary[1:10,]
# mcmcplot(model_out_eggs_constant)
# 
# str(model_out_eggs_constant)
# 
# # Save model output 
# save(model_out_eggs_constant, file = "saved-posteriors/constant_eggs_briere_uniform.Rdata")
# 
# # Plot trait data, model mean and CIs
# plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_constant_eggs, 
#      ylab = "Lifetime eggs", xlab = expression(paste("Temperature (",degree,"C)")))
# lines(model_out_eggs_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
# lines(model_out_eggs_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
# lines(model_out_eggs_constant$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
# 

# Model fitting: Lifetime eggs at DTR9 ------------------------------------


############ Lifetime eggs at DTR 9C
# Get data
data_eggs_dtr9 <- with(data.fluc9, data.frame('T' = temp, 'trait' = lifetime.eggs)) # subset specific trait data from complete list
# data <- data.specific # assign trait data to variable 'data'

# Organize Data for JAGS
trait <- data_eggs_dtr9$trait
N.obs <- length(trait)
temp <- data_eggs_dtr9$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
data_eggs_dtr9
model_out_eggs_dtr9 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
							n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# # Examine model output & run diagnostics
# model_out_eggs_dtr9$BUGSoutput$summary[1:10,]
# mcmcplot(model_out_eggs_dtr9)
# 
# # Save model output 
# save(model_out_eggs_dtr9, file = "saved-posteriors/dtr9_eggs_briere_uniform.Rdata")
# 
# # Plot trait data, model mean and CIs
# plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_eggs_dtr9, 
#      ylab = "Lifetime eggs", xlab = expression(paste("Temperature (",degree,"C)")))
# lines(model_out_eggs_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
# lines(model_out_eggs_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
# lines(model_out_eggs_dtr9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


# Model fitting: Lifetime eggs DTR12 --------------------------------------


############ Lifetime eggs at dtr12
# Get data
data_eggs_dtr12 <- with(data.fluc12, data.frame('T' = temp, 'trait' = lifetime.eggs)) # subset specific trait data from complete list
# data <- data.specific # assign trait data to variable 'data'

# Organize Data for JAGS
trait <- data_eggs_dtr12$trait
N.obs <- length(trait)
temp <- data_eggs_dtr12$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
data_eggs_dtr12
model_out_eggs_dtr12 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
							 n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# # Examine model output & run diagnostics
# model_out_eggs_dtr12$BUGSoutput$summary[1:10,]
# mcmcplot(model_out_eggs_dtr12)
# 
# # Save model output 
# save(model_out_eggs_dtr12, file = "saved-posteriors/dtr12_eggs_briere_uniform.Rdata")
# 
# # Plot trait data, model mean and CIs
# plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_eggs_dtr12, 
#      ylab = "Lifetime eggs", xlab = expression(paste("Temperature (",degree,"C)")))
# lines(model_out_eggs_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
# lines(model_out_eggs_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
# lines(model_out_eggs_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



#### The parameters of S(T) include: bite rate (a), vector competence (bc; the proportion of infectious mosquitoes),
# adult mosquito lifespan (lf), lifetime egg production (B), probability of egg to adult survival (pEA),
# mosquito development rate (MDR)


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
data_eip50
model_out_eip50 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
						n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# # Examine model output & run diagnostics
# model_out_eip50$BUGSoutput$summary[1:10,]
# mcmcplot(model_out_eip50)
# 
# # Save model output 
# save(model_out_eip50, file = "saved-posteriors/EIP50_briere_uniform.Rdata")
# 
# # Plot trait data, model mean and CIs
# plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_eip50, 
#      ylab = "EIP-50", xlab = expression(paste("Temperature (",degree,"C)")))
# lines(model_out_eip50$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
# lines(model_out_eip50$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
# lines(model_out_eip50$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


# Model fitting: Vector competence ----------------------------------------



############ Vector competence from Shapiro et al. 2017 Plos Biology
# Get data
data_bc_constant <- with(data.bc.EIP, data.frame('T' = temp, 'trait' = bc)) # subset specific trait data from complete list
# data <- data.specific # assign trait data to variable 'data'

# Organize Data for JAGS
trait <- data_bc_constant$trait
N.obs <- length(trait)
temp <- data_bc_constant$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
model_out_bc <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
					 n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# # Examine model output & run diagnostics
# model_out_bc$BUGSoutput$summary[1:10,]
# mcmcplot(model_out_bc)
# 
# # Save model output 
# save(model_out_bc, file = "saved-posteriors/bc_briere_uniform.Rdata")
# 
# # Plot trait data, model mean and CIs
# plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_bc, 
#      ylab = "Vector competence", xlab = expression(paste("Temperature (",degree,"C)")))
# lines(model_out_bc$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
# lines(model_out_bc$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
# lines(model_out_bc$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


# Model fitting: Pea ------------------------------------------------------



############ Pea - Quadratic from Paaijmans 2013 Global Climate Change
# Get data
data_pea <- with(data.pea.MDR, data.frame('T' = temp, 'trait' = Pea)) # subset specific trait data from complete list
# data <- data.specific # assign trait data to variable 'data'

# Organize Data for JAGS
trait <- data_pea$trait
N.obs <- length(trait)
temp <- data_pea$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
data_pea
model_out_pea <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
					  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# # Examine model output & run diagnostics
# model_out_pea$BUGSoutput$summary[1:10,]
# mcmcplot(model_out_pea)
# 
# # Save model output 
# save(model_out_pea, file = "saved-posteriors/pea_briere_uniform_Krijn2013.Rdata")
# 
# # Plot trait data, model mean and CIs
# plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_pea, 
#      ylab = "Pea", xlab = expression(paste("Temperature (",degree,"C)")))
# lines(model_out_pea$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
# lines(model_out_pea$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
# lines(model_out_pea$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


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
data_mdr
model_out_mdr <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
					  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# # Examine model output & run diagnostics
# model_out_mdr$BUGSoutput$summary[1:10,]
# mcmcplot(model_out_mdr)
# 
# # Save model output 
# save(model_out_mdr, file = "saved-posteriors/MDR_briere_uniform_Krijn2013.Rdata")
# 
# # Plot trait data, model mean and CIs
# plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_mdr, 
#      ylab = "MDR", xlab = expression(paste("Temperature (",degree,"C)")))
# lines(model_out_mdr$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
# lines(model_out_mdr$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
# lines(model_out_mdr$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


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


#### ok since I can't get the RDatas to work, let's just export the dataframes with model outputs

## ok there are some things we need:
# 1) We need the predicted values for the trait from the fitted model, to make the graphs
# 2) We need the parameter estimates for Tmin, Topt and Tmax




# Bring all model outputs together ----------------------------------------


# 1. vector competence (bc) -------------------------------------------------------


View(b_params_bc)

predictions_bc_constant <- as.data.frame(model_out_bc$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)   ### columns are temperatures, rows are iterations
colnames(predictions_bc_constant) <- Temp.xs

predictions_bc_constant_summary <- predictions_bc_constant %>%
	mutate(iteration = rownames(.)) %>% 
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>%
	dplyr::group_by(temperature) %>%  
	summarise(`2.5%`=quantile(growth_rate, probs=0.025),
			  `97.5%`=quantile(growth_rate, probs=0.975),
			  mean = mean(growth_rate)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "bc_constant")

write_csv(predictions_bc_constant_summary, "data-processed/predictions_bc_constant_summary.csv")

topt_bc_constant <- predictions_bc_constant %>% 
	mutate(iteration = rownames(.)) %>%
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>% 
	group_by(iteration) %>% 
	top_n(n = 1, wt = growth_rate) %>% 
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(`2.5%` =quantile(temperature, probs=0.025),
			  `97.5%`=quantile(temperature, probs=0.975),
			  mean = mean(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "bc_constant")

b_params_bc_constant <- as.data.frame(model_out_bc$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>% 
	mutate(treatment = "bc_constant")

params_bc_constant_all <- bind_rows(b_params_bc_constant, topt_bc_constant)


write_csv(params_bc_constant_all, "data-processed/params_bc_constant_all.csv")

### raw data to plot
data_bc_constant_sum <- data_bc_constant %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = std.error(trait)) %>% 
	mutate(trait = "vector competence") %>% 
	mutate(treatment = "bc_constant")

write_csv(data_bc_constant_sum, "data-processed/data_bc_constant_sum.csv")



ggplot() + 
	geom_ribbon(aes(x = temperature, ymax = `97.5%`, ymin = `2.5%`), data = predictions_bc_constant_summary, fill = "#8da0cb") +
	# geom_line(aes(x = temperature, y = growth_rate, group = iteration), alpha = 0.05, size = 1, data = predictions_long2) +
	# geom_point(aes(x = T, y = trait), data = data_bite_rate_dtr12, color = "#8da0cb", alpha = 0.7, size = 2) + geom_hline(yintercept = 0) +
	geom_line(aes(x = temperature, y = mean), data = predictions_bc_constant_summary, color = "black") +
	geom_line(aes(x = temperature, y = `2.5%`), data = predictions_bc_constant_summary, color = "black", linetype = "dashed") +
	geom_line(aes(x = temperature, y = `97.5%`), data = predictions_bc_constant_summary, color = "black", linetype = "dashed") +
	geom_pointrange(aes(x = temperature, y = mean, ymin = mean -std_error, ymax = mean + std_error), data = data_bc_constant_sum) +
	ylab("Vector competence, bc") + 
	# xlab("Temperature (°C)") +
	xlab("") +
	ylim(0, 1) + xlim(0, 45) 


# 2. Bite rate ---------------------------------------------------------------

# model_out_bite_rate_constant_briere

predictions_bite_rate_constant <- as.data.frame(model_out_bite_rate_constant_briere$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)   ### columns are temperatures, rows are iterations
colnames(predictions_bite_rate_constant) <- Temp.xs

predictions_bite_rate_constant_summary <- predictions_bite_rate_constant %>%
	mutate(iteration = rownames(.)) %>% 
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>%
	dplyr::group_by(temperature) %>%  
	summarise(`2.5%`=quantile(growth_rate, probs=0.025),
			  `97.5%`=quantile(growth_rate, probs=0.975),
			  mean = mean(growth_rate)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "bite_rate_constant")

write_csv(predictions_bite_rate_constant_summary, "data-processed/predictions_bite_rate_constant_summary.csv")

topt_bite_rate_constant <- predictions_bite_rate_constant %>% 
	mutate(iteration = rownames(.)) %>%
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>% 
	group_by(iteration) %>% 
	top_n(n = 1, wt = growth_rate) %>% 
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(`2.5%` =quantile(temperature, probs=0.025),
			  `97.5%`=quantile(temperature, probs=0.975),
			  mean = mean(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "bite_rate_constant")

b_params_bite_rate_constant <- as.data.frame(model_out_bite_rate_constant_briere$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>% 
	mutate(treatment = "bite_rate_constant")

params_bite_rate_constant_all <- bind_rows(b_params_bite_rate_constant, topt_bite_rate_constant)


write_csv(params_bite_rate_constant_all, "data-processed/params_bite_rate_constant_all.csv")

### raw data to plot
data_bite_rate_constant_sum <- data_constant_bite_rate_2 %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = std.error(trait)) %>% 
	mutate(trait = "bite rate") %>% 
	mutate(treatment = "bite_rate_constant")

write_csv(data_bite_rate_constant_sum, "data-processed/data_bite_rate_constant_sum.csv")


# 3. Bite rate dtr 9 ------------------------------------------------------

data.specific_bite_rate_dtr9
model.out_bite_rate_dtr9_briere


predictions_bite_rate_dtr9 <- as.data.frame(model.out_bite_rate_dtr9_briere$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)   ### columns are temperatures, rows are iterations
colnames(predictions_bite_rate_dtr9) <- Temp.xs

predictions_bite_rate_dtr9_summary <- predictions_bite_rate_dtr9 %>%
	mutate(iteration = rownames(.)) %>% 
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>%
	dplyr::group_by(temperature) %>%  
	summarise(`2.5%`=quantile(growth_rate, probs=0.025),
			  `97.5%`=quantile(growth_rate, probs=0.975),
			  mean = mean(growth_rate)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "bite_rate_dtr9")

write_csv(predictions_bite_rate_dtr9_summary, "data-processed/predictions_bite_rate_dtr9_summary.csv")

topt_bite_rate_dtr9 <- predictions_bite_rate_dtr9 %>% 
	mutate(iteration = rownames(.)) %>%
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>% 
	group_by(iteration) %>% 
	top_n(n = 1, wt = growth_rate) %>% 
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(`2.5%` =quantile(temperature, probs=0.025),
			  `97.5%`=quantile(temperature, probs=0.975),
			  mean = mean(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "bite_rate_dtr9")

b_params_bite_rate_dtr9 <- as.data.frame(model.out_bite_rate_dtr9_briere$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>% 
	mutate(treatment = "bite_rate_dtr9")

params_bite_rate_dtr9_all <- bind_rows(b_params_bite_rate_dtr9, topt_bite_rate_dtr9)
# View(params_bite_rate_dtr9_all)

write_csv(params_bite_rate_dtr9_all, "data-processed/params_bite_rate_dtr9_all.csv")

### raw data to plot
data_bite_rate_dtr9_sum <- data.specific_bite_rate_dtr9 %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = std.error(trait)) %>% 
	mutate(trait = "bite rate") %>% 
	mutate(treatment = "bite_rate_dtr9")

write_csv(data_bite_rate_dtr9_sum, "data-processed/data_bite_rate_dtr9_sum.csv")


# 4. Bite rate dtr 12 -----------------------------------------------------

data_bite_rate_dtr12
model.out_bite_rate_dtr12_briere

predictions_bite_rate_dtr12 <- as.data.frame(model.out_bite_rate_dtr12_briere$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)   ### columns are temperatures, rows are iterations
colnames(predictions_bite_rate_dtr12) <- Temp.xs

predictions_bite_rate_dtr12_summary <- predictions_bite_rate_dtr12 %>%
	mutate(iteration = rownames(.)) %>% 
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>%
	dplyr::group_by(temperature) %>%  
	summarise(`2.5%`=quantile(growth_rate, probs=0.025),
			  `97.5%`=quantile(growth_rate, probs=0.975),
			  mean = mean(growth_rate)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "bite_rate_dtr12")

write_csv(predictions_bite_rate_dtr12_summary, "data-processed/predictions_bite_rate_dtr12_summary.csv")

topt_bite_rate_dtr12 <- predictions_bite_rate_dtr12 %>% 
	mutate(iteration = rownames(.)) %>%
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>% 
	group_by(iteration) %>% 
	top_n(n = 1, wt = growth_rate) %>% 
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(`2.5%` =quantile(temperature, probs=0.025),
			  `97.5%`=quantile(temperature, probs=0.975),
			  mean = mean(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "bite_rate_dtr12")

b_params_bite_rate_dtr12 <- as.data.frame(model.out_bite_rate_dtr12_briere$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>% 
	mutate(treatment = "bite_rate_dtr12")

params_bite_rate_dtr12_all <- bind_rows(b_params_bite_rate_dtr12, topt_bite_rate_dtr12)
# View(params_bite_rate_dtr12_all)

write_csv(params_bite_rate_dtr12_all, "data-processed/params_bite_rate_dtr12_all.csv")

### raw data to plot
data_bite_rate_dtr12_sum <- data_bite_rate_dtr12 %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = std.error(trait)) %>% 
	mutate(trait = "bite rate") %>% 
	mutate(treatment = "bite_rate_dtr12")

write_csv(data_bite_rate_dtr12_sum, "data-processed/data_bite_rate_dtr12_sum.csv")


# 5. Lifespan constant temperature ---------------------------------------
data_constant_lifespan
model_out_constant_lifespan_briere
model_out_constant_lifespan_briere

predictions_lifespan_constant <- as.data.frame(model_out_constant_lifespan_briere$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)   ### columns are temperatures, rows are iterations
colnames(predictions_lifespan_constant) <- Temp.xs

predictions_lifespan_constant_summary <- predictions_lifespan_constant %>%
	mutate(iteration = rownames(.)) %>% 
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>% 
	dplyr::group_by(temperature) %>%  
	dplyr::summarise(`2.5%`=quantile(growth_rate, probs=0.025),
					 `97.5%`=quantile(growth_rate, probs=0.975),
					 mean = mean(growth_rate)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "lifespan_constant")

write_csv(predictions_lifespan_constant_summary, "data-processed/predictions_lifespan_constant_summary.csv")

### ok let's dig in here and see if we can figure out why we are getting these weird shaped confidence intervals


predictions_lifespan_constant %>% 
	mutate(iteration = rownames(.)) %>% 
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>%
	mutate(temperature = as.numeric(temperature)) %>% 
	ggplot(aes(x = temperature, y = growth_rate, group = iteration)) + geom_line()
ggsave("figures/lifespan_constant_predictions.pdf", width = 8, height = 6)

topt_lifespan_constant <- predictions_lifespan_constant %>% 
	mutate(iteration = rownames(.)) %>%
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>% 
	group_by(iteration) %>% 
	top_n(n = 1, wt = growth_rate) %>% 
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(`2.5%` =quantile(temperature, probs=0.025),
			  `97.5%`=quantile(temperature, probs=0.975),
			  mean = mean(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "lifespan_constant")

b_params_lifespan_constant <- as.data.frame(model_out_constant_lifespan_briere$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>% 
	mutate(treatment = "lifespan_constant")

params_lifespan_constant_all <- bind_rows(b_params_lifespan_constant, topt_lifespan_constant)
# View(params_lifespan_constant_all)

write_csv(params_lifespan_constant_all, "data-processed/params_lifespan_constant_all.csv")

### raw data to plot
data_lifespan_constant_sum <- data_constant_lifespan %>% 
	dplyr::rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = std.error(trait)) %>% 
	mutate(trait = "bite rate") %>% 
	mutate(treatment = "lifespan_constant")

write_csv(data_lifespan_constant_sum, "data-processed/data_lifespan_constant_sum.csv")



# 6. Lifespan DTR9 --------------------------------------------------------
### ok we need to come back to these results because they are showing this two grouped thing
data_dtr9_lifespan
model_out_lifespan_dtr9_briere



predictions_lifespan_dtr9 <- as.data.frame(model_out_lifespan_dtr9_briere$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)   ### columns are temperatures, rows are iterations
colnames(predictions_lifespan_dtr9) <- Temp.xs

# predictions_lifespan_dtr9 %>% 
# 	mutate(iteration = rownames(.)) %>% 
# 	dplyr::select(iteration, everything()) %>% 
# 	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>%
# 	mutate(temperature = as.numeric(temperature)) %>% 
# 	ggplot(aes(x = temperature, y = growth_rate, group = iteration)) + geom_line()


predictions_lifespan_dtr9_summary <- predictions_lifespan_dtr9 %>%
	mutate(iteration = rownames(.)) %>% 
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>%
	dplyr::group_by(temperature) %>%  
	summarise(`2.5%`=quantile(growth_rate, probs=0.025),
			  `97.5%`=quantile(growth_rate, probs=0.975),
			  mean = mean(growth_rate)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "lifespan_dtr9")

write_csv(predictions_lifespan_dtr9_summary, "data-processed/predictions_lifespan_dtr9_summary.csv")

topt_lifespan_dtr9 <- predictions_lifespan_dtr9 %>% 
	mutate(iteration = rownames(.)) %>%
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>% 
	group_by(iteration) %>% 
	top_n(n = 1, wt = growth_rate) %>% 
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(`2.5%` =quantile(temperature, probs=0.025),
			  `97.5%`=quantile(temperature, probs=0.975),
			  mean = mean(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "lifespan_dtr9")

b_params_lifespan_dtr9 <- as.data.frame(model_out_lifespan_dtr9_briere$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>% 
	mutate(treatment = "lifespan_dtr9")

params_lifespan_dtr9_all <- bind_rows(b_params_lifespan_dtr9, topt_lifespan_dtr9)
# View(params_lifespan_dtr9_all)

write_csv(params_lifespan_dtr9_all, "data-processed/params_lifespan_dtr9_all.csv")

### raw data to plot
data_lifespan_dtr9_sum <- data_dtr9_lifespan %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = std.error(trait)) %>% 
	mutate(trait = "bite rate") %>% 
	mutate(treatment = "lifespan_dtr9")

write_csv(data_lifespan_dtr9_sum, "data-processed/data_lifespan_dtr9_sum.csv")



# 7. Lifespan DTR12 -------------------------------------------------------

data_dtr12_lifespan
model_out_lifespan_dtr12

predictions_lifespan_dtr12 <- as.data.frame(model_out_lifespan_dtr12$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)   ### columns are temperatures, rows are iterations
colnames(predictions_lifespan_dtr12) <- Temp.xs

predictions_lifespan_dtr12_summary <- predictions_lifespan_dtr12 %>%
	mutate(iteration = rownames(.)) %>% 
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>%
	dplyr::group_by(temperature) %>%  
	summarise(`2.5%`=quantile(growth_rate, probs=0.025),
			  `97.5%`=quantile(growth_rate, probs=0.975),
			  mean = mean(growth_rate)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "lifespan_dtr12")

write_csv(predictions_lifespan_dtr12_summary, "data-processed/predictions_lifespan_dtr12_summary.csv")

topt_lifespan_dtr12 <- predictions_lifespan_dtr12 %>% 
	mutate(iteration = rownames(.)) %>%
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>% 
	group_by(iteration) %>% 
	top_n(n = 1, wt = growth_rate) %>% 
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(`2.5%` =quantile(temperature, probs=0.025),
			  `97.5%`=quantile(temperature, probs=0.975),
			  mean = mean(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "lifespan_dtr12")

b_params_lifespan_dtr12 <- as.data.frame(model_out_lifespan_dtr12$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>% 
	mutate(treatment = "lifespan_dtr12")

params_lifespan_dtr12_all <- bind_rows(b_params_lifespan_dtr12, topt_lifespan_dtr12)
# View(params_lifespan_dtr12_all)

write_csv(params_lifespan_dtr12_all, "data-processed/params_lifespan_dtr12_all.csv")

### raw data to plot
data_lifespan_dtr12_sum <- data_dtr12_lifespan %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = std.error(trait)) %>% 
	mutate(trait = "bite rate") %>% 
	mutate(treatment = "lifespan_dtr12")

write_csv(data_lifespan_dtr12_sum, "data-processed/data_lifespan_dtr12_sum.csv")

View(data_lifespan_dtr12_sum)

# 8. Lifetime eggs constant --------------------------------------------------
## ok come back here on January 18
### where do we get this data_constant_eggs from?

data_constant_eggs
model_out_eggs_constant

predictions_eggs_constant <- as.data.frame(model_out_eggs_constant$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)   ### columns are temperatures, rows are iterations
colnames(predictions_eggs_constant) <- Temp.xs

predictions_eggs_constant_summary <- predictions_eggs_constant %>%
	mutate(iteration = rownames(.)) %>% 
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>%
	dplyr::group_by(temperature) %>%  
	summarise(`2.5%`=quantile(growth_rate, probs=0.025),
			  `97.5%`=quantile(growth_rate, probs=0.975),
			  mean = mean(growth_rate)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "eggs_constant")

write_csv(predictions_eggs_constant_summary, "data-processed/predictions_eggs_constant_summary.csv")

topt_eggs_constant <- predictions_eggs_constant %>% 
	mutate(iteration = rownames(.)) %>%
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>% 
	group_by(iteration) %>% 
	top_n(n = 1, wt = growth_rate) %>% 
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(`2.5%` =quantile(temperature, probs=0.025),
			  `97.5%`=quantile(temperature, probs=0.975),
			  mean = mean(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "eggs_constant")

b_params_eggs_constant <- as.data.frame(model_out_eggs_constant$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>% 
	mutate(treatment = "eggs_constant")

params_constant_eggs_all <- bind_rows(b_params_eggs_constant, topt_eggs_constant)
# View(params_lifespan_dtr12_all)

write_csv(params_constant_eggs_all, "data-processed/params_constant_eggs_all.csv")

### raw data to plot
data_constant_eggs_sum <- data_constant_eggs %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = std.error(trait)) %>% 
	mutate(trait = "eggs") %>% 
	mutate(treatment = "eggs_constant")

write_csv(data_constant_eggs_sum, "data-processed/data_constant_eggs_sum.csv")




# 9. Lifetime eggs DTR9 ---------------------------------------------------

data_eggs_dtr9
model_out_eggs_dtr9

predictions_eggs_dtr9 <- as.data.frame(model_out_eggs_dtr9$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)   ### columns are temperatures, rows are iterations
colnames(predictions_eggs_dtr9) <- Temp.xs

predictions_eggs_dtr9_summary <- predictions_eggs_dtr9 %>%
	mutate(iteration = rownames(.)) %>% 
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>%
	dplyr::group_by(temperature) %>%  
	summarise(`2.5%`=quantile(growth_rate, probs=0.025),
			  `97.5%`=quantile(growth_rate, probs=0.975),
			  mean = mean(growth_rate)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "eggs_dtr9")

write_csv(predictions_eggs_dtr9_summary, "data-processed/predictions_eggs_dtr9_summary.csv")

topt_eggs_dtr9 <- predictions_eggs_dtr9 %>% 
	mutate(iteration = rownames(.)) %>%
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>% 
	group_by(iteration) %>% 
	top_n(n = 1, wt = growth_rate) %>% 
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(`2.5%` =quantile(temperature, probs=0.025),
			  `97.5%`=quantile(temperature, probs=0.975),
			  mean = mean(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "eggs_dtr9")

b_params_eggs_dtr9 <- as.data.frame(model_out_eggs_dtr9$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>% 
	mutate(treatment = "eggs_dtr9")

params_dtr9_eggs_all <- bind_rows(b_params_eggs_dtr9, topt_eggs_dtr9)
# View(params_lifespan_dtr12_all)

write_csv(params_dtr9_eggs_all, "data-processed/params_dtr9_eggs_all.csv")

### raw data to plot
data_dtr9_eggs_sum <- data_eggs_dtr9 %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = std.error(trait)) %>% 
	mutate(trait = "eggs") %>% 
	mutate(treatment = "eggs_dtr9")

write_csv(data_dtr9_eggs_sum, "data-processed/data_dtr9_eggs_sum.csv")


# 10. Lifetime eggs DTR12 -------------------------------------------------

data_eggs_dtr12
model_out_eggs_dtr12


predictions_eggs_dtr12 <- as.data.frame(model_out_eggs_dtr12$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)   ### columns are temperatures, rows are iterations
colnames(predictions_eggs_dtr12) <- Temp.xs

predictions_eggs_dtr12_summary <- predictions_eggs_dtr12 %>%
	mutate(iteration = rownames(.)) %>% 
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>%
	dplyr::group_by(temperature) %>%  
	summarise(`2.5%`=quantile(growth_rate, probs=0.025),
			  `97.5%`=quantile(growth_rate, probs=0.975),
			  mean = mean(growth_rate)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "eggs_dtr12")

write_csv(predictions_eggs_dtr12_summary, "data-processed/predictions_eggs_dtr12_summary.csv")

topt_eggs_dtr12 <- predictions_eggs_dtr12 %>% 
	mutate(iteration = rownames(.)) %>%
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>% 
	group_by(iteration) %>% 
	top_n(n = 1, wt = growth_rate) %>% 
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(`2.5%` =quantile(temperature, probs=0.025),
			  `97.5%`=quantile(temperature, probs=0.975),
			  mean = mean(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "eggs_dtr12")

b_params_eggs_dtr12 <- as.data.frame(model_out_eggs_dtr12$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>% 
	mutate(treatment = "eggs_dtr12")

params_dtr12_eggs_all <- bind_rows(b_params_eggs_dtr12, topt_eggs_dtr12)
# View(params_lifespan_dtr12_all)

write_csv(params_dtr12_eggs_all, "data-processed/params_dtr12_eggs_all.csv")

### raw data to plot
data_dtr12_eggs_sum <- data_eggs_dtr12 %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = std.error(trait)) %>% 
	mutate(trait = "eggs") %>% 
	mutate(treatment = "eggs_dtr12")

write_csv(data_dtr12_eggs_sum, "data-processed/data_dtr12_eggs_sum.csv")



# 11. EIP50 ---------------------------------------------------------------

data_eip50
model_out_eip50


predictions_eip50 <- as.data.frame(model_out_eip50$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)   ### columns are temperatures, rows are iterations
colnames(predictions_eip50) <- Temp.xs

predictions_eip50_summary <- predictions_eip50 %>%
	mutate(iteration = rownames(.)) %>% 
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>%
	dplyr::group_by(temperature) %>%  
	summarise(`2.5%`=quantile(growth_rate, probs=0.025),
			  `97.5%`=quantile(growth_rate, probs=0.975),
			  mean = mean(growth_rate)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "eip50")

write_csv(predictions_eip50_summary, "data-processed/predictions_eip50_summary.csv")

topt_eip50 <- predictions_eip50 %>% 
	mutate(iteration = rownames(.)) %>%
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>% 
	group_by(iteration) %>% 
	top_n(n = 1, wt = growth_rate) %>% 
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(`2.5%` =quantile(temperature, probs=0.025),
			  `97.5%`=quantile(temperature, probs=0.975),
			  mean = mean(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "eip50")

b_params_eip50 <- as.data.frame(model_out_eip50$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>% 
	mutate(treatment = "eip50")

params_eip50_all <- bind_rows(b_params_eip50, topt_eip50)
# View(params_lifespan_dtr12_all)

write_csv(params_eip50_all, "data-processed/params_eip50_all.csv")

### raw data to plot
data_eip50_sum <- data_eip50 %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = std.error(trait)) %>% 
	mutate(trait = "eip50") %>% 
	mutate(treatment = "eip50")

write_csv(data_eip50_sum, "data-processed/data_eip50_sum.csv")





# 12. PEA -----------------------------------------------------------------

data_pea
model_out_pea

predictions_pea <- as.data.frame(model_out_pea$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)   ### columns are temperatures, rows are iterations
colnames(predictions_pea) <- Temp.xs

predictions_pea_summary <- predictions_pea %>%
	mutate(iteration = rownames(.)) %>% 
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>%
	dplyr::group_by(temperature) %>%  
	summarise(`2.5%`=quantile(growth_rate, probs=0.025),
			  `97.5%`=quantile(growth_rate, probs=0.975),
			  mean = mean(growth_rate)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "pea")

write_csv(predictions_pea_summary, "data-processed/predictions_pea_summary.csv")

topt_pea <- predictions_pea %>% 
	mutate(iteration = rownames(.)) %>%
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>% 
	group_by(iteration) %>% 
	top_n(n = 1, wt = growth_rate) %>% 
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(`2.5%` =quantile(temperature, probs=0.025),
			  `97.5%`=quantile(temperature, probs=0.975),
			  mean = mean(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "pea")

b_params_pea <- as.data.frame(model_out_pea$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>% 
	mutate(treatment = "pea")

params_pea_all <- bind_rows(b_params_pea, topt_pea)
# View(params_lifespan_dtr12_all)

write_csv(params_pea_all, "data-processed/params_pea_all.csv")

### raw data to plot
data_pea_sum <- data_pea %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = std.error(trait)) %>% 
	mutate(trait = "pea") %>% 
	mutate(treatment = "pea")

write_csv(data_pea_sum, "data-processed/data_pea_sum.csv")




# 13. MDR -----------------------------------------------------------------

data_mdr
model_out_mdr


predictions_mdr <- as.data.frame(model_out_mdr$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)   ### columns are temperatures, rows are iterations
colnames(predictions_mdr) <- Temp.xs

predictions_mdr_summary <- predictions_mdr %>%
	mutate(iteration = rownames(.)) %>% 
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>%
	dplyr::group_by(temperature) %>%  
	summarise(`2.5%`=quantile(growth_rate, probs=0.025),
			  `97.5%`=quantile(growth_rate, probs=0.975),
			  mean = mean(growth_rate)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "mdr")

write_csv(predictions_mdr_summary, "data-processed/predictions_mdr_summary.csv")

topt_mdr <- predictions_mdr %>% 
	mutate(iteration = rownames(.)) %>%
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>% 
	group_by(iteration) %>% 
	top_n(n = 1, wt = growth_rate) %>% 
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(`2.5%` =quantile(temperature, probs=0.025),
			  `97.5%`=quantile(temperature, probs=0.975),
			  mean = mean(temperature)) %>% 
	mutate(term = "Topt") %>% 
	mutate(treatment = "mdr")

b_params_mdr <- as.data.frame(model_out_mdr$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>% 
	mutate(treatment = "mdr")

params_mdr_all <- bind_rows(b_params_mdr, topt_mdr)
# View(params_lifespan_dtr12_all)

write_csv(params_mdr_all, "data-processed/params_mdr_all.csv")

### raw data to plot
data_mdr_sum <- data_mdr %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = std.error(trait)) %>% 
	mutate(trait = "mdr") %>% 
	mutate(treatment = "mdr")

write_csv(data_mdr_sum, "data-processed/data_mdr_sum.csv")
