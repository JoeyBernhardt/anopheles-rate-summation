


library('R2jags') # Fits Bayesian models
library('mcmcplots') # Diagnostic plots for fits
library('MASS') # Fits distributions for informative priors
library(plyr) # Slices and dices data
library(plotrix) # For standard error function
library(tidyverse)



data.fluc9 <- read.csv("data-raw/fluc9.individual.trait.csv")
names(data.fluc9)[names(data.fluc9) == 'Treatment'] <- 'temp'

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
model.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
				  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
# model.out$BUGSoutput$summary[1:10,]
# mcmcplot(model.out)

# Save model output 
# save(model.out, file = "saved posteriors/dtr9_lifespan_quadratic_uniform.Rdata")

# Plot trait data, model mean and CIs


plot(trait ~ T, xlim = c(0, 45), data = data.specific, 
	 ylab = "Lifespan", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


