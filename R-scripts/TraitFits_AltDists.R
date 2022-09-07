#### Marta Shocket, May 2022
#### Count models (Poisson and negative binomial) based on: https://georgederpa.github.io/teaching/countModels.html

##########
###### 1. Set up
##########

# Load packages
library(tidyverse)
library(R2jags)
library(coda)
library(mcmcplots)

# Set working directory
#mainDir = "C:/Users/Kerri/Desktop/Chapter2 InProgress"
mainDir = "~/Dropbox/Research Savage Lab/anopheles-rate-summation"
setwd(mainDir)

#  Load raw trait data
cdata <- data.frame(read.csv("data-raw/constant.individual.trait.csv"))
f9data <- data.frame(read.csv("data-raw/fluc9.individual.trait.csv"))
f12data <- data.frame(read.csv("data-raw/fluc12.individual.trait.csv"))

# Calculate blocks means for plotting
# Constant temps
cdata.blockmeans <- cdata %>%
	dplyr::group_by(Treatment,Block)%>%
	dplyr::summarise(lifespan = mean(lifespan),
					 bite.rate = mean(bite.rate),
					 lifetime.eggs = mean(lifetime.eggs)) %>%
	ungroup()

# Fluctuating 9 temps
f9data.blockmeans <- f9data %>%
	dplyr::group_by(Treatment,Block)%>%
	dplyr::summarise(lifespan = mean(lifespan),
					 bite.rate = mean(bite.rate),
					 lifetime.eggs = mean(lifetime.eggs)) %>%
	ungroup()

# Fluctuating 12 temps
f12data.blockmeans <- f12data %>%
	dplyr::group_by(Treatment,Block)%>%
	dplyr::summarise(lifespan = mean(lifespan),
					 bite.rate = mean(bite.rate),
					 lifetime.eggs = mean(lifetime.eggs)) %>%
	ungroup()


##########
###### 2. Alternate Models
##########

################### Poisson distributed Briere

sink("briere_pois.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 20)
    cf.Tm ~ dunif(28, 45)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait[i] ~ dpois(trait.mu[i])
    trait.mu[i] <- cf.q * temp[i] * (temp[i] - cf.T0) * sqrt((cf.Tm - temp[i]) * (cf.Tm > temp[i])) * (cf.T0 < temp[i])
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i])
    }
    
    } # close model
    ",fill=T)
sink()

################### Negative Binomial distributed Briere

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

################### Gamma distributed Briere

sink("briere_gamma.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 20)
    cf.Tm ~ dunif(28, 45)
    cf.ra ~ dunif(0.0001, 1)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait[i] ~ dgamma(sh[i], cf.ra)
    sh[i] <- cf.ra * trait.mu[i]
    trait.mu[i] <- cf.q * temp[i] * (temp[i] - cf.T0) * sqrt((cf.Tm - temp[i]) * (cf.Tm > temp[i])) * (cf.T0 < temp[i])
    
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i])
    }
    
    } # close model
    ",fill=T)
sink()

################### Truncated normal Briere

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


################### Poisson distributed Quadratic

sink("quad_pois.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 20)
    cf.Tm ~ dunif(28, 45)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait[i] ~ dpois(trait.mu[i])
    trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])
    }
    
    } # close model
    ",fill=T)
sink()

################### Negative Binomial distributed Quadratic

sink("quad_negbin.txt")
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
    trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
    
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])
    }
    
    } # close model
    ",fill=T)
sink()


################### Gamma distributed Quadratic

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


##########
###### 3. Shared Settings
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


##########
###### 3. Lifetime Egg Fits - Norm, Poisson, Negative Binomial Distributions
##########

################### Data for Lifetime eggs - individual

# Save as vectors for JAGS
data <- cdata
data <- f9data
data <- f12data
data <- cdata[which(cdata$Treatment <= 32),]

trait <- data$lifetime.eggs
N.obs <- length(trait)
temp <- data$Treatment

# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

################### Data for Lifetime eggs - block means (normal distribution only)

# Save as vectors for JAGS
data <- cdata.blockmeans
trait <- data$lifetime.eggs
N.obs <- length(trait)
temp <- data$Treatment

# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

################### Normal distribution

# inits Function
inits<-function(){list(
	cf.q = 0.01,
	cf.Tm = 35,
	cf.T0 = 5,
	cf.sigma = rlnorm(1))}

# Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")

# Fit the model - individual data, constant temp
eggs.norm.c <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
					model.file="briere.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
					n.iter=ni, DIC=T, working.directory=getwd())

# Fit the model - block means data, constant temp
eggs.norm.bm.c <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
					model.file="briere.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
					n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
eggs.norm.c$BUGSoutput$summary[1:5,]
eggs.norm.bm.c$BUGSoutput$summary[1:5,]

################### Poisson distribution

##### inits Function
inits <- function(){list(
	cf.q = 0.01,
	cf.Tm = 40,
	cf.T0 = 5)}

##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm", "z.trait.mu.pred")

# Fit the model
eggs.pois.c <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
					model.file="briere_pois.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
					n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
eggs.pois.c$BUGSoutput$summary[1:5,]

################### Negative binomial distribution

##### inits Function
inits <- function(){list(
	cf.q = 0.01,
	cf.Tm = 40,
	cf.T0 = 5,
	cf.r = 1)}

##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm", "cf.r", "z.trait.mu.pred")

# Fit the model - constant temps
eggs.negbin.c <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
					model.file="briere_negbin.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
					n.iter=ni, DIC=T, working.directory=getwd())

# Fit the model - fluctuating 9 temps
eggs.negbin.f9 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
					  model.file="briere_negbin.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
					  n.iter=ni, DIC=T, working.directory=getwd())	

# Fit the model - fluctuatin 12 temps
eggs.negbin.f12 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
					  model.file="briere_negbin.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
					  n.iter=ni, DIC=T, working.directory=getwd())	

# Fit the model - constant temps no 36 degree data
eggs.negbin.cno36 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
					  model.file="briere_negbin.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
					  n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
eggs.negbin.c$BUGSoutput$summary[1:5,]
eggs.negbin.f9$BUGSoutput$summary[1:5,]
eggs.negbin.f12$BUGSoutput$summary[1:5,]
mcmcplot(eggs.negbin.c)
mcmcplot(eggs.negbin.f9)
mcmcplot(eggs.negbin.f12)
	
################### Gamma distribution - DOES NOT WORK - "Error in node cf.q - Slicer stuck at value with infinite density"

##### inits Function
inits <- function(){list(
	cf.q = 0.01,
	cf.Tm = 40,
	cf.T0 = 5,
	cf.ra = 0.003)}

##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm", "cf.ra", "z.trait.mu.pred")

# Fit the model
eggs.gamma.c <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
					  model.file="briere_gamma.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
					  n.iter=ni, DIC=T, working.directory=getwd())	

# Examine Output
eggs.gamma.c$BUGSoutput$summary[1:5,]

fitdist(subset(data$lifetime.eggs, data$Treatment == 36), distr = "gamma", method = "mme")
fitdist(data$lifespan, distr = "gamma", method = "mme")

################### Plot model output

par(mfrow = c(2,2))

# Plot individual data
plot(lifetime.eggs ~ jitter(Treatment, 0.5), xlim = c(0, 45), ylim = c(0, 1400), data = cdata,
     ylab = "Lifetime eggs", xlab = expression(paste("Temperature (",degree,"C)")))

# Plot block means
plot(lifetime.eggs ~ Treatment, xlim = c(0, 45), ylim = c(0, 500), data = cdata.blockmeans,
	 ylab = "Lifetime eggs", xlab = expression(paste("Temperature (",degree,"C)")))

lines(eggs.norm.bm.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "darkgrey")
lines(eggs.norm.bm.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "darkgrey")
lines(eggs.norm.bm.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "darkgrey")

lines(eggs.norm.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(eggs.norm.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(eggs.norm.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

lines(eggs.pois.c$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(eggs.pois.c$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(eggs.pois.c$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")

lines(eggs.negbin.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "blue")
lines(eggs.negbin.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "blue")
lines(eggs.negbin.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "blue")

legend("topleft", bty = "n", legend = c("normal dist - ind. data", "poisson dist", "negative binomial dist"), lwd = 1, col = c("black", "red", "blue"))
legend("topleft", bty = "n", legend = c("poisson dist", "negative binomial dist"), lwd = 1, col = c("red", "blue"))
legend("topleft", bty = "n", legend = c("normal dist - block means", "negative binomial dist"), lwd = 1, col = c("darkgrey", "blue"))
legend("to10pleft", bty = "n", legend = c("normal dist - ind. data", "normal dist - block means"), lwd = 1, col = c("black", "darkgrey"))

# Plot block means
plot(lifetime.eggs ~ Treatment, xlim = c(0, 45), ylim = c(0, 500), data = cdata.blockmeans, pch = 19, col = "plum3",
	 ylab = "Lifetime eggs", xlab = expression(paste("Temperature (",degree,"C)")))
points(lifetime.eggs ~ Treatment, data = f9data.blockmeans, pch = 19, col = "cyan3")
points(lifetime.eggs ~ Treatment, data = f12data.blockmeans, pch = 19, col = "chartreuse3")

lines(eggs.negbin.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "plum3")
lines(eggs.negbin.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "plum3")
lines(eggs.negbin.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "plum3")

lines(eggs.negbin.f9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "cyan3")
lines(eggs.negbin.f9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "cyan3")
lines(eggs.negbin.f9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "cyan3")

lines(eggs.negbin.f12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "chartreuse3")
lines(eggs.negbin.f12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "chartreuse3")
lines(eggs.negbin.f12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "chartreuse3")

lines(eggs.negbin.cno36$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "orchid4")
lines(eggs.negbin.cno36$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "orchid4")
lines(eggs.negbin.cno36$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "orchid4")

legend("topleft", bty = "n", legend = c("constant", "DTR 9", "DTR 12"), lwd = 1, pch = 19, col = c("plum3", "cyan3", "chartreuse3"))
legend("topleft", bty = "n", legend = c("constant", "constant - no 36"), lwd = 1, pch = 19, col = c("plum3", "orchid4"))

##########
###### 4. Lifespan Fits - Norm, Poisson, Negative Binomial Distributions
##########

################### Data for Lifespan - individual

# Save as vectors for JAGS
data <- cdata
data <- f9data
data <- f12data
data <- cdata[which(cdata$Treatment <= 32),]

trait <- data$lifespan
N.obs <- length(trait)
temp <- data$Treatment

# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

################### Data for Lifespan - block means (normal distribution only)

# Save as vectors for JAGS
data <- cdata.blockmeans
trait <- data$lifespan
N.obs <- length(trait)
temp <- data$Treatment

# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

################### Normal distribution

# inits Function
inits<-function(){list(
	cf.q = 0.01,
	cf.Tm = 35,
	cf.T0 = 5,
	cf.sigma = rlnorm(1))}

# Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")

# Fit the model - individual data
lf.norm.c <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
					model.file="quad.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
					n.iter=ni, DIC=T, working.directory=getwd())

# Fit the model - block means data
lf.norm.bm.c <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
					   model.file="quad.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
					   n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
lf.norm.c$BUGSoutput$summary[1:5,]
lf.norm.bm.c$BUGSoutput$summary[1:5,]

################### Poisson distribution

##### inits Function
inits <- function(){list(
	cf.q = 0.01,
	cf.Tm = 40,
	cf.T0 = 5)}

##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm", "z.trait.mu.pred")

# Fit the model
lf.pois.c <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
					model.file="quad_pois.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
					n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
lf.pois.c$BUGSoutput$summary[1:5,]

################### Negative binomial distribution

##### inits Function
inits <- function(){list(
	cf.q = 0.01,
	cf.Tm = 40,
	cf.T0 = 5,
	cf.r = 1)}

##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm", "cf.r", "z.trait.mu.pred")

# Fit the model
lf.negbin.c <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
					  model.file="quad_negbin.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
					  n.iter=ni, DIC=T, working.directory=getwd())	

# Fit the model
lf.negbin.f9 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
					model.file="quad_negbin.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
					n.iter=ni, DIC=T, working.directory=getwd())	

# Fit the model
lf.negbin.f12 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
					model.file="quad_negbin.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
					n.iter=ni, DIC=T, working.directory=getwd())	

# Fit the model
lf.negbin.cno36 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
					model.file="quad_negbin.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
					n.iter=ni, DIC=T, working.directory=getwd())	

# Examine Output
lf.negbin.c$BUGSoutput$summary[1:5,]
mcmcplot(lf.negbin.c)

lf.negbin.c$BUGSoutput$DIC
lf.pois.c$BUGSoutput$DIC
lf.norm.c$BUGSoutput$DIC


################### Gamma distribution

##### inits Function
inits <- function(){list(
	cf.q = 0.01,
	cf.Tm = 40,
	cf.T0 = 5,
	cf.ra = 0.003)}

##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm", "cf.ra", "z.trait.mu.pred")

# Fit the model
lf.gamma.c <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
					 model.file="quad_gamma.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
					 n.iter=ni, DIC=T, working.directory=getwd())

# Fit the model
lf.gamma.f9 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
				   model.file="quad_gamma.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
				   n.iter=ni, DIC=T, working.directory=getwd())

# Fit the model
lf.gamma.f12 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
				   model.file="quad_gamma.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
				   n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
lf.gamma.c$BUGSoutput$summary[1:5,]
lf.gamma.f9$BUGSoutput$summary[1:5,]
lf.gamma.f12$BUGSoutput$summary[1:5,]


################### Plot model output

par(mfrow = c(2,2))

# Plot individual data
plot(lifespan ~ jitter(Treatment, 0.5), xlim = c(0, 45), ylim = c(0, 70), data = cdata,
	 ylab = "Lifespan", xlab = expression(paste("Temperature (",degree,"C)")))

# Plot block means
plot(lifespan ~ Treatment, xlim = c(0, 45), ylim = c(0, 50), data = cdata.blockmeans,
	 ylab = "Lifespan", xlab = expression(paste("Temperature (",degree,"C)")))

lines(lf.norm.bm.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "darkgrey")
lines(lf.norm.bm.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "darkgrey")
lines(lf.norm.bm.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "darkgrey")

lines(lf.norm.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(lf.norm.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(lf.norm.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

lines(lf.pois.c$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(lf.pois.c$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(lf.pois.c$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")

lines(lf.negbin.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "blue")
lines(lf.negbin.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "blue")
lines(lf.negbin.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "blue")

lines(lf.gamma.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "purple")
lines(lf.gamma.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "purple")
lines(lf.gamma.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "purple")

legend("topleft", bty = "n", legend = c("normal dist - ind. data", "poisson dist", "negative binomial dist"), lwd = 1, col = c("black", "red", "blue"))
legend("topleft", bty = "n", legend = c("poisson dist", "negative binomial dist"), lwd = 1, col = c("red", "blue"))
legend("topleft", bty = "n", legend = c("normal dist - block means", "negative binomial dist"), lwd = 1, col = c("darkgrey", "blue"))
legend("topleft", bty = "n", legend = c("normal dist - ind. data", "normal dist - block means"), lwd = 1, col = c("black", "darkgrey"))

par(mfrow = c(1,2))

# Plot block means
plot(lifespan ~ Treatment, xlim = c(0, 45), ylim = c(0, 60), data = cdata.blockmeans, pch = 19, col = "plum3",
	 ylab = "Lifespan", xlab = expression(paste("Temperature (",degree,"C)"), main = "Negative binomial"))
points(lifespan ~ Treatment, data = f9data.blockmeans, pch = 19, col = "cyan3")
points(lifespan ~ Treatment, data = f12data.blockmeans, pch = 19, col = "chartreuse3")

lines(lf.negbin.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "plum3")
lines(lf.negbin.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "plum3")
lines(lf.negbin.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "plum3")

lines(lf.negbin.f9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "cyan3")
lines(lf.negbin.f9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "cyan3")
lines(lf.negbin.f9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "cyan3")

lines(lf.negbin.f12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "chartreuse3")
lines(lf.negbin.f12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "chartreuse3")
lines(lf.negbin.f12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "chartreuse3")

#lines(lf.negbin.cno36$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "orchid4")
#lines(lf.negbin.cno36$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "orchid4")
#lines(lf.negbin.cno36$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "orchid4")

legend("topleft", bty = "n", legend = c("constant", "DTR 9", "DTR 12"), lwd = 1, pch = 19, col = c("plum3", "cyan3", "chartreuse3"))
#legend("topleft", bty = "n", legend = c("constant", "constant - no 36"), lwd = 1, pch = 19, col = c("plum3", "orchid4"))

##### Gamma dists

# Plot block means
plot(lifespan ~ Treatment, xlim = c(0, 45), ylim = c(0, 60), data = cdata.blockmeans, pch = 19, col = "plum3",
	 ylab = "Lifespan", xlab = expression(paste("Temperature (",degree,"C)"), main = "Gamma"))
points(lifespan ~ Treatment, data = f9data.blockmeans, pch = 19, col = "cyan3")
points(lifespan ~ Treatment, data = f12data.blockmeans, pch = 19, col = "chartreuse3")

lines(lf.gamma.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "plum3")
lines(lf.gamma.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "plum3")
lines(lf.gamma.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "plum3")

lines(lf.gamma.f9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "cyan3")
lines(lf.gamma.f9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "cyan3")
lines(lf.gamma.f9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "cyan3")

lines(lf.gamma.f12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "chartreuse3")
lines(lf.gamma.f12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "chartreuse3")
lines(lf.gamma.f12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "chartreuse3")

#lines(lf.gamma.cno36$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "orchid4")
#lines(lf.gamma.cno36$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "orchid4")
#lines(lf.gamma.cno36$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "orchid4")

legend("topleft", bty = "n", legend = c("constant", "DTR 9", "DTR 12"), lwd = 1, pch = 19, col = c("plum3", "cyan3", "chartreuse3"))
#legend("topleft", bty = "n", legend = c("constant", "constant - no 36"), lwd = 1, pch = 19, col = c("plum3", "orchid4"))


##########
###### 5. Biting rate Fits - Norm, Poisson, Negative Binomial Distributions
##########

################### Data for Lifetime eggs - individual

# Save as vectors for JAGS
data <- cdata
data <- f9data
data <- f12data
data <- cdata[which(cdata$Treatment <= 32),]

trait <- data$bite.rate
N.obs <- length(trait)
temp <- data$Treatment

# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

################### Truncated normal distribution

# inits Function
inits<-function(){list(
	cf.q = 0.01,
	cf.Tm = 35,
	cf.T0 = 5,
	cf.sigma = rlnorm(1))}

# Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")

# Fit the model - individual data
br.tnorm.c <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
				  model.file="briere_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
				  n.iter=ni, DIC=T, working.directory=getwd())

br.tnorm.f9 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
				  model.file="briere_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
				  n.iter=ni, DIC=T, working.directory=getwd())

br.tnorm.f12 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
				  model.file="briere_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
				  n.iter=ni, DIC=T, working.directory=getwd())

# Fit the model - block means data
br.norm.bm.c <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
					 model.file="briere.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
					 n.iter=ni, DIC=T, working.directory=getwd())


# Examine Output
br.tnorm.c$BUGSoutput$summary[1:5,]
br.norm.bm.c$BUGSoutput$summary[1:5,]


################### Plot model output

par(mfrow = c(2,2))

# Plot individual data
plot(bite.rate ~ jitter(Treatment, 0.5), xlim = c(0, 45), ylim = c(0, 1), data = cdata,
	 ylab = "Lifespan", xlab = expression(paste("Temperature (",degree,"C)")))

# Plot block means
plot(bite.rate ~ Treatment, xlim = c(0, 45), ylim = c(0, 0.6), data = cdata.blockmeans, pch = 19, col = "plum3",
	 ylab = "Lifespan", xlab = expression(paste("Temperature (",degree,"C)"), main = "Gamma"))
points(bite.rate ~ Treatment, data = f9data.blockmeans, pch = 19, col = "cyan3")
points(bite.rate ~ Treatment, data = f12data.blockmeans, pch = 19, col = "chartreuse3")

lines(br.tnorm.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "plum3")
lines(br.tnorm.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "plum3")
lines(br.tnorm.c$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "plum3")

lines(br.tnorm.f9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "cyan3")
lines(br.tnorm.f9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "cyan3")
lines(br.tnorm.f9$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "cyan3")

lines(br.tnorm.f12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "chartreuse3")
lines(br.tnorm.f12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "chartreuse3")
lines(br.tnorm.f12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "chartreuse3")
