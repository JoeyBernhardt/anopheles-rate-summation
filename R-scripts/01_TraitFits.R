################################# Contents
# 1. Set-up + data vis 
# 2. Bayesian Settings (code to write model files is in 01_TraitFit_Functions.R)
# 3. Bayesian Fitting - lifespan
# 4. Bayesian Fitting - biting rate
# 5. Bayesian Fitting - fecundity

################################# 1. Set-up + data vis

###### Load packages
library(tidyverse)
library(R2jags)
library(coda)
library(cowplot)
library(mcmcplots)

######  Load raw trait data
cdata <- read_csv("data-raw/constant.individual.trait.csv")
f9data <- read_csv("data-raw/fluc9.individual.trait.csv")
f12data <- read_csv("data-raw/fluc12.individual.trait.csv")

###### Summarize + plot data for all 3 traits: lifespan, biting rate, fecundity

# Add in DTR Treatment columns, combine in a single dataframe
cdata$DTR <- "constant"
f9data$DTR <- "+/- 9"
f12data$DTR <- "+/- 12"
alldata <- rbind(cdata, f9data, f12data)

alldata <- alldata %>% mutate(TempK = Treatment + 273.15)

# Summarize dataframe, add jitter
alldata.sum <- alldata %>%
	group_by(DTR, Treatment) %>%
	summarise(lifespan.mean = mean(lifespan), lifespan.SE = sd(lifespan)/sqrt(n()),
			  bite.rate.mean = mean(bite.rate), bite.rate.SE = sd(bite.rate)/sqrt(n()),
			  lifetime.eggs.mean = mean(lifetime.eggs), lifetime.eggs.SE = sd(lifetime.eggs)/sqrt(n()))
alldata.sum$jitter <- c(rep(0.4,5), rep(0.2,5), rep(0,6))

# Plot lifespan
p.lf <- ggplot() + 
	geom_point(data = alldata.sum, aes(x = Treatment + jitter, y = lifespan.mean, color = DTR)) +
	geom_errorbar(data = alldata.sum, aes(x = Treatment + jitter, ymin = lifespan.mean - lifespan.SE, ymax = lifespan.mean + lifespan.SE, color = DTR), width = 0) +
	scale_color_manual(values = c("dodgerblue", "blue", "black")) + theme(legend.position = c(0.6,0.8)) +
	ylab("Lifespan (days)") + xlab("Temperature (°C)") 

# Plot biting rate
p.br <- ggplot() + 
	geom_point(data = alldata.sum, aes(x = Treatment + jitter, y = bite.rate.mean, color = DTR)) +
	geom_errorbar(data = alldata.sum, aes(x = Treatment + jitter, ymin = bite.rate.mean - bite.rate.SE, ymax = bite.rate.mean + bite.rate.SE, color = DTR), width = 0) +
	scale_color_manual(values = c("dodgerblue", "blue", "black")) + theme(legend.position = "none") +
	ylab(expression(paste("Biting rate (days"^-1,")"))) + xlab("Temperature (°C)") 

# Plot fecundity
p.f <- ggplot() + 
	geom_point(data = alldata.sum, aes(x = Treatment + jitter, y = lifetime.eggs.mean, color = DTR)) +
	geom_errorbar(data = alldata.sum, aes(x = Treatment + jitter, ymin = lifetime.eggs.mean - lifetime.eggs.SE, ymax = lifetime.eggs.mean + lifetime.eggs.SE, color = DTR), width = 0) +
	scale_color_manual(values = c("dodgerblue", "blue", "black")) + theme(legend.position = "none") +
	ylab("Lifetime fecundity (eggs)") + xlab("Temperature (°C)") 

plot_grid(p.lf, p.br, p.f, nrow = 1)
ggsave("figures/SummarizedTraitData.pdf")


################################# 2. Bayesian Settings

### MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 25000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains

### Temperature sequence for derived quantity calculations
Temp.xs <- seq(5, 45, 1)
Temp.xs.K <- Temp.xs + 273.15 # Needed for Sharpe-Schoolfied type models


### Norberg function parameters and inits
# parameters to estimate
parameters <- c("cf.z", "cf.w", "cf.a", "cf.b", "cf.tau")

# Initial values for the parameters
inits<-function(){list(
	cf.z = 25,
	cf.w = 22,
	cf.a = 12,
	cf.b = 0.01,
	cf.tau = 100)}


# ### Modified Sharpe-Schoolfield function parameters and inits
# # parameters to estimate
# parameters <- c("cf.c", "cf.Topt", "cf.ED", "cf.tau")
# 
# # Initial values for the parameters
# inits<-function(){list(
# 	cf.c = 8,
# 	cf.Topt = 293,
# 	cf.ED = 1.5,
# 	cf.sigma = rlnorm(1))}


### Quadratic/Briere function parameters and inits
# Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma")

# Initial values for the parameters
inits<-function(){list(
	cf.q = 0.01,
	cf.Tm = 35,
	cf.T0 = 5,
	cf.sigma = rlnorm(1))}


################################# 3. Bayesian Fitting - lifespan

###### Lifespan - constant temperature

# Pull out data - *** Make sure Temperature units are correct (C or K as needed) ***
data <- alldata %>%
	filter(DTR == "constant") %>%
	rename(trait = lifespan,
		   T = TempK)

# Save as vectors for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)

# Fit the model
# lf.c.nor <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
# 				   model.file="norberg_tnorm.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
# 				   n.iter=ni, DIC=T, working.directory=getwd())
# 
# lf.c.quad <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
# 			 model.file="quad.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
# 			 n.iter=ni, DIC=T, working.directory=getwd())	
# 
# lf.c.quad_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
# 				  model.file="quad_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
# 				  n.iter=ni, DIC=T, working.directory=getwd())	

lf.c.SS_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
					model.file="SS_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
					n.iter=ni, DIC=T, working.directory=getwd())	
lf.c.SS_T$BUGSoutput$summary

lf.c.SS <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
				model.file="SS.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
				n.iter=ni, DIC=T, working.directory=getwd())	
lf.c.SS$BUGSoutput$summary

mcmcplot(lf.c.SS)
lf.c.SS$BUGSoutput$DIC
lf.c.SS_T$BUGSoutput$DIC

# Calculate summary statistics of trait trajectories for plotting
#lf.c.nor.summary <- trait.trajs.norberg(lf.c.nor, Temp.xs, summary = TRUE)
lf.c.quad.summary <- trait.trajs.quad(lf.c.quad, Temp.xs, summary = TRUE)
lf.c.quad_T.summary <- trait.trajs.quad(lf.c.quad_T, Temp.xs, summary = TRUE)
lf.c.SS_T.summary <- trait.trajs.SS(lf.c.SS_T, Temp.xs.K, summary = TRUE)
lf.c.SS.summary <- trait.trajs.SS(lf.c.SS, Temp.xs.K, summary = TRUE)


# Plot data and trait fit
p.lf.fit.c <- alldata.sum %>%
	filter(DTR == "constant") %>%
	ggplot() + 
	geom_point(aes(x = Treatment, y = lifespan.mean), color = "blue", size = 3) +
	geom_errorbar(aes(x = Treatment, ymin = lifespan.mean - lifespan.SE, ymax = lifespan.mean + lifespan.SE), width = 0, color = "blue", size = 1) +
	geom_line(data = lf.c.SS.summary, aes(x = temp.C, y = mean)) +
	geom_line(data = lf.c.SS.summary, aes(x = temp.C, y = upper), linetype = "dashed") +
	geom_line(data = lf.c.SS.summary, aes(x = temp.C, y = lower), linetype = "dashed") +
	geom_line(data = lf.c.SS_T.summary, aes(x = temp.C, y = mean), color = "red") +
	geom_line(data = lf.c.SS_T.summary, aes(x = temp.C, y = upper), linetype = "dashed", color = "red") +
	geom_line(data = lf.c.SS_T.summary, aes(x = temp.C, y = lower), linetype = "dashed", color = "red") +
	ylab("Lifespan (days)") + xlab("Temperature (°C)") 
p.lf.fit.c 
ggsave("figures/lifespan_SSfit.pdf")


# Plot data and trait fit
p.lf.fit.c <- alldata.sum %>%
	filter(DTR == "constant") %>%
	ggplot() + 
	geom_point(aes(x = Treatment, y = lifespan.mean), color = "blue", size = 3) +
	geom_errorbar(aes(x = Treatment, ymin = lifespan.mean - lifespan.SE, ymax = lifespan.mean + lifespan.SE), width = 0, color = "blue", size = 1) +
	geom_line(data = lf.c.quad.summary, aes(x = temp, y = mean)) +
	geom_line(data = lf.c.quad.summary, aes(x = temp, y = upper), linetype = "dashed") +
	geom_line(data = lf.c.quad.summary, aes(x = temp, y = lower), linetype = "dashed") +
	geom_line(data = lf.c.quad_T.summary, aes(x = temp, y = mean), color = "red") +
	geom_line(data = lf.c.quad_T.summary, aes(x = temp, y = upper), linetype = "dashed", color = "red") +
	geom_line(data = lf.c.quad_T.summary, aes(x = temp, y = lower), linetype = "dashed", color = "red") +
	ylab("Lifespan (days)") + xlab("Temperature (°C)") 
p.lf.fit.c + annotate(geom = "text", x = 40, y = 40, label = "Model: No truncation", size = 6, color = "black") +
	annotate(geom = "text", x = 40, y = 35, label = "Model: Truncation", size = 6, color = "red")
ggsave("figures/lifespan_norbergfit.pdf")
ggsave("figures/lifespan_quadfit_truncvsnotrunc.pdf")



################################# 4. Bayesian Fitting - biting rate

###### bite rate - constant temperature

# Pull out data - *** Make sure Temperature units are correct (C or K as needed) ***
data <- alldata %>%
	filter(DTR == "constant") %>%
	rename(trait = bite.rate,
		   T = Treatment)

# Save as vectors for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)

# Fit the model
br.c.briere <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
				  model.file="briere.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
				  n.iter=ni, DIC=T, working.directory=getwd())
br.c.briere$BUGSoutput$summary

br.c.briere_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
					model.file="briere_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
					n.iter=ni, DIC=T, working.directory=getwd())
br.c.briere$BUGSoutput$summary

br.c.nor <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
				 model.file="norberg_tnorm.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
				 n.iter=ni, DIC=T, working.directory=getwd())
br.c.nor$BUGSoutput$summary

br.c.quad_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
					model.file="quad_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
					n.iter=ni, DIC=T, working.directory=getwd())

br.c.quad <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
					model.file="quad.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
					n.iter=ni, DIC=T, working.directory=getwd())


# br.c.SS_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
# 				  model.file="SS_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
# 				  n.iter=ni, DIC=T, working.directory=getwd())	
# 
# br.c.SS <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
# 				model.file="SS.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
# 				n.iter=ni, DIC=T, working.directory=getwd())	

mcmcplot(br.c.briere)


# Calculate summary statistics of trait trajectories for plotting
br.c.briere.summary <- trait.trajs.briere(br.c.briere, Temp.xs, summary = TRUE)
br.c.briere_T.summary <- trait.trajs.briere(br.c.briere_T, Temp.xs, summary = TRUE)

br.c.nor.summary <- trait.trajs.norberg(br.c.nor, Temp.xs, summary = TRUE)
br.c.quad_T.summary <- trait.trajs.quad(br.c.quad_T, Temp.xs, summary = TRUE)

# Plot data and trait fit
p.br.fit.c <- alldata.sum %>%
	filter(DTR == "constant") %>%
	ggplot() + 
	geom_point(aes(x = Treatment, y = bite.rate.mean), color = "blue", size = 3) +
	geom_errorbar(aes(x = Treatment, ymin = bite.rate.mean - bite.rate.SE, ymax = bite.rate.mean + bite.rate.SE), width = 0, color = "blue", size = 1) +
	# geom_line(data = br.c.quad_T.summary, aes(x = temp, y = mean)) +
	# geom_line(data = br.c.quad_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
	# geom_line(data = br.c.quad_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
	geom_line(data = br.c.briere.summary, aes(x = temp, y = mean)) +
	geom_line(data = br.c.briere.summary, aes(x = temp, y = upper), linetype = "dashed") +
	geom_line(data = br.c.briere.summary, aes(x = temp, y = lower), linetype = "dashed") +
	geom_line(data = br.c.briere_T.summary, aes(x = temp, y = mean), color = "red") +
	geom_line(data = br.c.briere_T.summary, aes(x = temp, y = upper), linetype = "dashed", color = "red") +
	geom_line(data = br.c.briere_T.summary, aes(x = temp, y = lower), linetype = "dashed", color = "red") +
	ylab("Lifespan (days)") + xlab("Temperature (°C)") 
p.br.fit.c + annotate(geom = "text", x = 11, y = 0.4, label = "Model: No truncation", size = 6, color = "black") +
	annotate(geom = "text", x = 11, y = 0.5, label = "Model: Truncation", size = 6, color = "red")
ggsave("figures/biterate_brierefits.pdf")
