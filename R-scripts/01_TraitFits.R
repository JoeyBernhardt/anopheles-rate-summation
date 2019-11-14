################################# Contents
# 1. Set-up + data vis 
# 2. Bayesian Settings
# 3. Bayesian Fitting
# 4. Rate Summation Calculation


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

# ###### Summarizing + plotting data - method #1 (lifespan only)
# 
# # Sumarize dataframes
# cdata.sum <- cdata %>%
# 	group_by(Treatment) %>%
# 	summarise(DTR = "constant", lifespan.mean = mean(lifespan), lifespan.SE = sd(lifespan)/sqrt(n()),
# 			  bite.rate.mean = mean(lifespan), bite.rate.SE = sd(bite.rate)/sqrt(n()),
# 			  lifetime.eggs.mean = mean(lifetime.eggs), lifetime.eggs.SE = sd(lifetime.eggs)/sqrt(n()))
# 
# f9data.sum <- f9data %>%
# 	group_by(Treatment) %>%
# 	summarise(DTR = "+/- 9", lifespan.mean = mean(lifespan), lifespan.SE = sd(lifespan)/sqrt(n()),
# 			  bite.rate.mean = mean(lifespan), bite.rate.SE = sd(bite.rate)/sqrt(n()),
# 			  lifetime.eggs.mean = mean(lifetime.eggs), lifetime.eggs.SE = sd(lifetime.eggs)/sqrt(n()))
# 
# f12data.sum <- f12data %>%
# 	group_by(Treatment) %>%
# 	summarise(DTR = "+/- 12", lifespan.mean = mean(lifespan), lifespan.SE = sd(lifespan)/sqrt(n()),
# 			  bite.rate.mean = mean(lifespan), bite.rate.SE = sd(bite.rate)/sqrt(n()),
# 			  lifetime.eggs.mean = mean(lifetime.eggs), lifetime.eggs.SE = sd(lifetime.eggs)/sqrt(n()))
# 
# # Plot the separate dataframes with manual jittering
# ggplot() + 
# 	geom_point(data = cdata.sum, aes(x = Treatment, y = lifespan.mean)) + ylim(0,55) + 
# 	geom_errorbar(data = cdata.sum, aes(x = Treatment, ymin = lifespan.mean - lifespan.SE, ymax = lifespan.mean + lifespan.SE), width = 0) +
# 	geom_point(data = f9data.sum, aes(x = Treatment + 0.2, y = lifespan.mean), color = "blue") + 
# 	geom_errorbar(data = f9data.sum, aes(x = Treatment + 0.2, ymin = lifespan.mean - lifespan.SE, ymax = lifespan.mean + lifespan.SE), width = 0, 
# 				  color = "blue") +
# 	geom_point(data = f12data.sum, aes(x = Treatment + 0.4, y = lifespan.mean), color = "dodgerblue") +
# 	geom_errorbar(data = f12data.sum, aes(x = Treatment + 0.4, ymin = lifespan.mean - lifespan.SE, ymax = lifespan.mean + lifespan.SE), width = 0, 
# 				  color = "dodgerblue")
# 
# ###### Plotting data - method #2 (aka, Marta slowly learning dplyr/ggplot stuff) (lifespan only)
# 
# # Combine into a single dataframe, add jittering, plot
# alldata.sum <- rbind(cdata.sum, f9data.sum, f12data.sum)
# alldata.sum$jitter <- c(rep(0,6), rep(0.2,5), rep(0.4,5))
# ggplot() + 
# 	geom_point(data = alldata.sum, aes(x = Treatment + jitter, y = lifespan.mean, color = DTR)) +
# 	geom_errorbar(data = alldata.sum, aes(x = Treatment + jitter, ymin = lifespan.mean - lifespan.SE, ymax = lifespan.mean + lifespan.SE, color = DTR), width = 0) +
# 	scale_color_manual(values = c("dodgerblue", "blue", "black"))

###### Summarizing + plotting data - method #3 (aka, Marta slowly learning dplyr/ggplot stuff) (all 3 traits: lifespan, biting rate, fecundity)

# Add in DTR Treatment columns, combine in a single dataframe
cdata$DTR <- "constant"
f9data$DTR <- "+/- 9"
f12data$DTR <- "+/- 12"
alldata <- rbind(cdata, f9data, f12data)

alldata <- alldata %>% mutate(TempK = Treatment + 273.15)

# Summarize dataframe, add jitter
alldata.sum.2 <- alldata %>%
	group_by(DTR, Treatment) %>%
	summarise(lifespan.mean = mean(lifespan), lifespan.SE = sd(lifespan)/sqrt(n()),
			  bite.rate.mean = mean(bite.rate), bite.rate.SE = sd(bite.rate)/sqrt(n()),
			  lifetime.eggs.mean = mean(lifetime.eggs), lifetime.eggs.SE = sd(lifetime.eggs)/sqrt(n()))
alldata.sum.2$jitter <- c(rep(0.4,5), rep(0.2,5), rep(0,6))

# Plot lifespan
p.lf <- ggplot() + 
	geom_point(data = alldata.sum.2, aes(x = Treatment + jitter, y = lifespan.mean, color = DTR)) +
	geom_errorbar(data = alldata.sum.2, aes(x = Treatment + jitter, ymin = lifespan.mean - lifespan.SE, ymax = lifespan.mean + lifespan.SE, color = DTR), width = 0) +
	scale_color_manual(values = c("dodgerblue", "blue", "black")) + theme(legend.position = c(0.6,0.8)) +
	ylab("Lifespan (days)") + xlab("Temperature (°C)") 

# Plot biting rate
p.br <- ggplot() + 
	geom_point(data = alldata.sum.2, aes(x = Treatment + jitter, y = bite.rate.mean, color = DTR)) +
	geom_errorbar(data = alldata.sum.2, aes(x = Treatment + jitter, ymin = bite.rate.mean - bite.rate.SE, ymax = bite.rate.mean + bite.rate.SE, color = DTR), width = 0) +
	scale_color_manual(values = c("dodgerblue", "blue", "black")) + theme(legend.position = "none") +
	ylab(expression(paste("Biting rate (days"^-1,")"))) + xlab("Temperature (°C)") 

# Plot fecundity
p.f <- ggplot() + 
	geom_point(data = alldata.sum.2, aes(x = Treatment + jitter, y = lifetime.eggs.mean, color = DTR)) +
	geom_errorbar(data = alldata.sum.2, aes(x = Treatment + jitter, ymin = lifetime.eggs.mean - lifetime.eggs.SE, ymax = lifetime.eggs.mean + lifetime.eggs.SE, color = DTR), width = 0) +
	scale_color_manual(values = c("dodgerblue", "blue", "black")) + theme(legend.position = "none") +
	ylab("Lifetime fecundity (eggs)") + xlab("Temperature (°C)") 

plot_grid(p.lf, p.br, p.f, nrow = 1)
ggsave("figures/SummarizedTraitData.pdf")


################################# 2. Bayesian Settings

# Save model files for JAGS to read - only need to do this once (or after editing the model)

# Model with log-normal distributed observations around trait mean - dlnorm() always positive, so no need to constrain values
sink("norberg_lnorm.txt") 
cat("model{ 
	## Priors
	cf.z ~ dunif(0, 100)
	cf.w ~ dunif(5, 30)
	cf.a ~ dunif(-0.02, 20)
	cf.b ~ dunif(0, 1)
	cf.tau ~ dnorm(1000, 1/(250^2))
	cf.sigma <- sqrt(1/cf.tau)
	# cf.tau <- 1/(cf.sigma^2)
	# cf.sigma ~ dexp(0.0001)

	## Likelihood 
	for(i in 1:N.obs){
	trait.mu[i] <- cf.a * exp(cf.b * temp[i]) * (1-((temp[i]-cf.z)/(cf.w/2))^2) 
	trait[i] ~ dlnorm(trait.mu[i], cf.tau) } 
	
	} # close model ",fill=T)
sink()

# Model with normally distributed observations around trait mean - dlnorm() can go negative, so need to constrain values
# The T() function truncates the normal distribution at zero (with +Inf as the upper limit because none is specified)
sink("norberg_tnorm.txt") 
cat("model{ 
	## Priors
	cf.z ~ dunif(0, 100)
	cf.w ~ dunif(5, 30)
	cf.a ~ dunif(-0.02, 20)
	cf.b ~ dunif(0, 1)
	cf.tau ~ dnorm(1000, 1/(250^2))
	cf.sigma <- sqrt(1/cf.tau)
	# cf.tau <- 1/(cf.sigma^2)
	# cf.sigma ~ dexp(0.0001)
	
	## Likelihood 
	for(i in 1:N.obs){
	trait.temporary[i] <- cf.a * exp(cf.b * temp[i]) * (1-((temp[i]-cf.z)/(cf.w/2))^2)
	trait.mu[i] <- trait.temporary[i] * ((cf.a * exp(cf.b * temp[i]) * (1-((temp[i]-cf.z)/(cf.w/2))^2)) > 0) 
	trait[i] ~ dnorm(trait.mu[i], cf.tau)T(0,) } 
	
	} # close model ",fill=T)
sink()

sink("quad_T.txt")
cat("
	model{
	
	## Priors
	cf.q ~ dunif(0, 1)
	cf.T0 ~ dunif(0, 24)
	cf.Tm ~ dunif(25, 50)
	cf.sigma ~ dunif(0, 1000)
	cf.tau <- 1 / (cf.sigma * cf.sigma)
	
	## Likelihood
	for(i in 1:N.obs){
	trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
	trait[i] ~ dnorm(trait.mu[i], cf.tau)T(0,) 
	}
	
	} # close model ",fill=T)
sink()

sink("quad.txt")
cat("
	model{
	
	## Priors
	cf.q ~ dunif(0, 1)
	cf.T0 ~ dunif(0, 24)
	cf.Tm ~ dunif(25, 50)
	cf.sigma ~ dunif(0, 1000)
	cf.tau <- 1 / (cf.sigma * cf.sigma)
	
	## Likelihood
	for(i in 1:N.obs){
	trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
	trait[i] ~ dnorm(trait.mu[i], cf.tau)
	}
	
	} # close model ",fill=T)
sink()


### NOTE THAT FOR SHARPE-SCHOOLFIEND MODEL TEMPERATURE MUST BE IN KELVIN!!!
sink("SS_T.txt")
cat("
	model{
	
	## Priors
	cf.c ~ dunif(0, 50)
	cf.Topt ~ dunif(274, 313)
	cf.ED ~ dunif(0, 10)
	cf.sigma ~ dunif(0, 1000)
	cf.tau <- 1 / (cf.sigma * cf.sigma)
	
	## Likelihood
	for(i in 1:N.obs){
	trait.mu[i] <- cf.c * 10^11 * exp(-0.6/(8.62*10^-5*temp[i])) / (1 + exp(-1/(8.62*10^-5*temp[i]) * (cf.ED - temp[i]*(cf.ED/cf.Topt + 8.62*10^-5*log(0.6/(cf.ED-0.6)) ) ) ) ) 
	trait[i] ~ dnorm(trait.mu[i], cf.tau)T(0,) 
	}
	
	} # close model ",fill=T)
sink()

### NOTE THAT FOR SHARPE-SCHOOLFIEND MODEL TEMPERATURE MUST BE IN KELVIN!!!
sink("SS.txt")
cat("
	model{
	
	## Priors
	cf.c ~ dunif(0, 50)
	cf.Topt ~ dunif(274, 313)
	cf.ED ~ dunif(0, 10)
	cf.sigma ~ dunif(0, 1000)
	cf.tau <- 1 / (cf.sigma * cf.sigma)
	
	## Likelihood
	for(i in 1:N.obs){
	trait.mu[i] <- cf.c * 10^11 * exp(-0.6/(8.62*10^-5*temp[i])) / (1 + exp(-1/(8.62*10^-5*temp[i]) * (cf.ED - temp[i]*(cf.ED/cf.Topt + 8.62*10^-5*log(0.6/(cf.ED-0.6)) ) ) ) ) 
	trait[i] ~ dnorm(trait.mu[i], cf.tau)
	}
	
	} # close model ",fill=T)
sink()


### Norberg function parameters and inits
# ## parameters to estimate
# parameters <- c("cf.z", "cf.w", "cf.a", "cf.b", "cf.tau")
# 
# # Initial values for the parameters
# inits<-function(){list(
# 	cf.z = 25,
# 	cf.w = 22,
# 	cf.a = 12,
# 	cf.b = 0.01,
# 	cf.tau = 100
# )}

### Modified Sharpe-Schoolfield function parameters and inits
## parameters to estimate
parameters <- c("cf.c", "cf.Topt", "cf.ED", "cf.tau")

# Initial values for the parameters
inits<-function(){list(
	cf.c = 8,
	cf.Topt = 293,
	cf.ED = 1.5,
	cf.sigma = rlnorm(1))}

# ### Quadratic/Briere function parameters and inits
# 
# ##### Parameters to Estimate
# parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma")
# 
# ### Quadratic function parameters and inits
# inits<-function(){list(
# 	cf.q = 0.01,
# 	cf.Tm = 35,
# 	cf.T0 = 5,
# 	cf.sigma = rlnorm(1))}




# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 25000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains

# Temperature sequence for derived quantity calculations
Temp.xs <- seq(5, 45, 1)
Temp.xs.K <- Temp.xs + 273.15

# Function to calculate traits as a function of temperature for a Norberg function - all MCMC steps or just summary
trait.trajs.norberg = function(jags.fitted.model, Temp.xs, summary){
	
	# Extract vectors of function coeffiencients
	mod <- jags.fitted.model
	cf.z <- mod$BUGSoutput$sims.list$cf.z
	cf.w <- mod$BUGSoutput$sims.list$cf.w
	cf.a <- mod$BUGSoutput$sims.list$cf.a
	cf.b <- mod$BUGSoutput$sims.list$cf.b
	
	# Loop through MCMC chains and calculate trait trajectories
	traj.df <- matrix(nrow = length(Temp.xs), ncol = nrow(cf.z))
	for(i in 1:nrow(cf.z)){
		traj.temporary <- cf.a[i] * exp(cf.b[i] * Temp.xs) * (1-((Temp.xs-cf.z[i])/(cf.w[i]/2))^2)
		traj.df[,i] <- traj.temporary * (traj.temporary > 0)
	}
	
	# Calculate mean, median, and CIs
	traj.CIs.df <- data.frame(temp = Temp.xs,
							  mean = numeric(length(Temp.xs)), 
							  median = numeric(length(Temp.xs)), 
							  upper = numeric(length(Temp.xs)), 
							  lower = numeric(length(Temp.xs)))
	for(j in 1:length(Temp.xs)){
		traj.CIs.df$mean[j] <- mean(traj.df[j,])
		traj.CIs.df$median[j] <- mean(traj.df[j,])
		traj.CIs.df$upper[j] <- quantile(traj.df[j,], 0.975)
		traj.CIs.df$lower[j] <- quantile(traj.df[j,], 0.025)
	}
	
	# Return either full set of trajectories for every step in the MCMC chain or just the summary statistics
	ifelse(summary == TRUE, return(traj.CIs.df), return(traj.df))
	
}

# Function to calculate traits as a function of temperature for a Norberg function - all MCMC steps or just summary
trait.trajs.quad = function(jags.fitted.model, Temp.xs, summary){
	
	# Extract vectors of function coeffiencients
	mod <- jags.fitted.model
	cf.q <- mod$BUGSoutput$sims.list$cf.q
	cf.Tm <- mod$BUGSoutput$sims.list$cf.Tm
	cf.T0 <- mod$BUGSoutput$sims.list$cf.T0
	
	# Loop through MCMC chains and calculate trait trajectories
	traj.df <- matrix(nrow = length(Temp.xs), ncol = nrow(cf.q))
	for(i in 1:nrow(cf.q)){
		traj.temporary <- -1 * cf.q[i] * (Temp.xs - cf.T0[i]) * (Temp.xs - cf.Tm[i])
		traj.df[,i] <- traj.temporary * (traj.temporary > 0)
	}
	
	# Calculate mean, median, and CIs
	traj.CIs.df <- data.frame(temp = Temp.xs,
							  mean = numeric(length(Temp.xs)), 
							  median = numeric(length(Temp.xs)), 
							  upper = numeric(length(Temp.xs)), 
							  lower = numeric(length(Temp.xs)))
	for(j in 1:length(Temp.xs)){
		traj.CIs.df$mean[j] <- mean(traj.df[j,])
		traj.CIs.df$median[j] <- mean(traj.df[j,])
		traj.CIs.df$upper[j] <- quantile(traj.df[j,], 0.975)
		traj.CIs.df$lower[j] <- quantile(traj.df[j,], 0.025)
	}
	
	# Return either full set of trajectories for every step in the MCMC chain or just the summary statistics
	ifelse(summary == TRUE, return(traj.CIs.df), return(traj.df))
	
}


# Function to calculate traits as a function of temperature for a Norberg function - all MCMC steps or just summary
trait.trajs.SS = function(jags.fitted.model, Temp.xs, summary){
	
	# Extract vectors of function coeffiencients
	mod <- jags.fitted.model
	cf.c <- mod$BUGSoutput$sims.list$cf.c
	cf.Topt <- mod$BUGSoutput$sims.list$cf.Topt
	cf.ED <- mod$BUGSoutput$sims.list$cf.ED
	E <- 0.6
	k <- 8.62*10^-5
	
	# Loop through MCMC chains and calculate trait trajectories
	traj.df <- matrix(nrow = length(Temp.xs), ncol = nrow(cf.c))
	for(i in 1:nrow(cf.c)){
		traj.temporary <- cf.c[i] * 10^11 * exp(-E/(k*Temp.xs)) / (1 + exp(-1/(k*Temp.xs) * (cf.ED[i] - Temp.xs*(cf.ED[i]/cf.Topt[i] + k*log(E/(cf.ED[i]-E)) ) ) ) ) 
		traj.df[,i] <- traj.temporary * (traj.temporary > 0)
	}
	
	# Calculate mean, median, and CIs
	traj.CIs.df <- data.frame(temp.K = Temp.xs,
							  temp.C = Temp.xs - 273.15,
							  mean = numeric(length(Temp.xs)), 
							  median = numeric(length(Temp.xs)), 
							  upper = numeric(length(Temp.xs)), 
							  lower = numeric(length(Temp.xs)))
	for(j in 1:length(Temp.xs)){
		traj.CIs.df$mean[j] <- mean(traj.df[j,])
		traj.CIs.df$median[j] <- mean(traj.df[j,])
		traj.CIs.df$upper[j] <- quantile(traj.df[j,], 0.975)
		traj.CIs.df$lower[j] <- quantile(traj.df[j,], 0.025)
	}
	
	# Return either full set of trajectories for every step in the MCMC chain or just the summary statistics
	ifelse(summary == TRUE, return(traj.CIs.df), return(traj.df))
	
}


################################# 3. Bayesian Fitting

###### Lifespan - constant temperature

# Pull out data
data <- alldata %>%
	filter(DTR == "constant") %>%
	rename(trait = lifespan,
		   T = TempK)

# Save as vectors for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle all data in a list for JAGS
#jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)

# Fit the model
lf.c.nor <- jags(data=jag.data, inits=inits, parameters.to.save=parameters,
				   model.file="norberg_tnorm.txt", n.thin=nt, n.chains=nc, n.burnin=nb,
				   n.iter=ni, DIC=T, working.directory=getwd())

lf.c.quad <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
			 model.file="quad.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
			 n.iter=ni, DIC=T, working.directory=getwd())	

lf.c.quad_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
				  model.file="quad_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
				  n.iter=ni, DIC=T, working.directory=getwd())	

lf.c.SS_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
					model.file="SS_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
					n.iter=ni, DIC=T, working.directory=getwd())	

lf.c.SS <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
				model.file="SS.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
				n.iter=ni, DIC=T, working.directory=getwd())	

mcmcplot(lf.c.SS)
lf.c.SS$BUGSoutput$summary
lf.c.SS$BUGSoutput$DIC
lf.c.SS_T$BUGSoutput$DIC

# Calculate summary statistics of trait trajectories for plotting
#lf.c.nor.summary <- trait.trajs.norberg(lf.c.nor, Temp.xs, summary = TRUE)
lf.c.quad.summary <- trait.trajs.quad(lf.c.quad, Temp.xs, summary = TRUE)
lf.c.quad_T.summary <- trait.trajs.quad(lf.c.quad_T, Temp.xs, summary = TRUE)
lf.c.SS_T.summary <- trait.trajs.SS(lf.c.SS_T, Temp.xs.K, summary = TRUE)
lf.c.SS.summary <- trait.trajs.SS(lf.c.SS, Temp.xs.K, summary = TRUE)


lf.c.SS_T.summary[1:10, 1:10]

# Plot data and trait fit
p.lf.fit.c <- alldata.sum.2 %>%
	filter(DTR == "constant") %>%
	ggplot() + 
	geom_point(aes(x = Treatment, y = lifespan.mean), color = "blue", size = 3) +
	geom_errorbar(aes(x = Treatment, ymin = lifespan.mean - lifespan.SE, ymax = lifespan.mean + lifespan.SE), width = 0, color = "blue", size = 1) +
	geom_line(data = lf.c.SS.summary, aes(x = temp.C, y = mean)) +
	geom_line(data = lf.c.SS.summary, aes(x = temp.C, y = upper), linetype = "dashed") +
	geom_line(data = lf.c.SS.summary, aes(x = temp.C, y = lower), linetype = "dashed") +
	#geom_line(data = lf.c.SS_T.summary, aes(x = temp.C, y = mean), color = "red") +
	#geom_line(data = lf.c.SS_T.summary, aes(x = temp.C, y = upper), linetype = "dashed", color = "red") +
	#geom_line(data = lf.c.SS_T.summary, aes(x = temp.C, y = lower), linetype = "dashed", color = "red") +
	ylab("Lifespan (days)") + xlab("Temperature (°C)") 
p.lf.fit.c 
ggsave("figures/lifespan_SSfit.pdf")


# Plot data and trait fit
p.lf.fit.c <- alldata.sum.2 %>%
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


lf.fit.nls = nlsLM(trait ~ quad(T, c, T0, Tm), data = data,
				start=list(c = 0.1, T0 = 10, Tm = 40))
coef.lf = coef(lf.fit.nls)
plot(lifespan.mean ~ Treatment, data = cdata.sum, ylab = "Lifespan", xlab = "Temperature", ylim = c(0,45))
plot.fit(T, "quad", coef.lf)


################################# 4. Rate Summation Calculation