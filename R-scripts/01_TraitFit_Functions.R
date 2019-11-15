################################# Functions to calculate trajectories

# Function to calculate traits as a function of temperature for a Norberg function - all MCMC steps or just summary
trait.trajs.norberg = function(jags.fitted.model, Temp.xs, summary){
	
	# Extract vectors of function coeffiencients
	mod <- jags.fitted.model
	cf.z <- mod$BUGSoutput$sims.list$cf.z
	cf.w <- mod$BUGSoutput$sims.list$cf.w
	cf.a <- mod$BUGSoutput$sims.list$cf.a
	cf.b <- mod$BUGSoutput$sims.list$cf.b
	
	# Loop through MCMC steps and calculate trait trajectories
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

# Function to calculate traits as a function of temperature for a quadratic function - all MCMC steps or just summary
trait.trajs.quad = function(jags.fitted.model, Temp.xs, summary){
	
	# Extract vectors of function coeffiencients
	mod <- jags.fitted.model
	cf.q <- mod$BUGSoutput$sims.list$cf.q
	cf.Tm <- mod$BUGSoutput$sims.list$cf.Tm
	cf.T0 <- mod$BUGSoutput$sims.list$cf.T0
	
	# Loop through MCMC steps and calculate trait trajectories
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

# Function to calculate traits as a function of temperature for a briere function - all MCMC steps or just summary
trait.trajs.briere = function(jags.fitted.model, Temp.xs, summary){
	
	# Extract vectors of function coeffiencients
	mod <- jags.fitted.model
	cf.q <- mod$BUGSoutput$sims.list$cf.q
	cf.Tm <- mod$BUGSoutput$sims.list$cf.Tm
	cf.T0 <- mod$BUGSoutput$sims.list$cf.T0
	
	# Loop through MCMC steps and calculate trait trajectories
	traj.df <- matrix(nrow = length(Temp.xs), ncol = nrow(cf.q))
	for(i in 1:nrow(cf.q)){
		
		# Takes care of NAs after the fact
		# traj.temporary <- cf.q[i] * Temp.xs * (Temp.xs - cf.T0[i]) * sqrt(cf.Tm[i] - Temp.xs) 
		# traj.temporary[is.na(traj.temporary)] <- 0
		# traj.df[,i] <- traj.temporary * (traj.temporary > 0)
		
		# This takes care of NAs during the calculation
		traj.temporary <- cf.q[i] * Temp.xs * (Temp.xs - cf.T0[i]) * sqrt((cf.Tm[i] - Temp.xs)*(cf.Tm[i] > Temp.xs))
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
	
	# Loop through MCMC steps and calculate trait trajectories
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


################################# Model Files - does not work with dplyr loaded?

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