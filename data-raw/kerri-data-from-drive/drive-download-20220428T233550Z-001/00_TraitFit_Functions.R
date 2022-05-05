#Set working directory
mainDir = "C:/Users/Kerri/Desktop/Chapter2 InProgress"
setwd(mainDir)


################################# Functions to calculate trajectories

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




# Save model files for JAGS to read - only need to do this once (or after editing the model)

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
	
	} # close model ",fill=T)
sink()

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
	
	} # close model ",fill=T)
sink()

