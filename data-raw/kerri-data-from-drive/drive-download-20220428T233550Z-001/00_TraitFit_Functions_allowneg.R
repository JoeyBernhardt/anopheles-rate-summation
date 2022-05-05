############Models and trajectories for functions that can become negative
################################# Functions to calculate trajectories
#Set working directory
mainDir = "C:/Users/Kerri/Desktop/Chapter2 InProgress"
setwd(mainDir)


# Function to calculate traits as a function of temperature for a quadratic function - all MCMC steps or just summary
trait.trajs.neg.quad = function(jags.fitted.model, Temp.xs, summary){
  
  # Extract vectors of function coeffiencients
  mod <- jags.fitted.model
  cf.q <- mod$BUGSoutput$sims.list$cf.q
  cf.Tm <- mod$BUGSoutput$sims.list$cf.Tm
  cf.T0 <- mod$BUGSoutput$sims.list$cf.T0
  
  # Loop through MCMC steps and calculate trait trajectories
  traj.df <- matrix(nrow = length(Temp.xs), ncol = nrow(cf.q))
  for(i in 1:nrow(cf.q)){
    traj.temporary <- -1 * cf.q[i] * (Temp.xs - cf.T0[i]) * (Temp.xs - cf.Tm[i])
    traj.df[,i] <- traj.temporary
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
trait.trajs.neg.briere = function(jags.fitted.model, Temp.xs, summary){
  
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
    traj.df[,i] <- traj.temporary
    
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




sink("quad_neg.txt")
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
	trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm)
	trait[i] ~ dnorm(trait.mu[i], cf.tau)T(0,)
	}
	
	} # close model ",fill=T)
sink()

sink("briere_neg_test.txt")
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
	trait.mu[i] <- ifelse(i<=cf.Tm, cf.q * temp[i] * (temp[i] - cf.T0) * sqrt((cf.Tm - temp[i]) * (cf.Tm > temp[i])),-2*(sqrt(((temp[i]-cf.Tm)/((cf.Tm-cf.T0)/2))*(cf.Tm < temp[i]))))
	trait[i] ~ dnorm(trait.mu[i], cf.tau)T(0,)
	}
	
	} # close model ",fill=T)
sink()


sink("SS_T_KMadjust.txt")
cat("
model{
	
	## Priors
	cf.c ~ dunif(0, 50)
	cf.Topt ~ dunif(274, 313)
	cf.ED ~ dunif(1, 10)
	cf.sigma ~ dunif(0, 1000)
	cf.tau <- 1 / (cf.sigma * cf.sigma)
	
	## Likelihood
	for(i in 1:N.obs){
	trait.mu[i] <- cf.c * 10^11 * exp(-0.6/(8.62*10^-5*temp[i])) / (1 + exp(-1/(8.62*10^-5*temp[i]) * (cf.ED - temp[i]*(cf.ED/cf.Topt + 8.62*10^-5*log(0.6/(cf.ED-0.6)) ) ) ) ) 
	trait[i] ~ dnorm(trait.mu[i], cf.tau)
	}
	
	} # close model ",fill=T)
sink()


# Function to calculate traits as a function of temperature for a briere function - all MCMC steps or just summary

trait.trajs.neg.test.briere = function(jags.fitted.model, Temp.xs, summary){
 
  # Extract vectors of function coeffiencients
  mod <- jags.fitted.model
  cf.q <- mod$BUGSoutput$sims.list$cf.q
  cf.Tm <- mod$BUGSoutput$sims.list$cf.Tm
  cf.T0 <- mod$BUGSoutput$sims.list$cf.T0
  
  # Loop through MCMC steps and calculate trait trajectories
  traj.df <- matrix(nrow = length(Temp.xs), ncol = nrow(cf.q)) #initializing the matrix
  traj.temporary <- matrix(nrow = length(Temp.xs), ncol = nrow(cf.q)) #initializing the matrix
  
  for(i in 1:nrow(cf.q)){
    for(j in 1: length(Temp.xs)){#I think I need to have another for statment here j representing the Temp.xs values....
    # This takes care of NAs during the calculation
    #Conditional statment added to change the equation used depending on what T is in relation to the cf.Tm
   traj.temporary[j,i] <- cf.q[i] * Temp.xs[j] * (Temp.xs[j] - cf.T0[i]) * sqrt((cf.Tm[i] - Temp.xs[j])*(cf.Tm[i] > Temp.xs[j]))
   if(traj.temporary[j,i] == 0){
      if(Temp.xs[j]> 20){
        traj.temporary[j,i] <- -2*(sqrt(((Temp.xs[j]-cf.Tm[i])/((cf.Tm[i]-cf.T0[i])/2))))
       }
    }
    traj.df[j,i] <- traj.temporary[j,i] #i correspondes to the full simulation output....
    
    }#end of j loop
  }#end of i loop
  # Calculate mean, median, and CIs
  traj.CIs.df <- data.frame(temp = Temp.xs,
                            mean = numeric(length(Temp.xs)),
                            median = numeric(length(Temp.xs)),
                            upper = numeric(length(Temp.xs)),
                            lower = numeric(length(Temp.xs)))
  
  #need to make anything negative below the T0 limit 0
  #need to make anything negative above the Tmax limit 0 for these calculations...
  for(k in 1:length(Temp.xs)){
    traj.CIs.df$mean[k] <- mean(traj.df[k,])
    traj.CIs.df$median[k] <- mean(traj.df[k,])
    traj.CIs.df$upper[k] <- quantile(traj.df[k,], 0.975)
    traj.CIs.df$lower[k] <- quantile(traj.df[k,], 0.025)
  }
  
  # Return either full set of trajectories for every step in the MCMC chain or just the summary statistics
  ifelse(summary == TRUE, return(traj.CIs.df), return(traj.df))
  
}


###trait trajectories for the modified SS (KMadjusted) outputs
# Function to calculate traits as a function of temperature for a briere function - all MCMC steps or just summary
#Generated 7/30/20

trait.trajs.SS_T_KMadjust = function(jags.fitted.model, Temp.xs, summary){
  
  # Extract vectors of function coeffiencients
  mod <- jags.fitted.model
  cf.c <- mod$BUGSoutput$sims.list$cf.c
  cf.Topt <- mod$BUGSoutput$sims.list$cf.Topt
  cf.ED <- mod$BUGSoutput$sims.list$cf.ED
  
    # Loop through MCMC steps and calculate trait trajectories
    traj.df <- matrix(nrow = length(Temp.xs), ncol = nrow(cf.c))
    for(i in 1:nrow(cf.c)){
      
      # Takes care of NAs after the fact
      # traj.temporary <- cf.q[i] * Temp.xs * (Temp.xs - cf.T0[i]) * sqrt(cf.Tm[i] - Temp.xs) 
      # traj.temporary[is.na(traj.temporary)] <- 0
      # traj.df[,i] <- traj.temporary * (traj.temporary > 0)
      
      traj.df[,i]  <- cf.c[i] * 10^11 * exp(-0.6/(8.62*10^-5*Temp.xs)) / (1 + exp(-1/(8.62*10^-5*Temp.xs) * (cf.ED[i] - Temp.xs*(cf.ED[i]/cf.Topt[i] + 8.62*10^-5*log(0.6/(cf.ED[i]-0.6)) ) ) ) ) 
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
  
  


#####################################
#####################################
#Generated 7/30/20
trait.trajs.Norberg_KMadjust = function(jags.fitted.model, Temp.xs, summary){
  
  # Extract vectors of function coeffiencients
  mod <- jags.fitted.model
  cf.z <- mod$BUGSoutput$sims.list$cf.z
  cf.w <- mod$BUGSoutput$sims.list$cf.w
  cf.a <- mod$BUGSoutput$sims.list$cf.a
  cf.b <- mod$BUGSoutput$sims.list$cf.b
  
  
  # Loop through MCMC steps and calculate trait trajectories
  traj.df <- matrix(nrow = length(Temp.xs), ncol = nrow(cf.z))
  for(i in 1:nrow(cf.z)){
    
    # Takes care of NAs after the fact
    # traj.temporary <- cf.q[i] * Temp.xs * (Temp.xs - cf.T0[i]) * sqrt(cf.Tm[i] - Temp.xs) 
    # traj.temporary[is.na(traj.temporary)] <- 0
    # traj.df[,i] <- traj.temporary * (traj.temporary > 0)
    
    traj.df[,i]  <- cf.a[i] * exp(cf.b[i] * Temp.xs) * (1-((Temp.xs-cf.z[i])/(cf.w[i]/2))^2)  
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

sink("norberg_KMadjust2.txt")
cat("
model{ 
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
    trait[i] ~ dnorm(trait.mu[i], cf.tau)T(0,) } 
	} # close model ",fill=T)
sink()