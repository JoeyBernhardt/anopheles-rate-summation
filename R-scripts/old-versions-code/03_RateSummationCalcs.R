
################################# 1. Get temperature sequences

##### Get the incubator temperature programs for the DTR treatments

# Read in Logan-Parton caculations for 0 < T < 40 for DTR = 9 and DTR = 12
tempdata <- read_csv("data-raw/temp.programs.continuous.csv")

# Select the columns for the treatments in the experiment: T = 16, 20, 24, 28, 32
col.list <- c(18, 22, 26, 30, 34, 64, 68, 72, 76, 80)
tempprofs <- tempdata[ ,col.list]

# Convert into Kelvin
tempprofs.K <- tempprofs + 273.15

##### Load the saved rjags model objects



################################# 4. Rate Summation Calculation - lifespan

# Function for the modified SS function
modSS.calc <- function(t, c, Topt, ED){
	
	k <- 8.62*10^-5 # Boltzman constant in eV/Kelvin
	E <- 0.6
	
	calc <- c * 10^11 * exp(-E/(k*t)) / (1 + exp(-1/(k*t) * (ED - t*(ED/Topt + k*log(E/(ED-E)) ) ) ) ) 
	calc # return
}

# Rate summation for a single set of parameters (posterior distribution means)
lf.hourly.calcs <- modSS.calc(tempprofs.K, 10.79757, 294.5079, 2.082657)
lf.daily.calcs <- colMeans(lf.hourly.calcs)

# Pull out full set of MCMC parameter values
modSS_params <- data.frame(cf.c = lf.c.SS$BUGSoutput$sims.list$cf.c, 
						   cf.Topt = lf.c.SS$BUGSoutput$sims.list$cf.Topt, 
						   cf.ED = lf.c.SS$BUGSoutput$sims.list$cf.ED)

# Function to do rate summation for full set of MCMC parameter values
modSS_RateSummationCalc = function(temp.DF, param.DF){
	
	# Create output dataframe
	DF.out <- data.frame(matrix(nrow = nrow(param.DF), ncol = ncol(temp.DF)))
	colnames(DF.out) <- colnames(temp.DF)
	
	# Loop through MCMC rows
	for(i in 1:nrow(param.DF)){
		
		# Calculate hourly rates, then take column means for daily values
		hourly.calcs <- modSS.calc(temp.DF, param.DF$cf.c[i],  param.DF$cf.Topt[i],  param.DF$cf.ED[i]) # this is a dataframe same dimensions as temp.DF
		DF.out[i, ] <- colMeans(hourly.calcs) # this is a vector with # elements = # col in temp.DF
		
	}
	
	DF.out # Return output dataframe
}


lf.modSS.RS <- modSS_RateSummationCalc(tempprofs.K, modSS_params)

lf.modSS.RS.long <- data.frame(Temp = c(rep(16, 7500), rep(20, 7500), rep(24, 7500), rep(28, 7500), rep(32, 7500), 
										rep(16, 7500), rep(20, 7500), rep(24, 7500), rep(28, 7500), rep(32, 7500)),
							   DTR = c(rep("+/- 9", 7500*5), rep("+/- 12", 7500*5)),
							   lifespan = c(lf.modSS.RS[,1], lf.modSS.RS[,2], lf.modSS.RS[,3], lf.modSS.RS[,4], lf.modSS.RS[,5],
							   			 lf.modSS.RS[,6], lf.modSS.RS[,7], lf.modSS.RS[,8], lf.modSS.RS[,9], lf.modSS.RS[,10]))

# Make summary df for plotting
lf.modSS.RS.summary <- lf.modSS.RS.long %>%
	group_by(DTR, Temp) %>%
	summarise(median = median(lifespan), lower = quantile(lifespan, 0.025), upper = quantile(lifespan, 0.975))

# Plot
p.lf.RS <- ggplot() + 
	geom_point(data = alldata.sum.2, aes(x = Treatment, y = lifespan.mean, color = DTR), size = 3) + ylim(7,50) + xlim(12,38.6) + 
	scale_color_manual(values = c("dodgerblue", "mediumblue", "black")) + theme(legend.position = c(0.6,0.8)) +
	geom_line(data = lf.c.SS.summary, aes(x = temp.C, y = mean)) +
	geom_line(data = lf.c.SS.summary, aes(x = temp.C, y = upper), linetype = "dashed") +
	geom_line(data = lf.c.SS.summary, aes(x = temp.C, y = lower), linetype = "dashed") +
	geom_point(data = lf.modSS.RS.summary, aes(x = Temp + 0.4, y = median, color = DTR), size = 3) +
	geom_errorbar(data = lf.modSS.RS.summary, aes(x = Temp + 0.4, ymin = lower, ymax = upper, color = DTR, width = 0.5)) +
	ylab("Lifespan (days)") + xlab("Temperature (°C)") 
p.lf.RS

# Add in additional data observations to make full plot
lf.modSS.RS.summary <- lf.modSS.RS.summary %>%
	mutate(type = "RS.calc")
alldata.summary <- data.frame(DTR = alldata.sum.2$DTR,
							  Temp = alldata.sum.2$Treatment - 0.3,
							  median = alldata.sum.2$lifespan.mean,
							  upper = rep(NA,nrow(alldata.sum.2)),
							  lower = rep(NA,nrow(alldata.sum.2)),
							  type = "data")
plot.df <- rbind(alldata.summary, as.data.frame(lf.modSS.RS.summary))

# Plot
p.lf.RS.2 <- plot.df %>%
	ggplot() + 
	geom_point(aes(x = Temp, y = median, color = DTR, shape = type), size = 3) + ylim(7,50) + xlim(12,38.6) + 
	scale_color_manual(values = c("dodgerblue", "mediumblue", "black")) + theme(legend.position = c(0.6,0.8)) +
	geom_errorbar(aes(x = Temp, ymin = lower, ymax = upper, color = DTR, width = 0)) +
	geom_line(data = lf.c.SS.summary, aes(x = temp.C, y = mean)) +
	geom_line(data = lf.c.SS.summary, aes(x = temp.C, y = upper), linetype = "dashed") +
	geom_line(data = lf.c.SS.summary, aes(x = temp.C, y = lower), linetype = "dashed") +
	ylab("Lifespan (days)") + xlab("Temperature (°C)") 
p.lf.RS.2


################################# 4. Rate Summation Calculation - bite rate

# Function for the modified SS function
briere.calc <- function(t, q, T0, Tm){
	
	calc <- q * t * (t - T0) * sqrt((Tm - t) * (Tm > t)) * (T0 < t)
	calc # return
	
}

briere.calc(Temp.xs, .00016, 1.6, 42)

# Rate summation for a single set of parameters (posterior distribution means)
lf.hourly.calcs <- briere.calc(tempprofs, 0.0001599283, 1.565493, 42.13702)
lf.daily.calcs <- colMeans(lf.hourly.calcs)

# Pull out full set of MCMC parameter values
briere.params <- data.frame(cf.q = br.c.briere$BUGSoutput$sims.list$cf.q, 
						   cf.T0 = br.c.briere$BUGSoutput$sims.list$cf.T0, 
						   cf.Tm = br.c.briere$BUGSoutput$sims.list$cf.Tm)

# Function to do rate summation for full set of MCMC parameter values
briere_RateSummationCalc = function(temp.DF, param.DF){
	
	# Create output dataframe
	DF.out <- data.frame(matrix(nrow = nrow(param.DF), ncol = ncol(temp.DF)))
	colnames(DF.out) <- colnames(temp.DF)
	
	# Loop through MCMC rows
	for(i in 1:nrow(param.DF)){
		
		# Calculate hourly rates, then take column means for daily values
		hourly.calcs <- briere.calc(temp.DF, param.DF$cf.q[i],  param.DF$cf.T0[i],  param.DF$cf.Tm[i]) # this is a dataframe same dimensions as temp.DF
		DF.out[i, ] <- colMeans(hourly.calcs) # this is a vector with # elements = # col in temp.DF
		
	}
	
	DF.out # Return output dataframe
}


br.briere.RS <- briere_RateSummationCalc(tempprofs, briere.params)

br.briere.RS.long <- data.frame(Temp = c(rep(16, 7500), rep(20, 7500), rep(24, 7500), rep(28, 7500), rep(32, 7500), 
										rep(16, 7500), rep(20, 7500), rep(24, 7500), rep(28, 7500), rep(32, 7500)),
							   DTR = c(rep("+/- 9", 7500*5), rep("+/- 12", 7500*5)),
							   lifespan = c(br.briere.RS[,1], br.briere.RS[,2], br.briere.RS[,3], br.briere.RS[,4], br.briere.RS[,5],
							   			 br.briere.RS[,6], br.briere.RS[,7], br.briere.RS[,8], br.briere.RS[,9], br.briere.RS[,10]))

# Make summary df for plotting
br.briere.RS.summary <- br.briere.RS.long %>%
	group_by(DTR, Temp) %>%
	summarise(median = median(lifespan), lower = quantile(lifespan, 0.025), upper = quantile(lifespan, 0.975))

# Plot
p.br.RS <- ggplot() + 
	geom_point(data = alldata.sum, aes(x = Treatment, y = bite.rate.mean, color = DTR), size = 3) + xlim(12,38.6) + 
	scale_color_manual(values = c("dodgerblue", "mediumblue", "black")) + theme(legend.position = c(0.6,0.2)) +
	geom_line(data = br.c.briere.summary, aes(x = temp, y = mean)) +
	geom_line(data = br.c.briere.summary, aes(x = temp, y = upper), linetype = "dashed") +
	geom_line(data = br.c.briere.summary, aes(x = temp, y = lower), linetype = "dashed") +
	geom_point(data = br.briere.RS.summary, aes(x = Temp + 0.4, y = median, color = DTR), size = 3, shape = 15) +
	geom_errorbar(data = br.briere.RS.summary, aes(x = Temp + 0.4, ymin = lower, ymax = upper, color = DTR, width = 0)) +
	ylab("Biting rate (1/day)") + xlab("Temperature (°C)") 
p.br.RS

# Add in additional data observations to make full plot
br.modSS.RS.summary <- br.modSS.RS.summary %>%
	mutate(type = "RS.calc")
alldata.summary <- data.frame(DTR = alldata.sum.2$DTR,
							  Temp = alldata.sum.2$Treatment - 0.3,
							  median = alldata.sum.2$lifespan.mean,
							  upper = rep(NA,nrow(alldata.sum.2)),
							  lower = rep(NA,nrow(alldata.sum.2)),
							  type = "data")
plot.df <- rbind(alldata.summary, as.data.frame(lf.modSS.RS.summary))

# Plot
p.lf.RS.2 <- plot.df %>%
	ggplot() + 
	geom_point(aes(x = Temp, y = median, color = DTR, shape = type), size = 3) + ylim(7,50) + xlim(12,38.6) + 
	scale_color_manual(values = c("dodgerblue", "mediumblue", "black")) + theme(legend.position = c(0.6,0.8)) +
	geom_errorbar(aes(x = Temp, ymin = lower, ymax = upper, color = DTR, width = 0)) +
	geom_line(data = lf.c.SS.summary, aes(x = temp.C, y = mean)) +
	geom_line(data = lf.c.SS.summary, aes(x = temp.C, y = upper), linetype = "dashed") +
	geom_line(data = lf.c.SS.summary, aes(x = temp.C, y = lower), linetype = "dashed") +
	ylab("Lifespan (days)") + xlab("Temperature (°C)") 
p.lf.RS.2
