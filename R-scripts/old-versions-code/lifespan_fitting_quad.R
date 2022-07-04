

#### 

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
model_out_constant_lifespan_quad <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
										   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Bundle Data
# jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
# data_constant_lifespan
# model_out_constant_lifespan_quad <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
#                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
# model.out$BUGSoutput$summary[1:10,]
# mcmcplot(model.out)

# Save model output 
save(model_out_constant_lifespan_quad, file = "saved-posteriors/constant_lifespan_quad_uniform.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_constant_lifespan,
	 ylab = "Lifespan", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_out_constant_lifespan_quad$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_out_constant_lifespan_quad$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_out_constant_lifespan_quad$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


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
model_out_lifespan_dtr9_quad <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
									   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model_out_lifespan_dtr9_quad$BUGSoutput$summary[1:10,]
mcmcplot(model_out_lifespan_dtr9_quad)

# Save model output 
save(model_out_lifespan_dtr9_quad, file = "saved-posteriors/dtr9_lifespan_quad_uniform.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_dtr9_lifespan, 
	 ylab = "Lifespan", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_out_lifespan_dtr9_quad$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_out_lifespan_dtr9_quad$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_out_lifespan_dtr9_quad$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



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
model_out_lifespan_dtr12 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
								 n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model_out_lifespan_dtr12$BUGSoutput$summary[1:10,]
mcmcplot(model_out_lifespan_dtr12)

# Save model output 
save(model_out_lifespan_dtr12, file = "saved-posteriors/dtr12_lifespan_quad_uniform.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_dtr12_lifespan, 
	 ylab = "Lifespan", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_out_lifespan_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_out_lifespan_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_out_lifespan_dtr12$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


# gathering results -------------------------------------------------------

predictions_lifespan_constant <- as.data.frame(model_out_constant_lifespan_quad$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)   ### columns are temperatures, rows are iterations
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

write_csv(predictions_lifespan_constant_summary, "data-processed/predictions_lifespan_constant_summary_quad.csv")

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

b_params_lifespan_constant <- as.data.frame(model_out_constant_lifespan_quad$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>% 
	mutate(treatment = "lifespan_constant")

params_lifespan_constant_all <- bind_rows(b_params_lifespan_constant, topt_lifespan_constant)
# View(params_lifespan_constant_all)

write_csv(params_lifespan_constant_all, "data-processed/params_lifespan_constant_all_quad.csv")

### raw data to plot
data_lifespan_constant_sum <- data_constant_lifespan %>% 
	dplyr::rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = std.error(trait)) %>% 
	mutate(trait = "bite rate") %>% 
	mutate(treatment = "lifespan_constant")

write_csv(data_lifespan_constant_sum, "data-processed/data_lifespan_constant_sum_quad.csv")



# 6. Lifespan DTR9 --------------------------------------------------------
### ok we need to come back to these results because they are showing this two grouped thing
data_dtr9_lifespan
model_out_lifespan_dtr9_quad



predictions_lifespan_dtr9 <- as.data.frame(model_out_lifespan_dtr9_quad$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)   ### columns are temperatures, rows are iterations
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

write_csv(predictions_lifespan_dtr9_summary, "data-processed/predictions_lifespan_dtr9_summary_quad.csv")

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

b_params_lifespan_dtr9 <- as.data.frame(model_out_lifespan_dtr9_quad$BUGSoutput$summary[1:5,]) %>%
	rownames_to_column(var = "term") %>% 
	mutate(treatment = "lifespan_dtr9")

params_lifespan_dtr9_all <- bind_rows(b_params_lifespan_dtr9, topt_lifespan_dtr9)
# View(params_lifespan_dtr9_all)

write_csv(params_lifespan_dtr9_all, "data-processed/params_lifespan_dtr9_all_quad.csv")

### raw data to plot
data_lifespan_dtr9_sum <- data_dtr9_lifespan %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = std.error(trait)) %>% 
	mutate(trait = "bite rate") %>% 
	mutate(treatment = "lifespan_dtr9")

write_csv(data_lifespan_dtr9_sum, "data-processed/data_lifespan_dtr9_sum_quad.csv")



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

write_csv(predictions_lifespan_dtr12_summary, "data-processed/predictions_lifespan_dtr12_summary_quad.csv")

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

write_csv(params_lifespan_dtr12_all, "data-processed/params_lifespan_dtr12_all_quad.csv")

### raw data to plot
data_lifespan_dtr12_sum <- data_dtr12_lifespan %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = std.error(trait)) %>% 
	mutate(trait = "bite rate") %>% 
	mutate(treatment = "lifespan_dtr12")

write_csv(data_lifespan_dtr12_sum, "data-processed/data_lifespan_dtr12_sum_quad.csv")



# now eggs quad -----------------------------------------------------------

Model fitting: Lifetime eggs constant -----------------------------------
	
	
	############ Lifetime eggs at constant temperature
	# Get data
data_constant_eggs <- with(data.constant, data.frame('T' = temp, 'trait' = lifetime.eggs)) # subset specific trait data from complete list
# data <- data.specific # assign trait data to variable 'data'

# Organize Data for JAGS
trait <- data_constant_eggs$trait
N.obs <- length(trait)
temp <- data_constant_eggs$T

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS - **select correct model file**
# data_constant_eggs
model_out_eggs_constant_quad <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
								n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())


str(model_out_eggs_constant_quad) 
str(model_out_eggs_constant) 
# Examine model output & run diagnostics
model_out_eggs_constant$BUGSoutput$summary[1:10,]
mcmcplot(model_out_eggs_constant)

# Save model output 
save(model_out_eggs_constant, file = "saved-posteriors/constant_eggs_quad_uniform.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_constant_eggs, 
	 ylab = "Lifetime eggs", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_out_eggs_constant_quad$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_out_eggs_constant_quad$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_out_eggs_constant_quad$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

str(model_out_eggs_constant)
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
model_out_eggs_dtr9_quad <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
							n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
# model_out_eggs_dtr9$BUGSoutput$summary[1:10,]
# mcmcplot(model_out_eggs_dtr9)

# Save model output 
# save(model_out_eggs_dtr9, file = "saved-posteriors/dtr9_eggs_quad_uniform.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_eggs_dtr9, 
	 ylab = "Lifetime eggs", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_out_eggs_dtr9_quad$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_out_eggs_dtr9_quad$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_out_eggs_dtr9_quad$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


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
model_out_eggs_dtr12_quad <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
							 n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
model_out_eggs_dtr12$BUGSoutput$summary[1:10,]
mcmcplot(model_out_eggs_dtr12)

# Save model output 
save(model_out_eggs_dtr12, file = "saved-posteriors/dtr12_eggs_quad_uniform.Rdata")

# Plot trait data, model mean and CIs
plot(trait ~ jitter(T, 0.5), xlim = c(0, 45), data = data_eggs_dtr12, 
	 ylab = "Lifetime eggs", xlab = expression(paste("Temperature (",degree,"C)")))
lines(model_out_eggs_dtr12_quad$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(model_out_eggs_dtr12_quad$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(model_out_eggs_dtr12_quad$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


#### Eggs model outputs quad

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

write_csv(predictions_eggs_constant_summary, "data-processed/predictions_eggs_constant_summary_quad.csv")

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

write_csv(params_constant_eggs_all, "data-processed/params_constant_eggs_all_quad.csv")

### raw data to plot
data_constant_eggs_sum <- data_constant_eggs %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = std.error(trait)) %>% 
	mutate(trait = "eggs") %>% 
	mutate(treatment = "eggs_constant")

write_csv(data_constant_eggs_sum, "data-processed/data_constant_eggs_sum_quad.csv")




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

write_csv(predictions_eggs_dtr9_summary, "data-processed/predictions_eggs_dtr9_summary_quad.csv")

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

write_csv(params_dtr9_eggs_all, "data-processed/params_dtr9_eggs_all_quad.csv")

### raw data to plot
data_dtr9_eggs_sum <- data_eggs_dtr9 %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = std.error(trait)) %>% 
	mutate(trait = "eggs") %>% 
	mutate(treatment = "eggs_dtr9")

write_csv(data_dtr9_eggs_sum, "data-processed/data_dtr9_eggs_sum_quad.csv")


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

write_csv(predictions_eggs_dtr12_summary, "data-processed/predictions_eggs_dtr12_summary_quad.csv")

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

write_csv(params_dtr12_eggs_all, "data-processed/params_dtr12_eggs_all_quad.csv")

### raw data to plot
data_dtr12_eggs_sum <- data_eggs_dtr12 %>%
	rename("temperature" = "T") %>% 
	group_by(temperature) %>% 
	summarise(mean =  mean(trait),
			  std_error = std.error(trait)) %>% 
	mutate(trait = "eggs") %>% 
	mutate(treatment = "eggs_dtr12")

write_csv(data_dtr12_eggs_sum, "data-processed/data_dtr12_eggs_sum_quad.csv")



