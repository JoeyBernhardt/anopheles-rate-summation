###Code written by Marta Shocket November 2019
#Later modified by Kerri Miazgowicz November 2019 for additional trait fits
##Later modified by KM April 24, 2020 to fit over the data means

################################# Contents
# 1. Set-up + data vis 
# 2. Bayesian Settings (code to write model files is in 01_TraitFit_Functions.R)
# 3. Bayesian Fitting - gamma

################################# 1. Set-up + data vis

###### Load packages
library(dplyr)
library(tidyverse)
library(R2jags)
library(coda)
library(cowplot)
library(mcmcplots)

#Set working directory
mainDir = "C:/Users/Kerri/Desktop/Chapter2 InProgress"
setwd(mainDir)

######  Load raw trait data
c.gamma <- data.frame(read_csv("data/constant.gamma.values.csv"))
f9.gamma <- data.frame(read_csv("data/dtr9.gamma.values.csv"))
f12.gamma <- data.frame(read_csv("data/dtr12.gamma.values.csv"))

#change default plot specifications
par(mar=c(1,1,1,1))
dev.off()
###### Summarize + plot data for all 3 traits: lifespan, biting rate, fecundity
############## First plot each trait to check functional forms
plot(gamma ~ temp, data = c.gamma) 
plot(gamma ~temp, data = f9.gamma)
plot(gamma~temp, data = f12.gamma)

# Add in DTR Treatment columns, combine in a single dataframe
c.gamma$DTR <- "constant"
f9.gamma$DTR <- "+/- 4.5"
f12.gamma$DTR <- "+/- 6"
all.data.block.temp <- rbind(c.gamma, f9.gamma, f12.gamma)

# Plot gamma
p.gamma <- ggplot() + 
	geom_point(data = all.data.block.temp, aes(x = jitter(temp,0.3), y = gamma, color = DTR), size =2) +
  scale_color_manual(values = c("dodgerblue", "blue", "black")) + theme(legend.position = "none") +
	ylab("Prop. surviving latancy period") + xlab("Temperature (°C)") 
p.gamma
par(mfrow = c(3,3))
plot_grid(p.gamma,NULL,NULL, nrow = 1)
ggsave("plots/gamma.data.pdf",   width = 8, height = 2.5, units = "in" ) #Saving 11.7 x 4.24 in image
dev.off()

################################# 2. Bayesian Settings

### MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 25000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains

### Temperature sequence for derived quantity calculations
#Temp.xs <- seq(5, 45, 1)  #Marta's original settings had at 1C increments
load("saved posteriors/temps.Rdata")
N.Temp.xs <-length(Temp.xs)



###################Load the set of parameters depending on the model to be run

### Quadratic/Briere function parameters and inits
# Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma")

# Initial values for the parameters
inits<-function(){list(
	cf.q = 0.01,
	cf.Tm = 35,
	cf.T0 = 5,
	cf.sigma = rlnorm(1))}


################################# 3. Bayesian Fitting - Gamma
#a - constant
#b - dtr9
#c - dtr12

#############################################3
######3a. Gamma - constant temperature
##############################################
# Pull out data - *** Make sure Temperature units are correct (C or K as needed) ***

data <- all.data.block.temp %>%
  dplyr::filter(DTR == "constant") %>%
  rename(trait = gamma,
         T = temp)

# Save as vectors for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)

# Fit the model
 gamma.c.quad_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
 				  model.file="quad_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
 				  n.iter=ni, DIC=T, working.directory=getwd())	

 gamma.c.briere_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                       model.file="briere_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                       n.iter=ni, DIC=T, working.directory=getwd())	
 
mcmcplot(gamma.c.quad_T)
mcmcplot(gamma.c.briere_T)

gamma.c.quad_T$BUGSoutput$summary
gamma.c.briere_T$BUGSoutput$summary

gamma.c.quad_T$BUGSoutput$DIC
gamma.c.briere_T$BUGSoutput$DIC

# Calculate summary statistics of trait trajectories for plotting
gamma.c.quad_T.summary <- trait.trajs.quad(gamma.c.quad_T, Temp.xs, summary = TRUE)
gamma.c.briere_T.summary <- trait.trajs.briere(gamma.c.briere_T, Temp.xs, summary = TRUE)

# Save model output 
save(gamma.c.quad_T, file = "saved posteriors/constant_gamma_c.quad_T.uniform.Rdata")
save(gamma.c.briere_T, file = "saved posteriors/constant_gamma.c.briere_T.uniform.Rdata")

# Plot data and trait fit
p.gamma.quad_T.fit.c <- all.data.block.temp %>%
	filter(DTR == "constant") %>%
	ggplot() + 
	geom_point(aes(x = temp, y = gamma), color = "black", size = 3) +
  	geom_line(data = gamma.c.quad_T.summary, aes(x = temp, y = mean)) +
  	geom_line(data = gamma.c.quad_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  	geom_line(data = gamma.c.quad_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  	ylab("Gamma") + xlab("Temperature (°C)") 
p.gamma.quad_T.fit.c 
ggsave("plots/gamma_c.lf_quad_T.fit_only.pdf")

# Plot data and trait fit
p.gamma.briere_T.fit.c <- all.data.block.temp %>%
  filter(DTR == "constant") %>%
  ggplot() + 
  geom_point(aes(x = temp, y = gamma), color = "black", size = 3) +
  geom_line(data = gamma.c.briere_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = gamma.c.briere_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = gamma.c.briere_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Gamma") + xlab("Temperature (°C)") 
p.gamma.briere_T.fit.c 
ggsave("plots/gamma_c.lf_briere_T.fit_only.pdf")

######################################
##########3b. lifespan dtr 9 models
#####################################

# Pull out data - *** Make sure Temperature units are correct (C or K as needed) ***

data <- all.data.block.temp %>%
  dplyr::filter(DTR == "+/- 4.5") %>%
  rename(trait = gamma,
         T = temp)

# Save as vectors for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)

# Fit the model
gamma.f9.quad_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                    model.file="quad_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                    n.iter=ni, DIC=T, working.directory=getwd())	

gamma.f9.briere_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                      model.file="briere_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                      n.iter=ni, DIC=T, working.directory=getwd())	

mcmcplot(gamma.f9.quad_T)
mcmcplot(gamma.f9.briere_T)

gamma.f9.quad_T$BUGSoutput$summary
gamma.f9.briere_T$BUGSoutput$summary

gamma.f9.quad_T$BUGSoutput$DIC
gamma.f9.briere_T$BUGSoutput$DIC

# Calculate summary statistics of trait trajectories for plotting
gamma.f9.quad_T.summary <- trait.trajs.quad(gamma.f9.quad_T, Temp.xs, summary = TRUE)
gamma.f9.briere_T.summary <- trait.trajs.briere(gamma.f9.briere_T, Temp.xs, summary = TRUE)

# Save model output 
save(gamma.f9.quad_T, file = "saved posteriors/dtr_gamma_f9.quad_T.uniform.Rdata")
save(gamma.f9.briere_T, file = "saved posteriors/dtr_gamma.f9.briere_T.uniform.Rdata")

# Plot data and trait fit
p.gamma.quad_T.fit.f9 <- all.data.block.temp %>%
  filter(DTR == "+/- 4.5") %>%
  ggplot() + 
  geom_point(aes(x = temp, y = gamma), color = "dodgerblue", size = 3) +
  geom_line(data = gamma.f9.quad_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = gamma.f9.quad_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = gamma.f9.quad_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Gamma") + xlab("Temperature (°C)") 
p.gamma.quad_T.fit.f9 
ggsave("plots/gamma_f9_quad_T.fit_only.pdf")

# Plot data and trait fit
p.gamma.briere_T.fit.f9 <- all.data.block.temp %>%
  filter(DTR == "+/- 4.5") %>%
  ggplot() + 
  geom_point(aes(x = temp, y = gamma), color = "dodgerblue", size = 3) +
  geom_line(data = gamma.f9.briere_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = gamma.f9.briere_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = gamma.f9.briere_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Gamma") + xlab("Temperature (°C)") 
p.gamma.briere_T.fit.f9 
ggsave("plots/gamma_f9_briere_T.fit_only.pdf")

##############################################################
#################3c - Lifespan - dtr 12
##########################################################

# Pull out data - *** Make sure Temperature units are correct (C or K as needed) ***

data <- all.data.block.temp %>%
  dplyr::filter(DTR == "+/- 6") %>%
  rename(trait = gamma,
         T = temp)

# Save as vectors for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)

# Fit the model
gamma.f12.quad_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                     model.file="quad_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                     n.iter=ni, DIC=T, working.directory=getwd())	

gamma.f12.briere_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                       model.file="briere_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                       n.iter=ni, DIC=T, working.directory=getwd())	

mcmcplot(gamma.f12.quad_T)
mcmcplot(gamma.f12.briere_T)

gamma.f12.quad_T$BUGSoutput$summary
gamma.f12.briere_T$BUGSoutput$summary

gamma.f12.quad_T$BUGSoutput$DIC
gamma.f12.briere_T$BUGSoutput$DIC

# Calculate summary statistics of trait trajectories for plotting
gamma.f12.quad_T.summary <- trait.trajs.quad(gamma.f12.quad_T, Temp.xs, summary = TRUE)
gamma.f12.briere_T.summary <- trait.trajs.briere(gamma.f12.briere_T, Temp.xs, summary = TRUE)

# Save model output 
save(gamma.f12.quad_T, file = "saved posteriors/dtr_gamma_f12.quad_T.uniform.Rdata")
save(gamma.f12.briere_T, file = "saved posteriors/dtr_gamma_f12.briere_T.uniform.Rdata")

# Plot data and trait fit
p.gamma.quad_T.fit.f12 <- all.data.block.temp %>%
  filter(DTR == "+/- 6") %>%
  ggplot() + 
  geom_point(aes(x = temp, y = gamma), color = "royalblue", size = 3) +
  geom_line(data = gamma.f12.quad_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = gamma.f12.quad_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = gamma.f12.quad_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Gamma") + xlab("Temperature (°C)") 
p.gamma.quad_T.fit.f12
ggsave("plots/gamma_f12_quad_T.fit_only.pdf")

# Plot data and trait fit
p.gamma.briere_T.fit.f12 <- all.data.block.temp %>%
  filter(DTR == "+/- 6") %>%
  ggplot() + 
  geom_point(aes(x = temp, y = gamma), color = "royalblue", size = 3) +
  geom_line(data = gamma.f12.briere_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = gamma.f12.briere_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = gamma.f12.briere_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Gamma") + xlab("Temperature (°C)") 
p.gamma.briere_T.fit.f12
ggsave("plots/gamma_f12_briere_T.fit_only.pdf")
