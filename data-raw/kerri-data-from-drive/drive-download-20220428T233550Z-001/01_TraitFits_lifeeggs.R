###Code written by Marta Shocket November 2019
#Later modified by Kerri Miazgowicz November 2019 for additional trait fits
##Later modified by KM April 23, 2020 to fit over the data means

################################# Contents
# 1. Set-up + data vis 
# 2. Bayesian Settings (code to write model files is in 01_TraitFit_Functions.R)
# 3. Bayesian Fitting - lifespan
# 4. Bayesian Fitting - biting rate
# 5. Bayesian Fitting - fecundity

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
cdata <- data.frame(read_csv("data/constant.individual.trait.csv"))
f9data <- data.frame(read_csv("data/fluc9.individual.trait.csv"))
f12data <- data.frame(read_csv("data/fluc12.individual.trait.csv"))

#change default plot specifications
par(mar=c(1,1,1,1))
dev.off()
###### Summarize + plot data for all 3 traits: lifespan, biting rate, fecundity
############## First plot each trait to check functional forms
boxplot(lifetime.eggs~Treatment, data = cdata)
boxplot(lifetime.eggs~Treatment, data = f9data)
boxplot(lifetime.eggs~Treatment, data = f12data)

#######Generate Block specific data values for fitting the curve over 
c.data.block.temp <- cdata %>%
  dplyr::group_by(Treatment,Block)%>%
  dplyr::summarise(lifetime.eggs = mean(lifetime.eggs)) %>%
  ungroup()

plot(lifetime.eggs ~ Treatment, data = c.data.block.temp) 

################dtr 9#################
f9.data.block.temp <- f9data %>%
  dplyr::group_by(Treatment,Block)%>%
  dplyr::summarise(lifetime.eggs = mean(lifetime.eggs)) %>%
  ungroup()

plot(lifetime.eggs ~ Treatment, data = f9.data.block.temp) 
###############dtr 12##################3
f12.data.block.temp <- f12data %>%
  dplyr::group_by(Treatment,Block)%>%
  dplyr::summarise(lifetime.eggs = mean(lifetime.eggs)) %>%
  ungroup()

plot(lifetime.eggs ~ Treatment, data = f12.data.block.temp) 

# Add in DTR Treatment columns, combine in a single dataframe
c.data.block.temp$DTR <- "constant"
f9.data.block.temp$DTR <- "+/- 4.5"
f12.data.block.temp$DTR <- "+/- 6"
all.data.block.temp <- rbind(c.data.block.temp, f9.data.block.temp, f12.data.block.temp)

# Plot fecundity
p.f <- ggplot() + 
	geom_point(data = all.data.block.temp, aes(x = jitter(Treatment,0.3), y = lifetime.eggs, color = DTR), size =2) +
  scale_color_manual(values = c("dodgerblue", "blue", "black")) + theme(legend.position = "none") +
	ylab("Lifetime fecundity (eggs)") + xlab("Temperature (°C)") 

p.f
################################# 2. Bayesian Settings

### MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 25000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains

### Temperature sequence for derived quantity calculations
#Temp.xs <- seq(5, 45, 1)  #Marta's original settings had at 1C increments
#load("saved posteriors/temps.Rdata")

Temp.xs <- seq(0,50,.1) #For better mapping needs later changed the increment to 0.1
N.Temp.xs <-length(Temp.xs)
Temp.xs.K <- Temp.xs + 273.15 # Needed for Sharpe-Schoolfied type models

# save the temperature sequence for future analyses
save(Temp.xs, file = "saved posteriors/temps.Rdata")#named identically as previous temps.Rdata

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


################################# 3. Bayesian Fitting - lifeeggs
#a - constant
#b - dtr9
#c - dtr12

#############################################3
######3a. Lifetime egg production - constant temperature
##############################################
# Pull out data - *** Make sure Temperature units are correct (C or K as needed) ***

data <- c.data.block.temp %>%
  dplyr::filter(DTR == "constant") %>%
  rename(trait = lifetime.eggs,
         T = Treatment)

# Save as vectors for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)

# Fit the model
 le.c.quad_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
 				  model.file="quad_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
 				  n.iter=ni, DIC=T, working.directory=getwd())	

 le.c.briere_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                       model.file="briere_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                       n.iter=ni, DIC=T, working.directory=getwd())	
 
mcmcplot(le.c.quad_T)
mcmcplot(le.c.briere_T)

le.c.quad_T$BUGSoutput$summary
le.c.briere_T$BUGSoutput$summary

le.c.quad_T$BUGSoutput$DIC
le.c.briere_T$BUGSoutput$DIC

# Calculate summary statistics of trait trajectories for plotting
le.c.quad_T.summary <- trait.trajs.quad(le.c.quad_T, Temp.xs, summary = TRUE)
le.c.briere_T.summary <- trait.trajs.briere(le.c.briere_T, Temp.xs, summary = TRUE)

# Save model output 
save(le.c.quad_T, file = "saved posteriors/constant_le_c.quad_T.uniform.Rdata")
save(le.c.briere_T, file = "saved posteriors/constant_le.c.briere_T.uniform.Rdata")

# Plot data and trait fit
p.le.quad_T.fit.c <- all.data.block.temp %>%
	filter(DTR == "constant") %>%
	ggplot() + 
	geom_point(aes(x = Treatment, y = lifetime.eggs), color = "black", size = 3) +
  	geom_line(data = le.c.quad_T.summary, aes(x = temp, y = mean)) +
  	geom_line(data = le.c.quad_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  	geom_line(data = le.c.quad_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  	ylab("Lifetime fecundity (eggs)") + xlab("Temperature (°C)") 
p.le.quad_T.fit.c 
ggsave("plots/lifeeggs_c.le_quad_T.fit_only.pdf")

# Plot data and trait fit
p.le.briere_T.fit.c <- all.data.block.temp %>%
  filter(DTR == "constant") %>%
  ggplot() + 
  geom_point(aes(x = Treatment, y = lifetime.eggs), color = "black", size = 3) +
  geom_line(data = le.c.briere_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = le.c.briere_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = le.c.briere_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Lifetime fecundity (eggs)") + xlab("Temperature (°C)") 
p.le.briere_T.fit.c 
ggsave("plots/lifeeggs_c.le_briere_T.fit_only.pdf")

######################################
##########3b. lifeeggs dtr 9 models
#####################################

# Pull out data - *** Make sure Temperature units are correct (C or K as needed) ***

data <- f9.data.block.temp %>%
  dplyr::filter(DTR == "+/- 4.5") %>%
  rename(trait = lifetime.eggs,
         T = Treatment)

# Save as vectors for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)

# Fit the model
le.f9.quad_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                    model.file="quad_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                    n.iter=ni, DIC=T, working.directory=getwd())	

le.f9.briere_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                      model.file="briere_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                      n.iter=ni, DIC=T, working.directory=getwd())	

mcmcplot(le.f9.quad_T)
mcmcplot(le.f9.briere_T)

le.f9.quad_T$BUGSoutput$summary
le.f9.briere_T$BUGSoutput$summary

le.f9.quad_T$BUGSoutput$DIC
le.f9.briere_T$BUGSoutput$DIC

# Calculate summary statistics of trait trajectories for plotting
le.f9.quad_T.summary <- trait.trajs.quad(le.f9.quad_T, Temp.xs, summary = TRUE)
le.f9.briere_T.summary <- trait.trajs.briere(le.f9.briere_T, Temp.xs, summary = TRUE)

# Save model output 
save(le.f9.quad_T, file = "saved posteriors/dtr_le_f9.quad_T.uniform.Rdata")
save(le.f9.briere_T, file = "saved posteriors/dtr_le.f9.briere_T.uniform.Rdata")

# Plot data and trait fit
p.le.quad_T.fit.f9 <- all.data.block.temp %>%
  filter(DTR == "+/- 4.5") %>%
  ggplot() + 
  geom_point(aes(x = Treatment, y = lifetime.eggs), color = "dodgerblue", size = 3) +
  geom_line(data = le.f9.quad_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = le.f9.quad_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = le.f9.quad_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Lifetime fecundity (eggs)") + xlab("Temperature (°C)") 
p.le.quad_T.fit.f9 
ggsave("plots/lifeeggs_f9.le_quad_T.fit_only.pdf")

# Plot data and trait fit
p.le.briere_T.fit.f9 <- all.data.block.temp %>%
  filter(DTR == "+/- 4.5") %>%
  ggplot() + 
  geom_point(aes(x = Treatment, y = lifetime.eggs), color = "dodgerblue", size = 3) +
  geom_line(data = le.f9.briere_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = le.f9.briere_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = le.f9.briere_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Lifetime fecundity (eggs)") + xlab("Temperature (°C)") 
p.le.briere_T.fit.f9 
ggsave("plots/lifeeggs_f9.le_briere_T.fit_only.pdf")

##############################################################
#################3c - Lifetime eggs - dtr 12
##########################################################

# Pull out data - *** Make sure Temperature units are correct (C or K as needed) ***

data <- f12.data.block.temp %>%
  dplyr::filter(DTR == "+/- 6") %>%
  rename(trait = lifetime.eggs,
         T = Treatment)

# Save as vectors for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)

# Fit the model
le.f12.quad_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                     model.file="quad_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                     n.iter=ni, DIC=T, working.directory=getwd())	


le.f12.briere_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                       model.file="briere_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                       n.iter=ni, DIC=T, working.directory=getwd())	

mcmcplot(le.f12.quad_T)
mcmcplot(le.f12.briere_T)

le.f12.quad_T$BUGSoutput$summary
le.f12.briere_T$BUGSoutput$summary

le.f12.quad_T$BUGSoutput$DIC
le.f12.briere_T$BUGSoutput$DIC

# Calculate summary statistics of trait trajectories for plotting
le.f12.quad_T.summary <- trait.trajs.quad(le.f12.quad_T, Temp.xs, summary = TRUE)
le.f12.briere_T.summary <- trait.trajs.briere(le.f12.briere_T, Temp.xs, summary = TRUE)

# Save model output 
save(le.f12.quad_T, file = "saved posteriors/dtr_le_f12.quad_T.uniform.Rdata")
save(le.f12.briere_T, file = "saved posteriors/dtr_le.f12.briere_T.uniform.Rdata")

# Plot data and trait fit
p.le.quad_T.fit.f12 <- all.data.block.temp %>%
  filter(DTR == "+/- 6") %>%
  ggplot() + 
  geom_point(aes(x = Treatment, y = lifetime.eggs), color = "royalblue", size = 3) +
  geom_line(data = le.f12.quad_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = le.f12.quad_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = le.f12.quad_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Lifetime fecundity (eggs)") + xlab("Temperature (°C)") 
p.le.quad_T.fit.f12
ggsave("plots/lifeeggs_f12.le_quad_T.fit_only.pdf")

# Plot data and trait fit
p.le.briere_T.fit.f12 <- all.data.block.temp %>%
  filter(DTR == "+/- 6") %>%
  ggplot() + 
  geom_point(aes(x = Treatment, y = lifetime.eggs), color = "royalblue", size = 3) +
  geom_line(data = le.f12.briere_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = le.f12.briere_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = le.f12.briere_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Lifetime fecundity (eggs)") + xlab("Temperature (°C)") 
p.le.briere_T.fit.f12
ggsave("plots/lifeeggs_f12.le_briere_T.fit_only.pdf")

