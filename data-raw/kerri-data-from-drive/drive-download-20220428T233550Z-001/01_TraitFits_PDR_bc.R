## Marta Shocket, Stanford University
## Updated by Erin Mordecai, August 6, 2018
## Updated by Kerri Miazgowicz between August 6, 2018 and August 5, 2019
### Modified by KM on April 12, 2020 to compare model fits using
## KM Aprial 23,2020
##            i. the full individual dataset vs. block|temp means
##            ii. normal distribution vs. truncated normal distributions
## Purpose: Use Bayesian Inference (JAGS) to fit temperature-dependent functions traits
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) JAGS models
##           3) Shared settings for all models
##           4) Code for fitting a trait - uniform priors
##           5) Example code for fitting a trait - informative priors fit from data
##           6) Example code for fitting prior distribution



##########
###### 1. Set up workspace, load packages, get data, etc.
##########

mainDir = "C:/Users/Kerri/Desktop/Chapter2 InProgress"
setwd(mainDir)

# Check whether there's a folder in the directory for saving plots
# If not, create one
subDir = "saved posteriors"
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)

# Load libraties for fitting traits
library(dplyr)
library(tidyverse)
library(R2jags)
library(coda)
library(cowplot)
library(mcmcplots)

#change default plot specifications
par(mar=c(1,1,1,1))

# Load raw trait data - for the code below to work, trait value should be called 'trait' & temperature should be called 'T'
data.bc.PDR <- read.csv("data/forErin_ShapiroData.csv")
#replace Na with 0s
data.bc.PDR[is.na(data.bc.PDR)] = 0

##########
###### 3. Shared settings for all models
##########

##### inits Function
### MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 25000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains

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


# save the temperature sequence for future analyses
load("saved posteriors/temps.Rdata")
##### Derived Quantity Settings
N.Temp.xs <-length(Temp.xs)


##########
###### 4. Fitting traits - uniform priors
##########

############## First plot each trait to check functional forms
plot(bc ~ temp, data = data.bc.PDR) 
plot(1/EIP50 ~ temp, data = data.bc.PDR) #now as PDR; fit as PDR

#####################################################
##########1. Pathogen developement rate
############## Trait name and parameter symbol - model type (e.g. quadratic, Briere)
# Get data
data.specific <- with(data.bc.PDR, data.frame('T' = temp, 'trait' = 1/EIP50)) # subset specific trait data from complete list
data <- data.specific # assign trait data to variable 'data'

# Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)

#############same data but now the truncated briere
# Run JAGS - **select correct model file**
PDR.c.briere_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_T.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
PDR.c.briere_T$BUGSoutput$summary
mcmcplot(PDR.c.briere_T)
PDR.c.briere_T.summary <- trait.trajs.briere(PDR.c.briere_T, Temp.xs, summary = TRUE)
DIC.PDR.briere.withT <- PDR.c.briere_T$BUGSoutput$DIC

# Save model output 
save(PDR.c.briere_T, file = "saved posteriors/PDR_briere_withT_uniform.Rdata")

# Plot data and trait fit
p.PDR.briere_T.fit.c <- data.bc.PDR %>%
  ggplot() + 
  geom_point(aes(x = temp, y = 1/EIP50), color = "black", size = 3) +
  geom_line(data = PDR.c.briere_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = PDR.c.briere_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = PDR.c.briere_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Development rate (days-1)") + xlab("Temperature (°C)") 
p.PDR.briere_T.fit.c 
ggsave("plots/PDR_c.pea_briere_T.fit_only.pdf")


#############same data but now the quadratic model
# Run JAGS - **select correct model file**
PDR.c.quad_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_T.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
PDR.c.quad_T$BUGSoutput$summary
mcmcplot(PDR.c.quad_T)
PDR.c.quad_T.summary <- trait.trajs.quad(PDR.c.quad_T, Temp.xs, summary = TRUE)
DIC.PDR.quad.withT <- PDR.c.quad_T$BUGSoutput$DIC

# Save model output 
save(PDR.c.quad_T, file = "saved posteriors/PDR_quad_withT_uniform.Rdata")

# Plot data and trait fit
p.PDR.quad_T.fit.c <- data.bc.PDR %>%
  ggplot() + 
  geom_point(aes(x = temp, y = 1/EIP50), color = "black", size = 3) +
  geom_line(data = PDR.c.quad_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = PDR.c.quad_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = PDR.c.quad_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Development rate (days-1)") + xlab("Temperature (°C)") 

p.PDR.quad_T.fit.c 
ggsave("plots/PDR_c.pea_quad_T.fit_only.pdf")


#############################################################################
#############
#############################################################################
######### Vector competence
###############
#####################################################################
############## Trait name and parameter symbol - model type (e.g. quadratic, Briere)
# Get data
data.specific <- with(data.bc.PDR, data.frame('T' = temp, 'trait' = bc)) # subset specific trait data from complete list
data <- data.specific # assign trait data to variable 'data'

# Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)

#############same data but now the truncated briere
# Run JAGS - **select correct model file**
bc.c.briere_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_T.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
bc.c.briere_T$BUGSoutput$summary
mcmcplot(bc.c.briere_T)
bc.c.briere_T.summary <- trait.trajs.briere(bc.c.briere_T, Temp.xs, summary = TRUE)
DIC.bc.briere.withT <- bc.c.briere_T$BUGSoutput$DIC

# Save model output 
save(bc.c.briere_T, file = "saved posteriors/bc_briere_withT_uniform.Rdata")

# Plot data and trait fit
p.bc.briere_T.fit.c <- data.bc.PDR %>%
  ggplot() + 
  geom_point(aes(x = temp, y = bc), color = "black", size = 3) +
  geom_line(data = bc.c.briere_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = bc.c.briere_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = bc.c.briere_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Vector competence") + xlab("Temperature (°C)") 
p.bc.briere_T.fit.c 
ggsave("plots/bc_c_briere_T.fit_only.pdf")


#############same data but now the quadratic model
# Run JAGS - **select correct model file**
bc.c.quad_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_T.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
bc.c.quad_T$BUGSoutput$summary
mcmcplot(bc.c.quad_T)
bc.c.quad_T.summary <- trait.trajs.quad(bc.c.quad_T, Temp.xs, summary = TRUE)
DIC.bc.quad.withT <- bc.c.quad_T$BUGSoutput$DIC

# Save model output 
save(bc.c.quad_T, file = "saved posteriors/bc_quad_withT_uniform.Rdata")

# Plot data and trait fit
p.bc.quad_T.fit.c <- data.bc.PDR %>%
  ggplot() + 
  geom_point(aes(x = temp, y = bc), color = "black", size = 3) +
  geom_line(data = bc.c.quad_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = bc.c.quad_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = bc.c.quad_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Vector competence") + xlab("Temperature (°C)") 

p.bc.quad_T.fit.c 
ggsave("plots/bc_c_quad_T.fit_only.pdf")

