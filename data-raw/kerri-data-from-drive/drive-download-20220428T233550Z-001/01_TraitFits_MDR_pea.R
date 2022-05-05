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
data.pea.MDR <- read.csv("data/Krijn_Raw_Data.csv")
#replace Na with 0s
data.pea.MDR[is.na(data.pea.MDR)] = 0

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
#create a new column for individual EFD.1stGC*lifespan to get an individual-based estimated lifetime egg production
plot(Pea ~ temp, data = data.pea.MDR) # quad
plot(MDR ~ temp, data = data.pea.MDR) # briere


############## Trait name and parameter symbol - model type (e.g. quadratic, Briere)
# Get data
data.specific <- with(data.pea.MDR, data.frame('T' = temp, 'trait' = Pea)) # subset specific trait data from complete list
data <- data.specific # assign trait data to variable 'data'

# Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)

#############same data but now the truncated briere
# Run JAGS - **select correct model file**
pea.c.briere_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_T.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
pea.c.briere_T$BUGSoutput$summary
mcmcplot(pea.c.briere_T)
pea.c.briere_T.summary <- trait.trajs.briere(pea.c.briere_T, Temp.xs, summary = TRUE)
DIC.pea.briere.withT <- pea.c.briere_T$BUGSoutput$DIC

# Save model output 
save(pea.c.briere_T, file = "saved posteriors/pea_briere_withT_uniform.Rdata")

# Plot data and trait fit
p.pea.briere_T.fit.c <- data.pea.MDR %>%
  ggplot() + 
  geom_point(aes(x = temp, y = Pea), color = "black", size = 3) +
  geom_line(data = pea.c.briere_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = pea.c.briere_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = pea.c.briere_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Prob. e2a survival") + xlab("Temperature (°C)") 
p.pea.briere_T.fit.c 
ggsave("plots/pea_c.pea_briere_T.fit_only.pdf")


#############same data but now the quadratic model
# Run JAGS - **select correct model file**
pea.c.quad_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_T.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
pea.c.quad_T$BUGSoutput$summary
mcmcplot(pea.c.quad_T)
pea.c.quad_T.summary <- trait.trajs.quad(pea.c.quad_T, Temp.xs, summary = TRUE)
DIC.pea.quad.withT <- pea.c.quad_T$BUGSoutput$DIC

# Save model output 
save(pea.c.quad_T, file = "saved posteriors/pea_quad_withT_uniform.Rdata")

# Plot data and trait fit
p.pea.quad_T.fit.c <- data.pea.MDR %>%
  ggplot() + 
  geom_point(aes(x = temp, y = Pea), color = "black", size = 3) +
  geom_line(data = pea.c.quad_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = pea.c.quad_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = pea.c.quad_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Prob. e2a survival") + xlab("Temperature (°C)") 

p.pea.quad_T.fit.c 
ggsave("plots/pea_c.pea_quad_T.fit_only.pdf")


###########same data but now the quadratic function can become negative in values
pea.c.quad_neg <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_neg.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
pea.c.quad_neg$BUGSoutput$summary
mcmcplot(pea.c.quad_neg)
pea.c.quad_neg.summary <- trait.trajs.neg.quad(pea.c.quad_neg, Temp.xs, summary = TRUE)
DIC.pea.quad.neg <- pea.c.quad_neg$BUGSoutput$DIC

# Save model output 
save(pea.c.quad_neg, file = "saved posteriors/pea_neg_quad_uniform.Rdata")

# Plot data and trait fit
p.pea.quad.neg.fit.c <- data.pea.MDR %>%
  ggplot() + 
  geom_point(aes(x = temp, y = Pea), color = "black", size = 3) +
  geom_line(data = pea.c.quad_neg.summary, aes(x = temp, y = mean)) +
  geom_line(data = pea.c.quad_neg.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = pea.c.quad_neg.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Prob. e2a survival") + xlab("Temperature (°C)") 

p.pea.quad.neg.fit.c 
ggsave("plots/pea_c.neg_pea_quad.fit_only.pdf")

plot(pea.c.quad_neg.summary$mean~pea.c.quad_neg.summary$temp, ylim=c(0,1),xlab = "Temp")



#############################################################################
#############
#############################################################################
######### Mosquito developement rate
###############
#####################################################################
############## Trait name and parameter symbol - model type (e.g. quadratic, Briere)
# Get data
data.specific <- with(data.pea.MDR, data.frame('T' = temp, 'trait' = MDR)) # subset specific trait data from complete list
data <- data.specific # assign trait data to variable 'data'

# Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)

#############same data but now the truncated briere
# Run JAGS - **select correct model file**
MDR.c.briere_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_T.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
MDR.c.briere_T$BUGSoutput$summary
mcmcplot(MDR.c.briere_T)
MDR.c.briere_T.summary <- trait.trajs.briere(MDR.c.briere_T, Temp.xs, summary = TRUE)
DIC.MDR.briere.withT <- MDR.c.briere_T$BUGSoutput$DIC

# Save model output 
save(MDR.c.briere_T, file = "saved posteriors/MDR_briere_withT_uniform.Rdata")

# Plot data and trait fit
p.MDR.briere_T.fit.c <- data.pea.MDR %>%
  ggplot() + 
  geom_point(aes(x = temp, y = MDR), color = "black", size = 3) +
  geom_line(data = MDR.c.briere_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = MDR.c.briere_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = MDR.c.briere_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Developement rate (days-1)") + xlab("Temperature (°C)") 
p.MDR.briere_T.fit.c 
ggsave("plots/MDR_c_briere_T.fit_only.pdf")


#############same data but now the quadratic model
# Run JAGS - **select correct model file**
MDR.c.quad_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_T.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
MDR.c.quad_T$BUGSoutput$summary
mcmcplot(MDR.c.quad_T)
MDR.c.quad_T.summary <- trait.trajs.quad(MDR.c.quad_T, Temp.xs, summary = TRUE)
DIC.MDR.quad.withT <- MDR.c.quad_T$BUGSoutput$DIC

# Save model output 
save(MDR.c.quad_T, file = "saved posteriors/MDR_quad_withT_uniform.Rdata")

# Plot data and trait fit
p.MDR.quad_T.fit.c <- data.pea.MDR %>%
  ggplot() + 
  geom_point(aes(x = temp, y = MDR), color = "black", size = 3) +
  geom_line(data = MDR.c.quad_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = MDR.c.quad_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = MDR.c.quad_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Development rate (days-1)") + xlab("Temperature (°C)") 

p.MDR.quad_T.fit.c 
ggsave("plots/MDR_c_quad_T.fit_only.pdf")


###Same data but allowing for non-truncation
# Run JAGS - **select correct model file**
MDR.c.briere_neg <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_neg.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine model output & run diagnostics
MDR.c.briere_neg$BUGSoutput$summary
mcmcplot(MDR.c.briere_neg)
MDR.c.briere_neg.summary <- trait.trajs.neg.briere(MDR.c.briere_neg, Temp.xs, summary = TRUE)
DIC.MDR.c.briere_neg <- MDR.c.briere_neg$BUGSoutput$DIC

# Save model output 
save(MDR.c.briere_neg, file = "saved posteriors/MDR_briere_neg_uniform.Rdata")

# Plot data and trait fit
p.MDR.briere_neg_T.fit.c <- data.pea.MDR %>%
  ggplot() + 
  geom_point(aes(x = temp, y = MDR), color = "black", size = 3) +
  geom_line(data = MDR.c.briere_neg.summary, aes(x = temp, y = mean)) +
  geom_line(data = MDR.c.briere_neg.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = MDR.c.briere_neg.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Development rate (days-1)") + xlab("Temperature (°C)") 

p.MDR.briere_neg_T.fit.c 
ggsave("plots/MDR_c_briere_neg_T.fit_only.pdf")

