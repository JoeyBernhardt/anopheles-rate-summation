###Code written by Marta Shocket November 2019
#Later modified by Kerri Miazgowicz November 2019 for additional trait fits
##Later modified by KM April 22, 2020 to fit over the data means

################################# Contents
# 1. Set-up + data vis 
# 2. Bayesian Settings (code to write model files is in 01_TraitFit_Functions.R)
# 3. Bayesian Fitting - lifespan
# 4. Bayesian Fitting - biting rate
# 5. Bayesian Fitting - fecundity

################################# 1. Set-up + data vis

###### Load packages
library(plyr)
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
boxplot(bite.rate~Treatment, data = cdata)
boxplot(bite.rate~Treatment, data = f9data)
boxplot(bite.rate~Treatment, data = f12data)

#######Generate Block specific data values for fitting the curve over 
c.data.block.temp <- cdata %>%
  dplyr::group_by(Treatment,Block)%>%
  dplyr::summarise(lifespan = mean(lifespan),
                   bite.rate = mean(bite.rate),
                   lifetime.eggs = mean(lifetime.eggs)) %>%
  ungroup()

plot(bite.rate ~ Treatment, data = c.data.block.temp) 

################dtr 9#################
f9.data.block.temp <- f9data %>%
  dplyr::group_by(Treatment,Block)%>%
  dplyr::summarise(lifespan = mean(lifespan),
                   bite.rate = mean(bite.rate),
                   lifetime.eggs = mean(lifetime.eggs)) %>%
  ungroup()

plot(bite.rate ~ Treatment, data = f9.data.block.temp) 
###############dtr 12##################3
f12.data.block.temp <- f12data %>%
  dplyr::group_by(Treatment,Block)%>%
  dplyr::summarise(lifespan = mean(lifespan),
                   bite.rate = mean(bite.rate),
                   lifetime.eggs = mean(lifetime.eggs)) %>%
  ungroup()

plot(bite.rate ~ Treatment, data = f12.data.block.temp) 




# Add in DTR Treatment columns, combine in a single dataframe
c.data.block.temp$DTR <- "constant"
f9.data.block.temp$DTR <- "+/- 4.5"
f12.data.block.temp$DTR <- "+/- 6"
all.data.block.temp <- rbind(c.data.block.temp, f9.data.block.temp, f12.data.block.temp)

all.data.block.temp <- all.data.block.temp %>% 
           dplyr::mutate(TempK = Treatment + 273.15)

# Plot biting rate
p.br <- ggplot() + 
	geom_point(data = all.data.block.temp, aes(x = jitter(Treatment,0.3), y = bite.rate, color = DTR), size =2) +
	scale_color_manual(values = c("dodgerblue", "blue", "black")) + theme(legend.position = "none") +
	ylab(expression(paste("Biting rate (days"^-1,")"))) + xlab("Temperature (°C)") 
p.br

##################################
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

Temp.xs.2 <- seq(-6,56,.1) #For better mapping needs later changed the increment to 0.1
N.Temp.xs.2 <-length(Temp.xs.2)


# save the temperature sequence for future analyses
save(Temp.xs, file = "saved posteriors/temps.Rdata")#named identically as previous temps.Rdata
save(Temp.xs.2, file ="saved posteriors/temps.expanded.Rdata") #name differently to distinguish


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


################################# 3. Bayesian Fitting - biterate
#a - constant
#b - dtr9
#c - dtr12

#############################################3
######3a. biterate- constant temperature
##############################################
# Pull out data - *** Make sure Temperature units are correct (C or K as needed) ***

data <- c.data.block.temp %>%
  dplyr::filter(DTR == "constant") %>%
  dplyr::rename(trait = bite.rate,
         T = Treatment)

# Save as vectors for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)

# Fit the model
 br.c.quad_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
 				  model.file="quad_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
 				  n.iter=ni, DIC=T, working.directory=getwd())	

 br.c.briere_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                       model.file="briere_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                       n.iter=ni, DIC=T, working.directory=getwd())	
 br.c.neg.briere <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                         model.file="briere_neg.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                         n.iter=ni, DIC=T, working.directory=getwd())	
 
 br.c.neg.test.briere <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                              model.file="briere_neg_test.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                              n.iter=ni, DIC=T, working.directory=getwd())	
 #Remember T in Kelvin
 # Save as vectors for JAGS
 trait <- data$trait
 N.obs <- length(trait)
 temp <- data$T +273.15
 
 # Bundle all data in a list for JAGS
 jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)
 
 
 br.c.mod.SS <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                     model.file="SS_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                     n.iter=ni, DIC=T, working.directory=getwd())
 
mcmcplot(br.c.quad_T)
mcmcplot(br.c.briere_T)
mcmcplot(br.c.neg.briere)
mcmcplot(br.c.neg.test.briere)

br.c.quad_T$BUGSoutput$summary
br.c.briere_T$BUGSoutput$summary

br.c.quad_T$BUGSoutput$DIC
br.c.briere_T$BUGSoutput$DIC
br.c.neg.test.briere$BUGSoutput$DIC

# Calculate summary statistics of trait trajectories for plotting
br.c.quad_T.summary <- trait.trajs.quad(br.c.quad_T, Temp.xs, summary = TRUE)
br.c.briere_T.summary <- trait.trajs.briere(br.c.briere_T, Temp.xs, summary = TRUE)
br.c.neg.briere.summary <- trait.trajs.neg.briere(br.c.neg.briere, Temp.xs, summary = TRUE)
br.c.neg.test.briere.summary <- trait.trajs.neg.test.briere(br.c.neg.test.briere, Temp.xs, summary = TRUE)
br.c.mod.SS.summary <- trait.trajs.SS(br.c.mod.SS,(Temp.xs+273.15) , summary = TRUE)

# Save model output 
save(br.c.quad_T, file = "saved posteriors/constant_br_c.quad_T.uniform.Rdata")
save(br.c.briere_T, file = "saved posteriors/constant_br.c.briere_T.uniform.Rdata")
save(br.c.neg.briere, file = "saved posteriors/constant_br.c.neg.briere.uniform.Rdata")
save(br.c.neg.test.briere, file = "saved posteriors/constant_br.c.neg.test.briere.uniform.Rdata")

# Plot data and trait fit
p.br.quad_T.fit.c <- all.data.block.temp %>%
	filter(DTR == "constant") %>%
	ggplot() + 
	geom_point(aes(x = Treatment, y = bite.rate), color = "black", size = 3) +
  	geom_line(data = br.c.quad_T.summary, aes(x = temp, y = mean)) +
  	geom_line(data = br.c.quad_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  	geom_line(data = br.c.quad_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  	ylab("Bite rate (days^-1)") + xlab("Temperature (°C)") 
p.br.quad_T.fit.c 
ggsave("plots/biterate_c.br_quad_T.fit_only.pdf")

# Plot data and trait fit
p.br.briere_T.fit.c <- all.data.block.temp %>%
  filter(DTR == "constant") %>%
  ggplot() + 
  geom_point(aes(x = Treatment, y = bite.rate), color = "black", size = 3) +
  geom_line(data = br.c.briere_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = br.c.briere_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = br.c.briere_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Bite rate (days^-1)") + xlab("Temperature (°C)") 
p.br.briere_T.fit.c 
ggsave("plots/biterate_c.br_briere_T.fit_only.pdf")

# Plot data and trait fit
p.br.neg.briere_T.fit.c <- all.data.block.temp %>%
  filter(DTR == "constant") %>%
  ggplot() + 
  geom_point(aes(x = Treatment, y = bite.rate), color = "black", size = 3) +
  geom_line(data = br.c.neg.briere.summary, aes(x = temp, y = mean)) +
  geom_line(data = br.c.neg.briere.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = br.c.neg.briere.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Bite rate (days^-1)") + xlab("Temperature (°C)") 
p.br.neg.briere_T.fit.c 
ggsave("plots/biterate_c.neg.br_briere_T.fit_only.pdf")

# Plot data and trait fit #see if this function works....function is wrong
p.br.neg.test.briere_T.fit.c <- all.data.block.temp %>%
  filter(DTR == "constant") %>%
  ggplot() + 
  geom_point(aes(x = Treatment, y = bite.rate), color = "black", size = 3) +
  scale_y_continuous(breaks = seq(-0.5,1,0.2), limits = c(-1,1))+
  geom_line(data = br.c.neg.test.briere.summary, aes(x = temp, y = median)) +
  geom_line(data = br.c.neg.test.briere.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = br.c.neg.test.briere.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Bite rate (days^-1)") + xlab("Temperature (°C)") 
p.br.neg.test.briere_T.fit.c 

#plot(br.c.neg.test.briere.summary[,2]~Temp.xs) #none summary version #good to go

######################################
##########3b. lifespan dtr 9 models
#####################################

# Pull out data - *** Make sure Temperature units are correct (C or K as needed) ***

data <- f9.data.block.temp %>%
  dplyr::filter(DTR == "+/- 4.5") %>%
  rename(trait = bite.rate,
         T = Treatment)

# Save as vectors for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)

# Fit the model
br.f9.quad_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                    model.file="quad_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                    n.iter=ni, DIC=T, working.directory=getwd())	

br.f9.briere_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                      model.file="briere_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                      n.iter=ni, DIC=T, working.directory=getwd())	

mcmcplot(br.f9.quad_T)
mcmcplot(br.f9.briere_T)

br.f9.quad_T$BUGSoutput$summary
br.f9.briere_T$BUGSoutput$summary

br.f9.quad_T$BUGSoutput$DIC
br.f9.briere_T$BUGSoutput$DIC

# Calculate summary statistics of trait trajectories for plotting
br.f9.quad_T.summary <- trait.trajs.quad(br.f9.quad_T, Temp.xs, summary = TRUE)
br.f9.briere_T.summary <- trait.trajs.briere(br.f9.briere_T, Temp.xs, summary = TRUE)

# Save model output 
save(br.f9.quad_T, file = "saved posteriors/dtr_br_f9.quad_T.uniform.Rdata")
save(br.f9.briere_T, file = "saved posteriors/dtr_br.f9.briere_T.uniform.Rdata")

# Plot data and trait fit
p.br.quad_T.fit.f9 <- all.data.block.temp %>%
  filter(DTR == "+/- 4.5") %>%
  ggplot() + 
  geom_point(aes(x = Treatment, y = bite.rate), color = "dodgerblue", size = 3) +
  geom_line(data = br.f9.quad_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = br.f9.quad_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = br.f9.quad_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Bite rate (days^-1)") + xlab("Temperature (°C)") 
p.br.quad_T.fit.f9 
ggsave("plots/Biterate_f9.br_quad_T.fit_only.pdf")

# Plot data and trait fit
p.br.briere_T.fit.f9 <- all.data.block.temp %>%
  filter(DTR == "+/- 4.5") %>%
  ggplot() + 
  geom_point(aes(x = Treatment, y = bite.rate), color = "dodgerblue", size = 3) +
  geom_line(data = br.f9.briere_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = br.f9.briere_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = br.f9.briere_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Bite rate (days^-1)") + xlab("Temperature (°C)") 
p.br.briere_T.fit.f9 
ggsave("plots/biterate_f9.br_briere_T.fit_only.pdf")

##############################################################
#################3c - biterate - dtr 12
##########################################################

# Pull out data - *** Make sure Temperature units are correct (C or K as needed) ***

data <- f12.data.block.temp %>%
  dplyr::filter(DTR == "+/- 6") %>%
  rename(trait = bite.rate,
         T = Treatment)

# Save as vectors for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)

# Fit the model
br.f12.quad_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                     model.file="quad_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                     n.iter=ni, DIC=T, working.directory=getwd())	

br.f12.briere_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                       model.file="briere_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                       n.iter=ni, DIC=T, working.directory=getwd())	

mcmcplot(br.f12.quad_T)
mcmcplot(br.f12.briere_T)

br.f12.quad_T$BUGSoutput$summary
br.f12.briere_T$BUGSoutput$summary

br.f12.quad_T$BUGSoutput$DIC
br.f12.briere_T$BUGSoutput$DIC

# Calculate summary statistics of trait trajectories for plotting
br.f12.quad_T.summary <- trait.trajs.quad(br.f12.quad_T, Temp.xs, summary = TRUE)
br.f12.briere_T.summary <- trait.trajs.briere(br.f12.briere_T, Temp.xs, summary = TRUE)

# Save model output 
save(br.f12.quad_T, file = "saved posteriors/dtr_br_f12.quad_T.uniform.Rdata")
save(br.f12.briere_T, file = "saved posteriors/dtr_br.f12.briere_T.uniform.Rdata")

# Plot data and trait fit
p.br.quad_T.fit.f12 <- all.data.block.temp %>%
  filter(DTR == "+/- 6") %>%
  ggplot() + 
  geom_point(aes(x = Treatment, y = bite.rate), color = "royalblue", size = 3) +
  geom_line(data = br.f12.quad_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = br.f12.quad_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = br.f12.quad_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Bite rate (days^-1)") + xlab("Temperature (°C)") 
p.br.quad_T.fit.f12
ggsave("plots/biterate_f12.br_quad_T.fit_only.pdf")

# Plot data and trait fit
p.br.briere_T.fit.f12 <- all.data.block.temp %>%
  filter(DTR == "+/- 6") %>%
  ggplot() + 
  geom_point(aes(x = Treatment, y = bite.rate), color = "royalblue", size = 3) +
  geom_line(data = br.f12.briere_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = br.f12.briere_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = br.f12.briere_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Bite rate (days^-1)") + xlab("Temperature (°C)") 
p.br.briere_T.fit.f12
ggsave("plots/biterate_f12.br_briere_T.fit_only.pdf")
