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
#add in new EFD column -> total eggs laid by female divided by her lifespan
cdata$EFD.by.individual <- cdata$lifetime.eggs / cdata$lifespan
f9data$EFD.by.individual <- f9data$lifetime.eggs / f9data$lifespan
f12data$EFD.by.individual <- f12data$lifetime.eggs / f12data$lifespan


#change default plot specifications
par(mar=c(1,1,1,1))
dev.off()
###### Summarize + plot data for all 3 traits: lifespan, biting rate, fecundity
############## First plot each trait to check functional forms
boxplot(lifespan ~ Treatment, data = cdata) 
boxplot(bite.rate~Treatment, data = cdata)
boxplot(lifetime.eggs~Treatment, data = cdata)
boxplot(lifespan ~ Treatment, data = f9data) 
boxplot(bite.rate~Treatment, data = f9data)
boxplot(lifetime.eggs~Treatment, data = f9data)
boxplot(lifespan ~ Treatment, data = f12data) 
boxplot(bite.rate~Treatment, data = f12data)
boxplot(lifetime.eggs~Treatment, data = f12data)

#######Generate Block specific data values for fitting the curve over 
c.data.block.temp <- cdata %>%
  dplyr::group_by(Treatment,Block)%>%
  dplyr::summarise(lifespan = mean(lifespan),
                   bite.rate = mean(bite.rate),
                   lifetime.eggs = mean(lifetime.eggs),
                   EFD = mean(EFD.by.individual)) %>%
  ungroup()

plot(lifespan ~ Treatment, data = c.data.block.temp) 
plot(bite.rate ~ Treatment, data = c.data.block.temp) 
plot(lifetime.eggs ~ Treatment, data = c.data.block.temp) 
plot( EFD ~ Treatment, data = c.data.block.temp)

################dtr 9#################
f9.data.block.temp <- f9data %>%
  dplyr::group_by(Treatment,Block)%>%
  dplyr::summarise(lifespan = mean(lifespan),
                   bite.rate = mean(bite.rate),
                   lifetime.eggs = mean(lifetime.eggs),
                   EFD = mean(EFD.by.individual)) %>%
  ungroup()

plot(lifespan ~ Treatment, data = f9.data.block.temp) 
plot(bite.rate ~ Treatment, data = f9.data.block.temp) 
plot(lifetime.eggs ~ Treatment, data = f9.data.block.temp) 
plot(EFD ~ Treatment, data = f9.data.block.temp)
###############dtr 12##################3
f12.data.block.temp <- f12data %>%
  dplyr::group_by(Treatment,Block)%>%
  dplyr::summarise(lifespan = mean(lifespan),
                   bite.rate = mean(bite.rate),
                   lifetime.eggs = mean(lifetime.eggs),
                   EFD = mean(EFD.by.individual)) %>%
  ungroup()

plot(lifespan ~ Treatment, data = f12.data.block.temp) 
plot(bite.rate ~ Treatment, data = f12.data.block.temp) 
plot(lifetime.eggs ~ Treatment, data = f12.data.block.temp) 
plot(EFD ~ Treatment, data = f12.data.block.temp)



# Add in DTR Treatment columns, combine in a single dataframe
c.data.block.temp$DTR <- "constant"
f9.data.block.temp$DTR <- "+/- 4.5"
f12.data.block.temp$DTR <- "+/- 6"
all.data.block.temp <- rbind(c.data.block.temp, f9.data.block.temp, f12.data.block.temp)

all.data.block.temp <- all.data.block.temp %>% 
           dplyr::mutate(TempK = Treatment + 273.15)

#write.csv(all.data.block.temp, "fluc.program.trait.means.csv", row.names = FALSE)

# Plot lifespan
plot(c.data.block.temp$Treatment, c.data.block.temp$lifespan, ylim = c(0,60), ylab = "Lifespan (days)", xlab = ("Temperature (°C)"), pch =16, col = "black")
points(f9.data.block.temp$Treatment, f9.data.block.temp$lifespan, pch =16, col = "dodgerblue")
points(f12.data.block.temp$Treatment, f12.data.block.temp$lifespan, pch =16, col = "blue")

p.lf <- ggplot() + 
	geom_point(data = all.data.block.temp, aes(x = jitter(Treatment,0.3), y = lifespan, color = DTR), size = 2) +
  scale_color_manual(values = c("dodgerblue", "blue", "black")) + theme(legend.position = c(0.7,0.8), legend.text=element_text(size=12)) +
	ylab("Lifespan (days)") + xlab("Temperature (°C)") 
p.lf

# Plot biting rate
p.br <- ggplot() + 
	geom_point(data = all.data.block.temp, aes(x = jitter(Treatment,0.3), y = bite.rate, color = DTR), size =2) +
	scale_color_manual(values = c("dodgerblue", "blue", "black")) + theme(legend.position = "none") +
	ylab(expression(paste("Biting rate (days"^-1,")"))) + xlab("Temperature (°C)") 
p.br

# Plot fecundity
p.f <- ggplot() + 
	geom_point(data = all.data.block.temp, aes(x = jitter(Treatment,0.3), y = lifetime.eggs, color = DTR), size =2) +
  scale_color_manual(values = c("dodgerblue", "blue", "black")) + theme(legend.position = "none") +
	ylab("Lifetime fecundity (eggs)") + xlab("Temperature (°C)") 

par(mfrow = c(3,3))
plot_grid(p.lf, p.br, p.f, nrow = 1)
ggsave("plots/TraitMeansData.pdf",   width = 8, height = 2.5, units = "in" ) #Saving 11.7 x 4.24 in image
dev.off()

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


### Mod SS_T function parameters and inits
#Parameters to Estimate
parameters <- c("cf.c", "cf.Topt", "cf.ED", "cf.sigma")
inits <- function(){list(
  cf.c = 10,
  cf.Topt = 300,
  cf.ED = 1.6,
  cf.sigma = rlnorm(1))}  

#### Norberg function parameters and inits
#Parameters to Estimate
parameters <- c("cf.z", "cf.w", "cf.a", "cf.b")
inits <- function(){list(
  cf.z = 20,
  cf.w = 20,
  cf.a = 5,
  cf.b =  .5 )}  

################################# 3. Bayesian Fitting - lifespan
#a - constant
#b - dtr9
#c - dtr12

#############################################3
######3a. Lifespan - constant temperature
##############################################
# Pull out data - *** Make sure Temperature units are correct (C or K as needed) ***

data <- c.data.block.temp %>%
  dplyr::filter(DTR == "constant") %>%
  rename(trait = lifespan,
         T = Treatment)

data <- c.data.block.temp %>%
  dplyr::filter(DTR == "constant") %>%
  rename(trait = bite.rate,
         T = Treatment)


# Save as vectors for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)

# Fit the model
 lf.c.quad_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
 				  model.file="quad_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
 				  n.iter=ni, DIC=T, working.directory=getwd())	

 lf.c.briere_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                       model.file="briere_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                       n.iter=ni, DIC=T, working.directory=getwd())	
 
 ############For the more thermally relevant equations
 c.data.block.temp$Temp.K <- c.data.block.temp$Treatment + 273.15
 data <- c.data.block.temp %>%
   dplyr::filter(DTR == "constant") %>%
   rename(trait = lifespan,
          T = Temp.K)
 
f9.data.block.temp$Temp.K <- f9.data.block.temp$Treatment + 273.15
 data <- f9.data.block.temp %>%
   dplyr::filter(DTR == "+/- 4.5") %>%
   rename(trait = lifespan,
          T = Temp.K)
 
 f12.data.block.temp$Temp.K <- f12.data.block.temp$Treatment + 273.15
 data <- f12.data.block.temp %>%
   dplyr::filter(DTR == "+/- 6") %>%
   rename(trait = lifespan,
          T = Temp.K)
 
 c.data.block.temp$Temp.K <- c.data.block.temp$Treatment + 273.15
 data <- c.data.block.temp %>%
   dplyr::filter(DTR == "constant") %>%
   rename(trait = bite.rate,
          T = Treatment)
 
 
 #Remember T in Kelvin
 # Save as vectors for JAGS
 trait <- data$trait
 N.obs <- length(trait)
 temp <- data$T
 
 # Bundle all data in a list for JAGS
 jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)
 
 
 lf.c.SS_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                     model.file="SS_T_KMadjust.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                     n.iter=ni, DIC=T, working.directory=getwd())
 
 lf.c.Norberg_KM <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                         model.file="norberg_KMadjust.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                         n.iter=ni, DIC=T, working.directory=getwd())
 lf.c.Norberg_KM2 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                         model.file="norberg_KMadjust2.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                         n.iter=ni, DIC=T, working.directory=getwd())
 
 br.c.Norberg_KM <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                         model.file="norberg_KMadjust.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                         n.iter=ni, DIC=T, working.directory=getwd())
 br.c.Norberg_KM2 <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                         model.file="norberg_KMadjust2.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                         n.iter=ni, DIC=T, working.directory=getwd())
 
 
 
 lf.f9.SS_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                   model.file="SS_T_KMadjust.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                   n.iter=ni, DIC=T, working.directory=getwd())
 
 lf.f12.SS_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                    model.file="SS_T_KMadjust.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                    n.iter=ni, DIC=T, working.directory=getwd())
 
mcmcplot(lf.c.SS_T) # worked with SS_T_KMadjust model spec = ED (1-10)*Prevents log term from going negative; inits as 1.6
lf.c.SS_T$BUGSoutput$summary #means ED = 2;Topt = 294; c = 10; sigma = 3.6
lf.c.SS_T$BUGSoutput$DIC #75.88615
lf.c.SS_T.summary <- trait.trajs.SS_T_KMadjust(lf.c.SS_T, Temp.xs.K, summary = TRUE)
save(lf.c.SS_T, file = "saved posteriors/constant_lf.c.SS_T.Rdata")

mcmcplot(lf.c.Norberg_KM2)
lf.c.Norberg_KM$BUGSoutput$summary #means a = 17, b = 0.03, w=28,z=21
lf.c.Norberg_KM2$BUGSoutput$DIC #98
lf.c.Norberg_KM2.summary <- trait.trajs.Norberg_KMadjust(lf.c.Norberg_KM2, Temp.xs, summary = TRUE)

mcmcplot(br.c.Norberg_KM2)
br.c.Norberg_KM2$BUGSoutput$summary #means a = -0.01, b = 0.003, w=17.9,z=44
br.c.Norberg_KM2$BUGSoutput$DIC #-35.277
#br.c.Norberg_KM.summary <- trait.trajs.Norberg_KMadjust(br.c.Norberg_KM, Temp.xs.K, summary = TRUE)
br.c.Norberg_KM2.summary <- trait.trajs.Norberg_KMadjust(br.c.Norberg_KM2, Temp.xs, summary = TRUE)




mcmcplot(lf.f9.SS_T)
lf.f9.SS_T$BUGSoutput$summary #Means ED = 1.45, Topt = 287, c = 24.9, sigma - 4.9
lf.f9.SS_T$BUGSoutput$DIC #65.24
lf.f9.SS_T.summary <- trait.trajs.SS_T_KMadjust(lf.f9.SS_T, Temp.xs.K, summary = TRUE)
save(lf.f9.SS_T, file = "saved posteriors/fluctuation_lf.f9.SS_T.Rdata")

mcmcplot(lf.f12.SS_T)
lf.f12.SS_T$BUGSoutput$summary #Means ED = 1.72, Topt = 287, c = 27.7, sigma - 8.3
lf.f12.SS_T$BUGSoutput$DIC #76.4
lf.f12.SS_T.summary <- trait.trajs.SS_T_KMadjust(lf.f12.SS_T, Temp.xs.K, summary = TRUE)
save(lf.f12.SS_T, file = "saved posteriors/fluctuation_lf.f12.SS_T.Rdata")


# Plot data and trait fit
p.lf.c.SS_T_KMadjust_fit <- all.data.block.temp %>%
  filter(DTR == "constant") %>%
  ggplot() + 
  geom_point(aes(x = Treatment, y = lifespan), color = "black", size = 3) +
  geom_line(data = lf.c.SS_T.summary, aes(x = temp-273.15, y = mean)) +
  geom_line(data = lf.c.SS_T.summary, aes(x = temp-273.15, y = upper), linetype = "dashed") +
  geom_line(data = lf.c.SS_T.summary, aes(x = temp-273.15, y = lower), linetype = "dashed") +
  ylab("Lifespan (days)") + xlab("Temperature (°C)") 
p.lf.c.SS_T_KMadjust_fit 
ggsave("plots/lifespan_c.lf_SS_T_KMadjust_only.pdf")

p.lf.c.Norberg_KM_fit <- all.data.block.temp %>%
  filter(DTR == "constant") %>%
  ggplot() + 
  geom_point(aes(x = Treatment, y = lifespan), color = "black", size = 3) +
  scale_y_continuous(limits= c(-30, 60))+
  geom_line(data = lf.c.Norberg_KM2.summary, aes(x = temp, y = mean)) +
  geom_line(data = lf.c.Norberg_KM2.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = lf.c.Norberg_KM2.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Lifespan (days)") + xlab("Temperature (°C)") 
p.lf.c.Norberg_KM_fit 
ggsave("plots/lifespan_c.lf_Norberg_KM_only.pdf")


p.br.c.Norberg_KM_fit2 <- all.data.block.temp %>%
  filter(DTR == "constant") %>%
  ggplot() + 
  geom_point(aes(x = Treatment, y = bite.rate), color = "black", size = 3) +
  scale_y_continuous(limits = c(-1,1))+
  geom_line(data = br.c.Norberg_KM2.summary, aes(x = temp, y = mean)) +
  geom_line(data = br.c.Norberg_KM2.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = br.c.Norberg_KM2.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Biterate (days)") + xlab("Temperature (°C)") 
p.br.c.Norberg_KM_fit2 
ggsave("plots/biterate_c.br_Norberg_KM_only.pdf")



# Plot data and trait fit
p.lf.f9.SS_T_KMadjust_fit <- all.data.block.temp %>%
  filter(DTR == "+/- 4.5") %>%
  ggplot() + 
  geom_point(aes(x = Treatment, y = lifespan), color = "black", size = 3) +
  geom_line(data = lf.f9.SS_T.summary, aes(x = temp-273.15, y = mean)) +
  geom_line(data = lf.f9.SS_T.summary, aes(x = temp-273.15, y = upper), linetype = "dashed") +
  geom_line(data = lf.f9.SS_T.summary, aes(x = temp-273.15, y = lower), linetype = "dashed") +
  ylab("Lifespan (days)") + xlab("Temperature (°C)") 
p.lf.f9.SS_T_KMadjust_fit 
ggsave("plots/lifespan_f9.lf_SS_T_KMadjust_only.pdf")

# Plot data and trait fit
p.lf.f12.SS_T_KMadjust_fit <- all.data.block.temp %>%
  filter(DTR == "+/- 6") %>%
  ggplot() + 
  geom_point(aes(x = Treatment, y = lifespan), color = "black", size = 3) +
  geom_line(data = lf.f12.SS_T.summary, aes(x = temp-273.15, y = mean)) +
  geom_line(data = lf.f12.SS_T.summary, aes(x = temp-273.15, y = upper), linetype = "dashed") +
  geom_line(data = lf.f12.SS_T.summary, aes(x = temp-273.15, y = lower), linetype = "dashed") +
  ylab("Lifespan (days)") + xlab("Temperature (°C)") 
p.lf.f12.SS_T_KMadjust_fit 
ggsave("plots/lifespan_f12.lf_SS_T_KMadjust_only.pdf")



mcmcplot(lf.c.quad_T)
mcmcplot(lf.c.briere_T)

lf.c.quad_T$BUGSoutput$summary
lf.c.briere_T$BUGSoutput$summary

lf.c.quad_T$BUGSoutput$DIC
lf.c.briere_T$BUGSoutput$DIC

# Calculate summary statistics of trait trajectories for plotting
lf.c.quad_T.summary <- trait.trajs.quad(lf.c.quad_T, Temp.xs, summary = TRUE)
lf.c.briere_T.summary <- trait.trajs.briere(lf.c.briere_T, Temp.xs, summary = TRUE)

# Save model output 
save(lf.c.quad_T, file = "saved posteriors/constant_lf_c.quad_T.uniform.Rdata")
save(lf.c.briere_T, file = "saved posteriors/constant_lf.c.briere_T.uniform.Rdata")

# Plot data and trait fit
p.lf.quad_T.fit.c <- all.data.block.temp %>%
	filter(DTR == "constant") %>%
	ggplot() + 
	geom_point(aes(x = Treatment, y = lifespan), color = "black", size = 3) +
  	geom_line(data = lf.c.quad_T.summary, aes(x = temp, y = mean)) +
  	geom_line(data = lf.c.quad_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  	geom_line(data = lf.c.quad_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  	ylab("Lifespan (days)") + xlab("Temperature (°C)") 
p.lf.quad_T.fit.c 
ggsave("plots/lifespan_c.lf_quad_T.fit_only.pdf")

# Plot data and trait fit
p.lf.briere_T.fit.c <- all.data.block.temp %>%
  filter(DTR == "constant") %>%
  ggplot() + 
  geom_point(aes(x = Treatment, y = lifespan), color = "black", size = 3) +
  geom_line(data = lf.c.briere_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = lf.c.briere_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = lf.c.briere_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Lifespan (days)") + xlab("Temperature (°C)") 
p.lf.briere_T.fit.c 
ggsave("plots/lifespan_c.lf_briere_T.fit_only.pdf")

######################################
##########3b. lifespan dtr 9 models
#####################################

# Pull out data - *** Make sure Temperature units are correct (C or K as needed) ***

data <- f9.data.block.temp %>%
  dplyr::filter(DTR == "+/- 4.5") %>%
  rename(trait = lifespan,
         T = Treatment)

# Save as vectors for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)

# Fit the model
lf.f9.quad_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                    model.file="quad_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                    n.iter=ni, DIC=T, working.directory=getwd())	

lf.f9.briere_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                      model.file="briere_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                      n.iter=ni, DIC=T, working.directory=getwd())	

mcmcplot(lf.f9.quad_T)
mcmcplot(lf.f9.briere_T)

lf.f9.quad_T$BUGSoutput$summary
lf.f9.briere_T$BUGSoutput$summary

lf.f9.quad_T$BUGSoutput$DIC
lf.f9.briere_T$BUGSoutput$DIC

# Calculate summary statistics of trait trajectories for plotting
lf.f9.quad_T.summary <- trait.trajs.quad(lf.f9.quad_T, Temp.xs, summary = TRUE)
lf.f9.briere_T.summary <- trait.trajs.briere(lf.f9.briere_T, Temp.xs, summary = TRUE)

# Save model output 
save(lf.f9.quad_T, file = "saved posteriors/dtr_lf_f9.quad_T.uniform.Rdata")
save(lf.f9.briere_T, file = "saved posteriors/dtr_lf.f9.briere_T.uniform.Rdata")

# Plot data and trait fit
p.lf.quad_T.fit.f9 <- all.data.block.temp %>%
  filter(DTR == "+/- 4.5") %>%
  ggplot() + 
  geom_point(aes(x = Treatment, y = lifespan), color = "dodgerblue", size = 3) +
  geom_line(data = lf.f9.quad_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = lf.f9.quad_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = lf.f9.quad_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Lifespan (days)") + xlab("Temperature (°C)") 
p.lf.quad_T.fit.f9 
ggsave("plots/lifespan_f9.lf_quad_T.fit_only.pdf")

# Plot data and trait fit
p.lf.briere_T.fit.f9 <- all.data.block.temp %>%
  filter(DTR == "+/- 4.5") %>%
  ggplot() + 
  geom_point(aes(x = Treatment, y = lifespan), color = "dodgerblue", size = 3) +
  geom_line(data = lf.f9.briere_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = lf.f9.briere_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = lf.f9.briere_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Lifespan (days)") + xlab("Temperature (°C)") 
p.lf.briere_T.fit.f9 
ggsave("plots/lifespan_f9.lf_briere_T.fit_only.pdf")

##############################################################
#################3c - Lifespan - dtr 12
##########################################################

# Pull out data - *** Make sure Temperature units are correct (C or K as needed) ***

data <- f12.data.block.temp %>%
  dplyr::filter(DTR == "+/- 6") %>%
  rename(trait = lifespan,
         T = Treatment)

# Save as vectors for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)

# Fit the model
lf.f12.quad_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                     model.file="quad_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                     n.iter=ni, DIC=T, working.directory=getwd())	

lf.f12.briere_T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                       model.file="briere_T.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
                       n.iter=ni, DIC=T, working.directory=getwd())	

mcmcplot(lf.f12.quad_T)
mcmcplot(lf.f12.briere_T)

lf.f12.quad_T$BUGSoutput$summary
lf.f12.briere_T$BUGSoutput$summary

lf.f12.quad_T$BUGSoutput$DIC
lf.f12.briere_T$BUGSoutput$DIC

# Calculate summary statistics of trait trajectories for plotting
lf.f12.quad_T.summary <- trait.trajs.quad(lf.f12.quad_T, Temp.xs, summary = TRUE)
lf.f12.briere_T.summary <- trait.trajs.briere(lf.f12.briere_T, Temp.xs, summary = TRUE)

# Save model output 
save(lf.f12.quad_T, file = "saved posteriors/dtr_lf_f12.quad_T.uniform.Rdata")
save(lf.f12.briere_T, file = "saved posteriors/dtr_lf.f12.briere_T.uniform.Rdata")

# Plot data and trait fit
p.lf.quad_T.fit.f12 <- all.data.block.temp %>%
  filter(DTR == "+/- 6") %>%
  ggplot() + 
  geom_point(aes(x = Treatment, y = lifespan), color = "royalblue", size = 3) +
  geom_line(data = lf.f12.quad_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = lf.f12.quad_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = lf.f12.quad_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Lifespan (days)") + xlab("Temperature (°C)") 
p.lf.quad_T.fit.f12
ggsave("plots/lifespan_f12.lf_quad_T.fit_only.pdf")

# Plot data and trait fit
p.lf.briere_T.fit.f12 <- all.data.block.temp %>%
  filter(DTR == "+/- 6") %>%
  ggplot() + 
  geom_point(aes(x = Treatment, y = lifespan), color = "royalblue", size = 3) +
  geom_line(data = lf.f12.briere_T.summary, aes(x = temp, y = mean)) +
  geom_line(data = lf.f12.briere_T.summary, aes(x = temp, y = upper), linetype = "dashed") +
  geom_line(data = lf.f12.briere_T.summary, aes(x = temp, y = lower), linetype = "dashed") +
  ylab("Lifespan (days)") + xlab("Temperature (°C)") 
p.lf.briere_T.fit.f12
ggsave("plots/lifespan_f12.lf_briere_T.fit_only.pdf")
