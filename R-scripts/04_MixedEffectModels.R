## Kerri Miazgowicz, University of Georgia
## October 1, 2018
## Stats included October 1, 2019
## Rerun 4/22/2020
## Additions by Marta Shocket Nov 2021
## 
## Purpose: Generate indiviudal-specific values for lifespan, daily biting rate, and lifetime egg production for constant and fluctuating temperatures
## ID / temp / a / lf / B <- from master raw csv files
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) Calculate a / lf / B across all datasets
##           3) Export files to be used in TraitFitting.R -> provided to Joey and Marta for fitting in code file ''.
##           4) Perform GLMM models on a / lf/ B across all datasets 
##        


#####################
###
#  1) Set-up,load packages, get data, etc.
###
######################
library(dplyr)

#Set working directory
#mainDir = "C:/Users/Kerri/Desktop/Chapter2 InProgress"
#setwd(mainDir)
getwd()

#Read in csv files for each dataset
constant.data <- read.csv("data/constant_master.csv")
constant.data[is.na(constant.data)] <- 0
fluc9.data <- read.csv("data/fluc_dtr9_master.csv")
fluc9.data[is.na(fluc9.data)] <- 0
fluc12.data <- read.csv("data/fluc_dtr12_master.csv")
fluc12.data[is.na(fluc12.data)] <- 0


#####################
###
#  2) Calculate a / lf / B for each individual across all datasets
###
######################

#create new dataframe for each temperature profile for the aggregated data

constant.individual.trait <- constant.data %>%
	dplyr::group_by(Female,Treatment,Block)%>%
	dplyr::summarise(lifespan = max(Day), bite.rate = sum(Feed)/ lifespan, lifetime.eggs = sum(Count))

fluc9.individual.trait <- fluc9.data %>%
	dplyr::group_by(Female,Treatment, Block)%>%
	dplyr::summarise(lifespan = max(Day), bite.rate = sum(Feed)/ lifespan, lifetime.eggs = sum(Count))


fluc12.individual.trait <- fluc12.data%>% 
	dplyr::group_by(Female,Treatment, Block)%>%
	dplyr::summarise(lifespan = max(Day), bite.rate = sum(Feed)/ lifespan, lifetime.eggs = sum(Count))




#####################
###
#  3) Export files to be used in TraitFitting.R
###
######################

#write.csv(constant.individual.trait, "data/constant.individual.trait.csv", row.names = F)
#write.csv(fluc9.individual.trait, "data/fluc9.individual.trait.csv", row.names = F)
#write.csv(fluc12.individual.trait, "data/fluc12.individual.trait.csv", row.names =F)

#####################
###
#  4) Perform GLMM models on a / lf/ B across all datasets 
### Note: Need to keep block variable for random effects; won't be able to include donor or day due to summary data structure
#####################
#

#Repeat above summary, but include Block for analysis
#create new dataframe for each temperature profile for the aggregated data

constant.individual.trait <- as.data.frame(summarise(group_by(constant.data, Female,Treatment, Block), lifespan = max(Day), bite.rate = sum(Feed)/ lifespan, lifetime.eggs = sum(Count)))
fluc9.individual.trait <- as.data.frame(summarise(group_by(fluc9.data, Female,Treatment, Block), lifespan = max(Day), bite.rate = sum(Feed)/ lifespan, lifetime.eggs = sum(Count)))
fluc12.individual.trait <- as.data.frame(summarise(group_by(fluc12.data, Female,Treatment, Block), lifespan = max(Day), bite.rate = sum(Feed)/ lifespan, lifetime.eggs = sum(Count)))

#replace "block" and "individual" values within each dataset with '.DTR0, 9, or 12' so that when I merge all data together those stay seperate; will become string value
#convert block to a string variable and then add a character segment to it
constant.individual.trait$Block <- paste0(as.character(constant.individual.trait$Block),".dtr0")
constant.individual.trait$Program <- as.character("dtr0")
fluc9.individual.trait$Block <- paste0(as.character(fluc9.individual.trait$Block), ".dtr9")
fluc9.individual.trait$Program <- as.character("dtr9")
fluc12.individual.trait$Block <- paste0(as.character(fluc12.individual.trait$Block), ".dtr12")
fluc12.individual.trait$Program <- as.character("dtr12")

#merge all data into one dataframe
all.individual.trait <- rbind(constant.individual.trait, fluc9.individual.trait, fluc12.individual.trait)

#perform GLMM analysis
#Main effects = Treatment(scaled & centered), Program (categorical)
#Random = Block (categorical), Female(categorical)
#Response variables = lifespan(continuous/count), bite.rate(rate), lifetime.eggs(continuous/count)

#Following Courtney's previous code[Tesla temperature], making Temperature a continuous variable
#feeddata$Treatment <- as.factor(feeddata$Treatment)
#all.individual.trait$Female <- as.factor(all.individual.trait$Female)
all.individual.trait$Block <- as.factor(all.individual.trait$Block)
all.individual.trait$Program <- as.factor(all.individual.trait$Program)

#Centering and scaling continuous variables of interest (Treatment)
all.individual.trait$Tscale <- scale(all.individual.trait$Treatment)

library(fitdistrplus)
library(lme4)
library(car)
library(glmmTMB) # can do zero-inflated GLMMs
library(broom.mixed) # Extracts residual deviance and df from GLMMs
library(insight) # Extracts variances from GLMMs
library(r2glmm) # Calculates residual r^2 for GLMMs
library(emmeans) # Estimates marginal means (for post hoc pairwise comparisons)
library(gridExtra)

#Looking at my data distribution
hist(all.individual.trait$lifespan, breaks=50, col="red") # right-skewed
hist(all.individual.trait$bite.rate, breaks=50, col="red") #normal between 0 and 1
hist(all.individual.trait$lifetime.eggs, breaks=50, col="red") #0-inflated nbinomial
#descdist(all.individual.trait$lifetime.eggs, discrete = FALSE) # This suggests to use a uniform distribution



#####################
#####################  Lifespan Analysis
#####################

##################### GLMMs

# Gamma distribution fits well
fit.gamma <- fitdist(all.individual.trait$lifespan, distr = "gamma", method = "mme")
plot(fit.gamma)

# Scale (divide by 10) so that model will fit without errors, and R^2 package will work without errors (doesn't work with formula in model call)
all.individual.trait$lifespan.t <- all.individual.trait$lifespan/10

# GLMMs for lifespan - divided by 10 to deal with scale issues - does not affect model validity, but will affect model coefficients
# See trouble-shooting methods below that did not work
lf.1 <- glmer(lifespan.t ~ (1|Block), family = Gamma(link = 'inverse'), data = all.individual.trait)
lf.2 <- glmer(lifespan.t ~ Tscale + (1|Block), family =Gamma (link = 'inverse'), data = all.individual.trait)
lf.3 <- glmer(lifespan.t ~ Tscale + I(Tscale^2) + (1|Block), family = Gamma(link = 'inverse'), data = all.individual.trait)
lf.4 <- glmer(lifespan.t ~ Program + (1|Block), family = Gamma(link = 'inverse'), data = all.individual.trait)
lf.5 <- glmer(lifespan.t ~ Tscale + Program + (1|Block), family = Gamma(link = 'inverse'), data = all.individual.trait)
lf.6 <- glmer(lifespan.t ~ Tscale + Program + I(Tscale^2) + (1|Block), family = Gamma(link = 'inverse'), data = all.individual.trait)
lf.7 <- glmer(lifespan.t ~ Tscale*Program + (1|Block), family = Gamma(link = 'inverse'), data = all.individual.trait)
lf.8 <- glmer(lifespan.t ~ Tscale*Program + I(Tscale^2) + (1|Block), family = Gamma(link = 'inverse'), data = all.individual.trait)
lf.9 <- glmer(lifespan.t ~ Tscale*Program + I(Tscale^2)*Program + (1|Block), family = Gamma(link = 'inverse'), data = all.individual.trait)

# Calculate delta AIC and AIC weights - remember to round AIC to 0.1 and AICw to 0.01 when reporting
lf.AIC.vec <- c(AIC(lf.1), AIC(lf.2), AIC(lf.3), AIC(lf.4), AIC(lf.5), AIC(lf.6), AIC(lf.7), AIC(lf.8), AIC(lf.9))
lf.dAIC.vec <- lf.AIC.vec - min(lf.AIC.vec)
lf.AICw.vec <- exp(-0.5*lf.dAIC.vec)/sum(exp(-0.5*lf.dAIC.vec))

# Get deviance and residual df for reporting stats according to: https://www.ashander.info/posts/2015/10/model-selection-glms-aic-what-to-report/
lf.table <- do.call(rbind, lapply(list(lf.1, lf.2, lf.3, lf.4, lf.5, lf.6, lf.7, lf.8, lf.9), broom::glance))

# Get marginal R^2 value (variance explained by fixed effects) - This doesn't work - need to troubleshoot further
lf.sR2.list <- lapply(list(lf.2, lf.3, lf.4, lf.5, lf.6, lf.7, lf.8, lf.9), r2beta, method = "nsj", partial = T)

# Examine values for manuscript table
lf.dAIC.vec
lf.AICw.vec
lf.table
lf.sR2.list

# Summary of best model
summary(lf.8)


# Troubleshooting original GLMM (with unscaled lifespan response variable) using:
#https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html

# Check for singularity
tt <- getME(lf.m2,"theta")
ll <- getME(lf.m2,"lower")
min(tt[ll==0])

# Check gradient calculations - this is below the typical tolerance (0.001), but so is the one for lifetime eggs that doesn't throw an error
derivs1 <- lf.m2@optinfo$derivs
grad1 <- with(derivs1,solve(Hessian,gradient))
max(abs(grad1))

# Start with previous fit - doesn't help
ss <- getME(lf.m2,c("theta","fixef"))
lf.m2b <- update(lf.m2,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))

# Try a different optimizer - doesn't help
lf.m2c <- update(lf.m2,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))


##################### Post-hoc Pairwise Comparisons

##### Constrasts from global model (averaged across all temperature groups)
emmeans(lf.8, pairwise ~ Program)

##### Constrasts with temperature treatments
# Subset data
individual.trait.16 <- subset(all.individual.trait, Treatment == 16)
individual.trait.20 <- subset(all.individual.trait, Treatment == 20)
individual.trait.24 <- subset(all.individual.trait, Treatment == 24)
individual.trait.28 <- subset(all.individual.trait, Treatment == 28)
individual.trait.32 <- subset(all.individual.trait, Treatment == 32)

# Fit new GLMMs for each mean temperature
lf.16 <- glm(lifespan.t ~ Program + Block, family = Gamma(link = 'inverse'), data = individual.trait.16)
lf.20 <- glm(lifespan.t ~ Program + Block, family = Gamma(link = 'inverse'), data = individual.trait.20)
lf.24 <- glm(lifespan.t ~ Program + Block, family = Gamma(link = 'inverse'), data = individual.trait.24)
lf.28 <- glm(lifespan.t ~ Program + Block, family = Gamma(link = 'inverse'), data = individual.trait.28)
lf.32 <- glm(lifespan.t ~ Program + Block, family = Gamma(link = 'inverse'), data = individual.trait.32)

summary(lf.16)
summary(lf.20)
summary(lf.24)
summary(lf.28)
summary(lf.32)

emmeans.lf.16 <- emmeans(lf.16, pairwise ~ Program)
emmeans.lf.20 <- emmeans(lf.20, pairwise ~ Program)
emmeans.lf.24 <- emmeans(lf.24, pairwise ~ Program)
emmeans.lf.28 <- emmeans(lf.28, pairwise ~ Program)
emmeans.lf.32 <- emmeans(lf.32, pairwise ~ Program)

p.lf.16 <- plot(emmeans.lf.16)
p.lf.20 <- plot(emmeans.lf.20)
p.lf.24 <- plot(emmeans.lf.24)
p.lf.28 <- plot(emmeans.lf.28)
p.lf.32 <- plot(emmeans.lf.32)
grid.arrange(p.lf.16 + labs(title = "16C"), p.lf.20 + labs(title = "20C"), p.lf.24 + labs(title = "24C"),
			 p.lf.28 + labs(title = "28C"), p.lf.32 + labs(title = "32C"))


#####################
#####################  Biting Rate Analysis
#####################

##################### GLMMs

# Gamma distribution fits well
fit.gamma2 <- fitdist(all.individual.trait$bite.rate, distr = "gamma", method = "mme")
plot(fit.gamma2)

# Add a small constant to create no zero values
all.individual.trait$bite.rate.t <- all.individual.trait$bite.rate + 0.0001

# GLMMs for biting rate
br.1 <- glmer(bite.rate.t ~ (1|Block), family = Gamma(link = 'inverse'), data = all.individual.trait)
br.2 <- glmer(bite.rate.t ~ Tscale + (1|Block), family =Gamma (link = 'inverse'), data = all.individual.trait)
br.3 <- glmer(bite.rate.t ~ Tscale + I(Tscale^2) + (1|Block), family = Gamma(link = 'inverse'), data = all.individual.trait)
br.4 <- glmer(bite.rate.t ~ Program + (1|Block), family = Gamma(link = 'inverse'), data = all.individual.trait)
br.5 <- glmer(bite.rate.t ~ Tscale + Program + (1|Block), family = Gamma(link = 'inverse'), data = all.individual.trait)
br.6 <- glmer(bite.rate.t ~ Tscale + Program + I(Tscale^2) + (1|Block), family = Gamma(link = 'inverse'), data = all.individual.trait)
br.7 <- glmer(bite.rate.t ~ Tscale*Program + (1|Block), family = Gamma(link = 'inverse'), data = all.individual.trait)
br.8 <- glmer(bite.rate.t ~ Tscale*Program + I(Tscale^2) + (1|Block), family = Gamma(link = 'inverse'), data = all.individual.trait)
br.9 <- glmer(bite.rate.t ~ Tscale*Program + I(Tscale^2)*Program + (1|Block), family = Gamma(link = 'inverse'), data = all.individual.trait)

# Calculate delta AIC and AIC weights - remember to round AIC to 0.1 and AICw to 0.01 when reporting
br.AIC.vec <- c(AIC(br.1), AIC(br.2), AIC(br.3), AIC(br.4), AIC(br.5), AIC(br.6), AIC(br.7), AIC(br.8), AIC(br.9))
br.dAIC.vec <- br.AIC.vec - min(br.AIC.vec)
br.AICw.vec <- exp(-0.5*br.dAIC.vec)/sum(exp(-0.5*br.dAIC.vec))

# Get deviance and residual df for reporting stats according to: https://www.ashander.info/posts/2015/10/model-selection-glms-aic-what-to-report/
br.table <- do.call(rbind, lapply(list(br.1, br.2, br.3, br.4, br.5, br.6, br.7, br.8, br.9), broom::glance))

# Get marginal R^2 value (variance explained by fixed effects)
br.sR2.list <- lapply(list(br.2, br.3, br.4, br.5, br.6, br.7, br.8, br.9), r2beta, method = "nsj", partial = T)

# Examine values for manuscript table
br.dAIC.vec
br.AICw.vec
br.table
br.sR2.list

# Summary of best model
summary(br.6)


##################### Post-hoc Pairwise 

# Making 15 total comparisons: 3 comparisons (DTR - 0 vs 9, 0 vs 12, 9 vs 12) within each of 5 mean temps
coef(br.6)
vcov(br.6)

emmeans(br.6, pairwise ~ Program)

# Subset data
individual.trait.16 <- subset(all.individual.trait, Treatment == 16)
individual.trait.20 <- subset(all.individual.trait, Treatment == 20)
individual.trait.24 <- subset(all.individual.trait, Treatment == 24)
individual.trait.28 <- subset(all.individual.trait, Treatment == 28)
individual.trait.32 <- subset(all.individual.trait, Treatment == 32)

# Fit new GLMMs for each mean temperature
br.16 <- glm(bite.rate.t ~ Program + Block, family = Gamma(link = 'inverse'), data = individual.trait.16)
br.20 <- glm(bite.rate.t ~ Program + Block, family = Gamma(link = 'inverse'), data = individual.trait.20)
br.24 <- glm(bite.rate.t ~ Program + Block, family = Gamma(link = 'inverse'), data = individual.trait.24)
br.28 <- glm(bite.rate.t ~ Program + Block, family = Gamma(link = 'inverse'), data = individual.trait.28)
br.32 <- glm(bite.rate.t ~ Program + Block, family = Gamma(link = 'inverse'), data = individual.trait.32)

summary(br.16)
summary(br.20)
summary(br.24)
summary(br.28)
summary(br.32)

emmeans.br.16 <- emmeans(br.16, pairwise ~ Program)
emmeans.br.20 <- emmeans(br.20, pairwise ~ Program)
emmeans.br.24 <- emmeans(br.24, pairwise ~ Program)
emmeans.br.28 <- emmeans(br.28, pairwise ~ Program)
emmeans.br.32 <- emmeans(br.32, pairwise ~ Program)

p.br.16 <- plot(emmeans.br.16)
p.br.20 <- plot(emmeans.br.20)
p.br.24 <- plot(emmeans.br.24)
p.br.28 <- plot(emmeans.br.28)
p.br.32 <- plot(emmeans.br.32)
grid.arrange(p.br.16 + labs(title = "16C"), p.br.20 + labs(title = "20C"), p.br.24 + labs(title = "24C"),
			 p.br.28 + labs(title = "28C"), p.br.32 + labs(title = "32C"))

#####################
#####################  Lifetime Eggs Analysis
#####################

hist(all.individual.trait$lifetime.eggs)

# Negative binomial fits better than poisson 
fit.nbinom <- fitdist(all.individual.trait$lifetime.eggs, distr = "nbinom", method = "mme") 
plot(fit.nbinom)

fit.pois <- fitdist(all.individual.trait$lifetime.eggs, distr = "pois", method = "mme") 
plot(fit.pois)

le.1 <- glm.nb(lifetime.eggs ~ Tscale*Program, data = all.individual.trait)
le.2 <- glm.nb(lifetime.eggs ~ I(Tscale^2) + Tscale*Program, data = all.individual.trait)
warnings()
summary(le.1)
summary(le.2)


AIC(le.1)
AIC(le.2)

lifeeggs.temp.model <- glm(lifetime.eggs ~ Tscale*Program, family = poisson(link ='log'), data = all.individual.trait)
lifeeggs.temp.model2 <- glm(lifetime.eggs ~ I(Tscale^2) + Tscale*Program, family = poisson(link ='log'), data = all.individual.trait)
Anova(lifeeggs.temp.model)
summary(lifeeggs.temp.model)
AIC(lifeeggs.temp.model)
AIC(lifeeggs.temp.model2)


le.m2 <- glmer(lifetime.eggs ~ Tscale*Program + (1|Block), family = poisson(link = "log"), data = all.individual.trait)
le.m3 <- glmer(lifetime.eggs ~ I(Tscale^2) + Tscale*Program + (1|Block), family = poisson(link = "log"), data = all.individual.trait)
AIC(le.m2)
AIC(le.m3)

Anova(lf.m2)
summary(le.m2)
AIC(le.m2) - AIC(le.m3)



le.1 <- glmer.nb(lifetime.eggs ~ (1|Block), data = all.individual.trait)
le.2 <- glmer.nb(lifetime.eggs ~ Tscale + (1|Block), data = all.individual.trait)
le.3 <- glmer.nb(lifetime.eggs ~ Tscale + I(Tscale^2) + (1|Block), data = all.individual.trait)
le.4 <- glmer.nb(lifetime.eggs ~ Program + (1|Block), data = all.individual.trait)
le.5 <- glmer.nb(lifetime.eggs ~ Tscale + Program + (1|Block), data = all.individual.trait)
le.6 <- glmer.nb(lifetime.eggs ~ Tscale + Program + I(Tscale^2) + (1|Block), data = all.individual.trait)
le.7 <- glmer.nb(lifetime.eggs ~ Tscale*Program + (1|Block), data = all.individual.trait)
le.8 <- glmer.nb(lifetime.eggs ~ Tscale*Program + I(Tscale^2) + (1|Block), data = all.individual.trait)
le.9 <- glmer.nb(lifetime.eggs ~ Tscale*Program + I(Tscale^2)*Program + (1|Block), data = all.individual.trait)

le.8b <- glm.nb(lifetime.eggs ~ Tscale*Program + I(Tscale^2), data = all.individual.trait)

# Calculate delta AIC and AIC weights - remember to round AIC to 0.1 and AICw to 0.01 when reporting
le.AIC.vec <- c(AIC(le.1), AIC(le.2), AIC(le.3), AIC(le.4), AIC(le.5), AIC(le.6), AIC(le.7), AIC(le.8), AIC(le.9))
le.dAIC.vec <- le.AIC.vec - min(le.AIC.vec)
le.AICw.vec <- exp(-0.5*le.dAIC.vec)/sum(exp(-0.5*le.dAIC.vec))

# Get deviance and residual df for reporting stats according to: https://www.ashander.info/posts/2015/10/model-selection-glms-aic-what-to-report/
le.table <- do.call(rbind, lapply(list(le.1, le.2, le.3, le.4, le.5, le.6, le.7, le.8, le.9), broom::glance))

# Get marginal R^2 value (variance explained by fixed effects)
le.sR2.list <- lapply(list(le.2, le.3, le.4, le.5, le.6, le.7, le.8, le.9), r2beta, method = "nsj", partial = T)

# Examine values for manuscript table
le.dAIC.vec
le.AICw.vec
le.table
le.sR2.list

# Summary of best model
summary(le.8)
summary(le.9)

r2beta(le.8)

summary(le.m2nb)
summary(le.m3nb)

AIC(le.m2nb)
AIC(le.m3nb)





library(pscl)
zero.infl1 <- zeroinfl(lifetime.eggs ~ Tscale*Program, data = all.individual.trait) #runs zero-inflated pois
zero.infl2 = zeroinfl(lifetime.eggs ~  Tscale*Program, data = all.individual.trait, dist = "negbin") #runs zero-inflated nb
zero.infl3 = zeroinfl(lifetime.eggs ~ I(Tscale^2) + Tscale*Program, data = all.individual.trait, dist = "negbin") #runs zero-inflated nb

summary(zero.infl1) #zero inflated seems to be working well
summary(zero.infl2) #zero inflated seems to be working well
summary(zero.infl3) #zero inflated seems to be working well
Anova(zero.infl2)
AIC(zero.infl2) - AIC(zero.infl3)


##################### Post-hoc Pairwise Comparisons

##### Constrasts from global model (averaged across all temperature groups)
emmeans(le.8, pairwise ~ Program)

##### Constrasts with temperature treatments
# Subset data
individual.trait.16 <- subset(all.individual.trait, Treatment == 16)
individual.trait.20 <- subset(all.individual.trait, Treatment == 20)
individual.trait.24 <- subset(all.individual.trait, Treatment == 24)
individual.trait.28 <- subset(all.individual.trait, Treatment == 28)
individual.trait.32 <- subset(all.individual.trait, Treatment == 32)

# Fit new GLMMs for each mean temperature
le.16 <- glm.nb(lifetime.eggs ~ Program + Block, data = individual.trait.16)
le.20 <- glm.nb(lifetime.eggs ~ Program + Block, data = individual.trait.20)
le.24 <- glm.nb(lifetime.eggs ~ Program + Block, data = individual.trait.24)
le.28 <- glm.nb(lifetime.eggs ~ Program + Block, data = individual.trait.28)
le.32 <- glm.nb(lifetime.eggs ~ Program + Block, data = individual.trait.32)

summary(le.16)
summary(le.20)
summary(le.24)
summary(le.28)
summary(le.32)

emmeans.le.16 <- emmeans(le.16, pairwise ~ Program)
emmeans.le.20 <- emmeans(le.20, pairwise ~ Program)
emmeans.le.24 <- emmeans(le.24, pairwise ~ Program)
emmeans.le.28 <- emmeans(le.28, pairwise ~ Program)
emmeans.le.32 <- emmeans(le.32, pairwise ~ Program)

p.le.16 <- plot(emmeans.le.16)
p.le.20 <- plot(emmeans.le.20)
p.le.24 <- plot(emmeans.le.24)
p.le.28 <- plot(emmeans.le.28)
p.le.32 <- plot(emmeans.le.32)
grid.arrange(p.le.16 + labs(title = "16C"), p.le.20 + labs(title = "20C"), p.le.24 + labs(title = "24C"),
			 p.le.28 + labs(title = "28C"), p.le.32 + labs(title = "32C"))


individual.trait.16.0 <- subset(individual.trait.16, Program == "dtr0")
individual.trait.16.9 <- subset(individual.trait.16, Program == "dtr9")
individual.trait.16.12 <- subset(individual.trait.16, Program == "dtr12")

individual.trait.20.0 <- subset(individual.trait.20, Program == "dtr0")
individual.trait.20.9 <- subset(individual.trait.20, Program == "dtr9")
individual.trait.20.12 <- subset(individual.trait.20, Program == "dtr12")

individual.trait.24.0 <- subset(individual.trait.24, Program == "dtr0")
individual.trait.24.9 <- subset(individual.trait.24, Program == "dtr9")
individual.trait.24.12 <- subset(individual.trait.24, Program == "dtr12")

individual.trait.28.0 <- subset(individual.trait.28, Program == "dtr0")
individual.trait.28.9 <- subset(individual.trait.28, Program == "dtr9")
individual.trait.28.12 <- subset(individual.trait.28, Program == "dtr12")

individual.trait.32.0 <- subset(individual.trait.32, Program == "dtr0")
individual.trait.32.9 <- subset(individual.trait.32, Program == "dtr9")
individual.trait.32.12 <- subset(individual.trait.32, Program == "dtr12")

bs.le.16.0 <- Bootstrap.le(individual.trait.16.0, 1000)
bs.le.16.9 <- Bootstrap.le(individual.trait.16.9, 1000)
bs.le.16.12 <- Bootstrap.le(individual.trait.16.12, 1000)

bs.le.20.0 <- Bootstrap.le(individual.trait.20.0, 1000)
bs.le.20.9 <- Bootstrap.le(individual.trait.20.9, 1000)
bs.le.20.12 <- Bootstrap.le(individual.trait.20.12, 1000)

bs.le.24.0 <- Bootstrap.le(individual.trait.24.0, 1000)
bs.le.24.9 <- Bootstrap.le(individual.trait.24.9, 1000)
bs.le.24.12 <- Bootstrap.le(individual.trait.24.12, 1000)

bs.le.28.0 <- Bootstrap.le(individual.trait.28.0, 1000)
bs.le.28.9 <- Bootstrap.le(individual.trait.28.9, 1000)
bs.le.28.12 <- Bootstrap.le(individual.trait.28.12, 1000)

bs.le.32.0 <- Bootstrap.le(individual.trait.32.0, 1000)
bs.le.32.9 <- Bootstrap.le(individual.trait.32.9, 1000)
bs.le.32.12 <- Bootstrap.le(individual.trait.32.12, 1000)

#bs.le.16.0  <- filter(individual.trait.16, Program == "dtr0") %>%
	Bootstrap.le(individual.trait.16, 10)

#c("dtr0", "dtr9", "dtr12")
bs.le.16 <- data.frame(dtr = c(1, 2, 3),
					  mean = c(mean(individual.trait.16.0$lifetime.eggs), mean(individual.trait.16.9$lifetime.eggs), mean(individual.trait.16.12$lifetime.eggs)),
					  lowerCI = c(quantile(bs.le.16.0, 0.025), quantile(bs.le.16.9, 0.025), quantile(bs.le.16.12, 0.025)),
					  upperCI = c(quantile(bs.le.16.0, 0.975), quantile(bs.le.16.9, 0.975), quantile(bs.le.16.12, 0.975)))

bs.le.20 <- data.frame(dtr = c(1, 2, 3),
					   mean = c(mean(individual.trait.20.0$lifetime.eggs), mean(individual.trait.20.9$lifetime.eggs), mean(individual.trait.20.12$lifetime.eggs)),
					   lowerCI = c(quantile(bs.le.20.0, 0.025), quantile(bs.le.20.9, 0.025), quantile(bs.le.20.12, 0.025)),
					   upperCI = c(quantile(bs.le.20.0, 0.975), quantile(bs.le.20.9, 0.975), quantile(bs.le.20.12, 0.975)))

bs.le.24 <- data.frame(dtr = c(1, 2, 3),
					   mean = c(mean(individual.trait.24.0$lifetime.eggs), mean(individual.trait.24.9$lifetime.eggs), mean(individual.trait.24.12$lifetime.eggs)),
					   lowerCI = c(quantile(bs.le.24.0, 0.025), quantile(bs.le.24.9, 0.025), quantile(bs.le.24.12, 0.025)),
					   upperCI = c(quantile(bs.le.24.0, 0.975), quantile(bs.le.24.9, 0.975), quantile(bs.le.24.12, 0.975)))

bs.le.28 <- data.frame(dtr = c(1, 2, 3),
					   mean = c(mean(individual.trait.28.0$lifetime.eggs), mean(individual.trait.28.9$lifetime.eggs), mean(individual.trait.28.12$lifetime.eggs)),
					   lowerCI = c(quantile(bs.le.28.0, 0.025), quantile(bs.le.28.9, 0.025), quantile(bs.le.28.12, 0.025)),
					   upperCI = c(quantile(bs.le.28.0, 0.975), quantile(bs.le.28.9, 0.975), quantile(bs.le.28.12, 0.975)))

bs.le.32 <- data.frame(dtr = c(1, 2, 3),
					   mean = c(mean(individual.trait.32.0$lifetime.eggs), mean(individual.trait.32.9$lifetime.eggs), mean(individual.trait.32.12$lifetime.eggs)),
					   lowerCI = c(quantile(bs.le.32.0, 0.025), quantile(bs.le.32.9, 0.025), quantile(bs.le.32.12, 0.025)),
					   upperCI = c(quantile(bs.le.32.0, 0.975), quantile(bs.le.32.9, 0.975), quantile(bs.le.32.12, 0.975)))

par(mfrow = c(3,2))

plot(mean ~ dtr, data = bs.le.16, ylim = c(50,500), pch = 19)
arrows(bs.le.16$dtr, bs.le.16$lowerCI, bs.le.16$dtr, bs.le.16$upperCI, length = 0, angle = 90, code = 3, lwd = 1)

plot(mean ~ dtr, data = bs.le.20, ylim = c(50,500), pch = 19)
arrows(bs.le.20$dtr, bs.le.20$lowerCI, bs.le.20$dtr, bs.le.20$upperCI, length = 0, angle = 90, code = 3, lwd = 1)

plot(mean ~ dtr, data = bs.le.24, ylim = c(50,500), pch = 19)
arrows(bs.le.24$dtr, bs.le.24$lowerCI, bs.le.24$dtr, bs.le.24$upperCI, length = 0, angle = 90, code = 3, lwd = 1)

plot(mean ~ dtr, data = bs.le.28, ylim = c(50,500), pch = 19)
arrows(bs.le.28$dtr, bs.le.28$lowerCI, bs.le.28$dtr, bs.le.28$upperCI, length = 0, angle = 90, code = 3, lwd = 1)

plot(mean ~ dtr, data = bs.le.32, ylim = c(50,500), pch = 19)
arrows(bs.le.32$dtr, bs.le.32$lowerCI, bs.le.32$dtr, bs.le.32$upperCI, length = 0, angle = 90, code = 3, lwd = 1)


par(mfrow = c(3,2))

plot(mean ~ dtr, data = bs.le.16, ylim = c(90,200), pch = 19, main = "16")
arrows(bs.le.16$dtr, bs.le.16$lowerCI, bs.le.16$dtr, bs.le.16$upperCI, length = 0, angle = 90, code = 3, lwd = 1)

plot(mean ~ dtr, data = bs.le.20, ylim = c(100,400), pch = 19, main = "20")
arrows(bs.le.20$dtr, bs.le.20$lowerCI, bs.le.20$dtr, bs.le.20$upperCI, length = 0, angle = 90, code = 3, lwd = 1)

plot(mean ~ dtr, data = bs.le.24, ylim = c(250,500), pch = 19, main = "24")
arrows(bs.le.24$dtr, bs.le.24$lowerCI, bs.le.24$dtr, bs.le.24$upperCI, length = 0, angle = 90, code = 3, lwd = 1)

plot(mean ~ dtr, data = bs.le.28, ylim = c(100,540), pch = 19, main = "28")
arrows(bs.le.28$dtr, bs.le.28$lowerCI, bs.le.28$dtr, bs.le.28$upperCI, length = 0, angle = 90, code = 3, lwd = 1)

plot(mean ~ dtr, data = bs.le.32, ylim = c(50,350), pch = 19, main = "32")
arrows(bs.le.32$dtr, bs.le.32$lowerCI, bs.le.32$dtr, bs.le.32$upperCI, length = 0, angle = 90, code = 3, lwd = 1)

Bootstrap.le = function(dataset, kIterations) {
	
	# create output vector
	out = numeric(kIterations)
	
	# Loop through kIterations
	for (i in 1:kIterations) {
		
		# Simulate a dataset
		simulated.data.set = SimulateDataSet(dataset)
		
		# calculate and store mean lifetime eggs for the simulated data set
		out[i] = mean(simulated.data.set$lifetime.eggs)
		
	}
	
	# return
	out 
}

# Create simulated data set from a treatment data pool ## Somehow need to only work with a subset of data - get row assignments for treatments
SimulateDataSet = function(dataset) {
	
	# Create a sample of random row numbers
	boot.rows = sample(nrow(dataset), size=nrow(dataset), replace=T)
	
	# Create a data frame and fill it with the row from random sample of row numbers
	simulated.data.set = dataset[boot.rows[1], ]
	
	#Finish filling in the data frame with the rest of the rows from random sample of row numbers
	for (iRow in 2:length(boot.rows)){
		
		simulated.data.set = rbind(simulated.data.set, dataset[boot.rows[iRow], ])
		
	}
	
	simulated.data.set #Return
	
}