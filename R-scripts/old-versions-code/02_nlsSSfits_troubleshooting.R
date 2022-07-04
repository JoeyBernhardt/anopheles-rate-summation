library("minpack.lm")


##########
####### Fitting Functions Using Non-Linear Least Squares: nls() function
##########

# Functions that we will fit to the data:
# Briere is an asymmetric, hump-shaped function with three parameters:
# c = rate, T0 = minimum temperature, Tm = maximum temperature
briere<-function(t, c, T0, Tm){
	b=c()
	for (i in 1:length(t)){
		if (t[i]>T0 && t[i]<Tm)  {b[i]<-(c*t[i]*(t[i]-T0)*sqrt(Tm-t[i]))} else {b[i]<-(0)}
	}
	b
}

# Quadratic is a symmetric, hump-shaped function with three parameters:
# c = rate, T0 = minimum temperature, Tm = maximum temperature
quad<-function(t, c, T0, Tm){
	b=c()
	for (i in 1:length(t)){
		if (t[i]>T0 && t[i]<Tm)  {b[i]<-(-c*(t[i]-T0)*(t[i]-Tm))} else {b[i]<-(0)}
	}
	b
}

# Sharp-Schoolfield from Dell et al. 2011 has four parameters:
# c, Topt = optimal temperature, E = rise below Topt, ED = decline above Topt, 
ss<-function(t, c, Topt, ED){
	b=c()
	k <- 8.62*10^-5 # Boltzman constant in eV/Kelvin
	E <- 0.6
	for (i in 1:length(t)){
		b[i] <- c * 10^11 * exp(-E/(k*t[i])) / (1 + exp(-1/(k*t[i]) * (ED - t[i]*(ED/Topt + k*log(E/(ED-E)) ) ) ) ) 
	}
	b
}


# Modified Briere function from Mauricio Cruz-Loya  (in prep), has five parameters:
# c, Topt = optimal temperature, E = rise below Topt, ED = decline above Topt, 
CLBriere<-function(t, gmax, Tmin, Tmax, alpha, s){
	b=c()
	for (i in 1:length(t)){
		b[i] <- gmax * ( ((t[i]-Tmin)/alpha)^alpha * ((Tmax-t[i])/(1-alpha))^(1-alpha) * (1/(Tmax-Tmin)) )^s
	}
	b
}


# This is a function that will take the fitted model coefficients and plot them with your data
plot.fit = function(t, fun, pars){
	if (fun=="quad") lines(t, quad(t, pars['c'], pars['T0'], pars['Tm'])) else{
		#lines(t, ss(t, pars['c'], pars['T0'], pars['Tm']))
		lines(t, briere(t, pars['c'], pars['T0'], pars['Tm']))
	}
}

T = seq(5, 45, by=0.5)
T.K = seq(278.15, 318.15, by=0.5)

cdata$Temp.K <- cdata$Treatment + 273.15
alldata.sum.2.c <- subset(alldata.sum.2, DTR == "constant")

fit.lf.q = nlsLM(lifespan ~ quad(Treatment, c, T0, Tm), data = cdata,
				start=list(c = 0.1, T0 = 10, Tm = 40))
coef.lf.q = coef(fit.lf.q)
par(mfrow = c(1,1))
plot(lifespan.mean ~ Treatment, data = alldata.sum.2, ylab = "Lifespan", xlab = "Temperature", ylim = c(0,50))
plot.fit(T, "quad", coef.lf.q)

# Try to fit E
fit.lf.ss = nlsLM(lifespan ~ ss(Temp.K, c, Topt, E, ED), data = cdata,
				 start=list(c = 12, Topt = 293, E = 0.6, ED = 2))

# Set E = 0.6
fit.lf.ss = nlsLM(lifespan ~ ss(Temp.K, c, Topt, ED), data = cdata,
				  start=list(c = 12, Topt = 293, ED = 2))
coef.lf.ss = coef(fit.lf.ss)
par(mfrow = c(1,1))
plot(lifespan.mean ~ Treatment, data = alldata.sum.2, ylab = "Lifespan", xlab = "Temperature", ylim = c(0,50))
plot.fit(T, "quad", coef.lf.q)

c <- 12
Topt <- 20 + 273
E <- 0.6
ED <- 2
k <- 8.62*10^-5

c <- coef.lf.ss[1]
Topt <- coef.lf.ss[2]
E <- coef.lf.ss[3]
ED <- coef.lf.ss[4]
k <- 8.62*10^-5

c <- coef.lf.ss[1]
Topt <- coef.lf.ss[2]
E <- 0.6
ED <- coef.lf.ss[3]
k <- 8.62*10^-5

ss.points <- c * 10^11 * exp(-E/(k*T.K)) / (1 + exp(-1/(k*T.K) * (ED - T.K*(ED/Topt + k*log(E/(ED-E)) ) ) ) ) 
plot(ss.points ~ T, ylim = c(0,50))
points(lifespan.mean ~ Treatment, data = alldata.sum.2.c, col = "red")

ss.points <- c * 10^11 * exp(-E/(k*T.K)) 
ss.points <- (1 + exp(-1/(k*T.K) * (ED - T.K*(ED/Topt + k*log(E/(ED-E)) ) ) ) ) 
ss.points <- exp(-1/(k*T.K) * (ED - T.K*(ED/Topt + k*log(E/(ED-E)) ) ) )
ss.points <- -1/(k*T.K)
ss.points <- (ED - T.K*(ED/Topt + k*log(E/(ED-E)) ) )
ss.points <- -1/(k*T.K) * (ED - T.K*(ED/Topt + k*log(E/(ED-E)) ) )
ss.points <- k*log(E/(ED-E))
ss.points <- E/(ED-E)
ss.points

#### Amarasekare and Savage 2012

EA <- 6387 #0.6/1.987
Tref <- 293.15
DTref <- 0.2
EH <- 28149 # 1.2/1.987
TH = 312

AS2012.points <- (T.K * Tref * DTref * exp(EA * ((1/Tref)-(1/T.K)) ) ) / (1 + exp(EH * ((1/TH)-(1/T.K)) ) )
plot(AS2012.points ~ T)																  


#### Amarasekare and Johnson 2017
Tref <- 290
DRref <- 0.016
EA <- 6387
EH <- 28149
EL <- -36734
TH <- 312
TL <- 290

AJ2017.points <- (T.K / Tref * DRref * exp(EA * ((1/Tref)-(1/T.K)) ) ) / (1 + exp(EH * ((1/TH)-(1/T.K)) ) )
plot(AJ2017.points ~ T)	

AJ2017.points <- (T.K / Tref * DRref * exp(EA * ((1/Tref)-(1/T.K)) ) ) / (1 + exp(EH * ((1/TH)-(1/T.K)) ) + exp(EL * ((1/TL)-(1/T.K)) ) )
plot(AJ2017.points ~ T)	

AJ2017.points <- (T.K / Tref * DRref * exp(EA * ((1/Tref)-(1/T.K)) ) )
plot(AJ2017.points ~ T)																			  


#### Cruz-Loya's Modified Briere (in prep)
gmax <- 1 # maximum value at the Topt
Tmin <- 10
Tmax <- 30
alpha <- 0.75 # skew parameter 0 < alpha < 1, 0.5 = symmetrical
s <- 3.2 # peakiness parameter, must be +

CLBriere.points <- gmax * ( ((T-Tmin)/alpha)^alpha * ((Tmax-T)/(1-alpha))^(1-alpha) * (1/(Tmax-Tmin)) )^s
plot(CLBriere.points ~ T, xlim = c(10,30), ylim = c(0,1.2))

# Try to fit CL Briere to Cpip MDR data
fit.MDR.CLBriere = nlsLM(1/trait ~ CLBriere(T.C, gmax, Tmin, Tmax, alpha, s), data = data.MDR.Cpip,
				  start=list(gmax = .1, Tmin = 5, Tmax = 45, alpha = 0.5, s = 5))
coef.Cpip.MDR.CLBriere = coef(fit.MDR.CLBriere)
gmax <- coef.Cpip.MDR.CLBriere[1]
Tmin <- coef.Cpip.MDR.CLBriere[2]
Tmax <- coef.Cpip.MDR.CLBriere[3]
alpha <- coef.Cpip.MDR.CLBriere[4]
s <- coef.Cpip.MDR.CLBriere[5]

fit.MDR.CLBriere.MLE = mle2(1/trait ~ dnorm(mean = gmax * ( ((T.C-Tmin)/alpha)^alpha * ((Tmax-T.C)/(1-alpha))^(1-alpha) * (1/(Tmax-Tmin)) )^s, sd = sd), 
			   start=list(gmax = .1, Tmin = 1, Tmax = 40, alpha = 0.75, s = 3, sd = 1), data = data.MDR.Cpip)
summary(fit.MDR.CLBriere.MLE)

CLBriere.points <- gmax * ( ((T-Tmin)/alpha)^alpha * ((Tmax-T)/(1-alpha))^(1-alpha) * (1/(Tmax-Tmin)) )^s

plot(1/trait ~ T.C, data = data.MDR.Cpip, col = "grey", xlim = c(0,50))
lines(CLBriere.points ~ T)

##### Bayesian fits of Briere / CLBriere

## parameters to estimate
parameters <- c("cf.gmax", "cf.Tmin", "cf.Tmax", "cf.alpha", "cf.s", "cf.sigma")

# Initial values for the parameters
inits<-function(){list(
	cf.gmax = 0.1,
	cf.Tmin = 2,
	cf.Tmax = 38,
	cf.alpha = 0.8,
	cf.s = 4,
	cf.sigma = .1
)}


##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma")

##### inits Function
inits<-function(){list(
	cf.q = 0.01,
	cf.Tm = 35,
	cf.T0 = 5,
	cf.sigma = rlnorm(1))}




# MCMC Settings: number of posterior dist elements = [(ni - nb) / nt ] * nc
ni <- 25000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains

# Temperature sequence for derived quantity calculations
Temp.xs <- seq(0, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


### Fitting the trait thermal response; Pull out data columns as vectors
data <- data.MDR.Cpip[-55,] # only works if I remove the 40 degree datapoint
data$T <- data$T.C

trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T
plot(trait ~ temp)

# Bundle all data in a list for JAGS
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)


CpipMDR.CLBriere <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
				   model.file="briere.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
				   n.iter=ni, DIC=T, working.directory=getwd())	
CpipMDR.CLBriere$BUGSoutput$summary

cf.T0 <- CpipMDR.CLBriere$BUGSoutput$summary[1,1]
cf.Tm <- CpipMDR.CLBriere$BUGSoutput$summary[2,1]
cf.q <- CpipMDR.CLBriere$BUGSoutput$summary[3,1]
cfs <- c(cf.q, cf.T0, cf.Tm)
briere.points <- briere(T, cf.q, cf.T0, cf.Tm)

plot(1/trait ~ T.C, data = data.MDR.Cpip, col = "grey", xlim = c(0,50))
lines(briere.points ~ T, col = "red")

CpipMDR.CLBriere <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
						 model.file="CLbriere.txt", n.thin=nt, n.chains=nc, n.burnin=nb, 
						 n.iter=ni, DIC=T, working.directory=getwd())	
CpipMDR.CLBriere$BUGSoutput$summary

Tmax <- CpipMDR.CLBriere$BUGSoutput$summary[1,1]
Tmin <- CpipMDR.CLBriere$BUGSoutput$summary[2,1]
alpha <- CpipMDR.CLBriere$BUGSoutput$summary[3,1]
gmax <- CpipMDR.CLBriere$BUGSoutput$summary[4,1]
s <- CpipMDR.CLBriere$BUGSoutput$summary[5,1]
CLBriere.points <- gmax * ( ((T-Tmin)/alpha)^alpha * ((Tmax-T)/(1-alpha))^(1-alpha) * (1/(Tmax-Tmin)) )^s

plot(1/trait ~ T.C, data = data.MDR.Cpip, col = "grey", xlim = c(0,50))
lines(CLBriere.points ~ T, col = "red")


trait.mu[i] <- cf.gmax * ( ((temp[i]-cf.Tmin)/cf.alpha)^cf.alpha * ((cf.Tmax-temp[i])/(1-cf.alpha))^(1-cf.alpha) * (1/(cf.Tmax-cf.Tmin)) )^cf.s
