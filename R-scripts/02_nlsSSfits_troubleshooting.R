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


# This is a function that will take the fitted model coefficients and plot them with your data
plot.fit = function(t, fun, pars){
	if (fun=="quad") lines(t, quad(t, pars['c'], pars['T0'], pars['Tm'])) else{
		lines(t, ss(t, pars['c'], pars['T0'], pars['Tm']))
		#lines(t, briere(t, pars['c'], pars['T0'], pars['Tm']))
	}
}

T = seq(10, 40, by=0.5)
T.K = seq(283, 313, by=0.5)

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