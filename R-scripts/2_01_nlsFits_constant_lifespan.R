###Code written by Kerri Miazgowicz July 30 2020
################################# Contents
# 1. Set-up + data vis 
# 2. NLS fits
################################# 1. Set-up + data vis

###### Load packages
library(plyr)
library(dplyr)
library(tidyverse)

#Set working directory
mainDir = "~/Dropbox/Research Savage Lab/Rate Summation Round 2"
setwd(mainDir)

######  Load raw trait data
alldata <- data.frame(read.csv("fluc.program.trait.means.csv", colClasses = c("DTR" ="character")))

######  Load temp data
tempdata <- data.frame(read.csv("temp.programs.continuous.csv"))


#change default plot specifications
par(mar=c(1,1,1,1))
dev.off()

#example for how to seperate out the data for each DTR
constant.data <- subset(alldata, DTR == "constant")
dtr9.data<- subset(alldata, DTR == "+/- 4.5")
dtr12.data<- subset(alldata, DTR == "+/- 6")




###################################################################

#Can I pre-define each function and call it in the nls fitting - YES
#Define each function of interest
QUAD <- function(cf.q,cf.T0,cf.Tm, Treatment){
  return(-1 * cf.q * (Treatment - cf.T0) * (Treatment - cf.Tm))
}
BRIERE <- function(cf.q, cf.T0,cf.Tm,Treatment){
  return( cf.q * Treatment * (Treatment - cf.T0) * sqrt((cf.Tm - Treatment)*(cf.Tm > Treatment)))
}
NORBERG <- function(cf.a, cf.b, cf.w, cf.z, Treatment){
  return( cf.a * exp(cf.b * Treatment) * (1-((Treatment-cf.z)/(cf.w/2))^2) )
}
MOD_SS <- function(cf.c, cf.Topt, cf.ED,TempK){ #Needs to be fit with temperature in Kelvin
  return(cf.c * 10^11 * exp(-0.6/(8.62*10^-5*TempK)) / (1 + exp(-1/(8.62*10^-5*TempK) * (cf.ED - TempK*(cf.ED/cf.Topt + 8.62*10^-5*log(0.6/(cf.ED-0.6)) ) ) ) ) )
}
CLBRIERE <- function(cf.gmax, cf.T0, cf.Tm, cf.alpha, cf.s, Treatment){
  return(cf.gmax * ( ((Treatment-cf.T0)/cf.alpha)^cf.alpha * ((cf.Tm-Treatment)/(1-cf.alpha))^(1-cf.alpha) * (1/(cf.Tm-cf.T0)) )^cf.s)
}


# clbriere<-function(t, gmax, Tmin, Tmax, alpha, s){
#   b=c()
#   for (i in 1:length(t)){
#     b[i] <- gmax * ( ((t[i]-Tmin)/alpha)^alpha * ((Tmax-t[i])/(1-alpha))^(1-alpha) * (1/(Tmax-Tmin)) )^s
#   }
#   b
# }


#call nls function over each dataset

#add in error around the parameter estimates
lf.quad <- nls(lifespan ~ QUAD(cf.q, cf.T0, cf.Tm, Treatment),
               data= constant.data,
               start = list(cf.q=0.15, cf.T0=4.7, cf.Tm=36.5))
plot(constant.data$Treatment, constant.data$lifespan, xlim= c(0,50), ylim = c(0,45))
lines(seq(0,50,0.1),predict(lf.quad,data.frame(Treatment= seq(0,50,0.1))))

lf.briere <- nls(lifespan ~ BRIERE(cf.q, cf.T0,cf.Tm,Treatment),
               data= constant.data,
               start = list(cf.q=0.03, cf.T0=3.8, cf.Tm=29.6))
plot(constant.data$Treatment, constant.data$lifespan, xlim= c(0,50), ylim = c(0,45))
lines(seq(0,50,0.1),predict(lf.briere,data.frame(Treatment= seq(0,50,0.1))))

lf.norberg <- nls(lifespan ~ NORBERG(cf.a, cf.b, cf.w, cf.z, Treatment),
                  data = constant.data,
                  start = list(cf.a=17, cf.b=0.03, cf.w=28, cf.z=21))
plot(constant.data$Treatment, constant.data$lifespan, xlim= c(0,50), ylim = c(0,45))
lines(seq(0,50,0.1),predict(lf.norberg,data.frame(Treatment= seq(0,50,0.1))))

#Fit this in temp Kelvin
lf.modSS <- nls(lifespan ~ MOD_SS(cf.c, cf.Topt, cf.ED, TempK),
                data = constant.data,
                start = list(cf.c=10 , cf.Topt=294, cf.ED = 2))
plot(constant.data$Treatment, constant.data$lifespan, xlim= c(0,50), ylim = c(0,45))
lines(seq(0,50,0.1),predict(lf.modSS,data.frame(TempK= seq(273.15,273.15+50,0.1))))

###FIT ANY ADDITIONAL FUNCTIONS HERE
lf.clbriere <- nls(lifespan ~ CLBRIERE(cf.gmax, cf.T0, cf.Tm, cf.alpha, cf.s, Treatment),
                 data= constant.data,
                 start = list(cf.gmax=40, cf.T0=5, cf.Tm=42, cf.alpha = 0.45, cf.s = 4))
plot(constant.data$Treatment, constant.data$lifespan, xlim= c(0,50), ylim = c(0,45))
lines(seq(0,50,0.1),predict(lf.briere,data.frame(Treatment= seq(0,50,0.1))))

T <- seq(0,50,0.1)
gmax <- 40
Tmin <- 5
Tmax <- 42
alpha <- 0.45
s <- 4
CLBriere.points <- gmax * ( ((T-Tmin)/alpha)^alpha * ((Tmax-T)/(1-alpha))^(1-alpha) * (1/(Tmax-Tmin)) )^s
lines(CLBriere.points ~ T)


# bite rate
br.briere <- nls(bite.rate ~ BRIERE(cf.q, cf.T0,cf.Tm,Treatment),
                 data= constant.data,
                 start = list(cf.q=0.03, cf.T0=3.8, cf.Tm=29.6))
plot(constant.data$Treatment, constant.data$bite.rate, xlim= c(0,50), ylim = c(0,45))
lines(seq(0,50,0.1),predict(br.briere,data.frame(Treatment= seq(0,50,0.1))))

# lifetime fecundity
f.quad <- nls(lifetime.eggs ~ QUAD(cf.q, cf.T0, cf.Tm, Treatment),
               data= constant.data,
               start = list(cf.q=0.15, cf.T0=4.7, cf.Tm=36.5))
plot(constant.data$Treatment, constant.data$lifetime.eggs, xlim= c(0,50), ylim = c(0,545))
lines(seq(0,50,0.1),predict(f.quad,data.frame(Treatment= seq(0,50,0.1))))


#extracting model fit parameters
lf.quad$m$getPars() #cf.q- 0.1186192, cf.T0- 1.5491849, cf.Tm- 37.571815

#determining the best fitting model
#Joey/Marta - How do we extract the AIC values from these fits to determine the 'best' functional form for each trait | dtr
#Joey- How do I get the parameter error interval???

summary(lf.norberg)
AIC(lf.norberg)
AIC(lf.briere)
AIC(lf.quad)
AIC(lf.modSS)


cf.norberg <- coef(lf.norberg)
cf.norberg

cf.quad <- coef(lf.quad)
cf.quad

cf.briere <- coef(br.briere)
cf.briere

f.quad <- coef(f.quad)
f.quad

# Get vector of temperatures
tempdata$X16.dtr9
tempdata$X20.dtr9
tempdata$X24.dtr9
tempdata$X28.dtr9
tempdata$X32.dtr9

tempdata$X16.dtr12
tempdata$X20.dtr12
tempdata$X24.dtr12
tempdata$X28.dtr12
tempdata$X32.dtr12



# Calculate estimated lifespan for given temperatures for norberg
rs.pred.9 <- data.frame(meantemp = c(16, 20, 24, 28, 32), rs.pred = numeric(5))
rs.pred.12 <- data.frame(meantemp = c(16, 20, 24, 28, 32), rs.pred = numeric(5))

rs.pred.9$rs.pred[1] <- mean(NORBERG(cf.norberg[[1]], cf.norberg[[2]], cf.norberg[[3]], cf.norberg[[4]], tempdata$X16.dtr9))
rs.pred.9$rs.pred[2] <- mean(NORBERG(cf.norberg[[1]], cf.norberg[[2]], cf.norberg[[3]], cf.norberg[[4]], tempdata$X20.dtr9))
rs.pred.9$rs.pred[3] <- mean(NORBERG(cf.norberg[[1]], cf.norberg[[2]], cf.norberg[[3]], cf.norberg[[4]], tempdata$X24.dtr9))
rs.pred.9$rs.pred[4] <- mean(NORBERG(cf.norberg[[1]], cf.norberg[[2]], cf.norberg[[3]], cf.norberg[[4]], tempdata$X28.dtr9))
rs.pred.9$rs.pred[5] <- mean(NORBERG(cf.norberg[[1]], cf.norberg[[2]], cf.norberg[[3]], cf.norberg[[4]], tempdata$X32.dtr9))

rs.pred.12$rs.pred[1] <- mean(NORBERG(cf.norberg[[1]], cf.norberg[[2]], cf.norberg[[3]], cf.norberg[[4]], tempdata$X16.dtr12))
rs.pred.12$rs.pred[2] <- mean(NORBERG(cf.norberg[[1]], cf.norberg[[2]], cf.norberg[[3]], cf.norberg[[4]], tempdata$X20.dtr12))
rs.pred.12$rs.pred[3] <- mean(NORBERG(cf.norberg[[1]], cf.norberg[[2]], cf.norberg[[3]], cf.norberg[[4]], tempdata$X24.dtr12))
rs.pred.12$rs.pred[4] <- mean(NORBERG(cf.norberg[[1]], cf.norberg[[2]], cf.norberg[[3]], cf.norberg[[4]], tempdata$X28.dtr12))
rs.pred.12$rs.pred[5] <- mean(NORBERG(cf.norberg[[1]], cf.norberg[[2]], cf.norberg[[3]], cf.norberg[[4]], tempdata$X32.dtr12))

# Calculate estimated lifespan for given temperatures for quad
rs.pred.9q <- data.frame(meantemp = c(16, 20, 24, 28, 32), rs.pred = numeric(5))
rs.pred.12q <- data.frame(meantemp = c(16, 20, 24, 28, 32), rs.pred = numeric(5))

rs.pred.9q$rs.pred[1] <- mean(QUAD(cf.quad[[1]], cf.quad[[2]], cf.quad[[3]], tempdata$X16.dtr9))
rs.pred.9q$rs.pred[2] <- mean(QUAD(cf.quad[[1]], cf.quad[[2]], cf.quad[[3]], tempdata$X20.dtr9))
rs.pred.9q$rs.pred[3] <- mean(QUAD(cf.quad[[1]], cf.quad[[2]], cf.quad[[3]], tempdata$X24.dtr9))
rs.pred.9q$rs.pred[4] <- mean(QUAD(cf.quad[[1]], cf.quad[[2]], cf.quad[[3]], tempdata$X28.dtr9))
rs.pred.9q$rs.pred[5] <- mean(QUAD(cf.quad[[1]], cf.quad[[2]], cf.quad[[3]], tempdata$X32.dtr9))

rs.pred.12q$rs.pred[1] <- mean(QUAD(cf.quad[[1]], cf.quad[[2]], cf.quad[[3]], tempdata$X16.dtr12))
rs.pred.12q$rs.pred[2] <- mean(QUAD(cf.quad[[1]], cf.quad[[2]], cf.quad[[3]], tempdata$X20.dtr12))
rs.pred.12q$rs.pred[3] <- mean(QUAD(cf.quad[[1]], cf.quad[[2]], cf.quad[[3]], tempdata$X24.dtr12))
rs.pred.12q$rs.pred[4] <- mean(QUAD(cf.quad[[1]], cf.quad[[2]], cf.quad[[3]], tempdata$X28.dtr12))
rs.pred.12q$rs.pred[5] <- mean(QUAD(cf.quad[[1]], cf.quad[[2]], cf.quad[[3]], tempdata$X32.dtr12))


# Calculate estimated bite for given temperatures for briere
br.pred.9 <- data.frame(meantemp = c(16, 20, 24, 28, 32), rs.pred = numeric(5))
br.pred.12 <- data.frame(meantemp = c(16, 20, 24, 28, 32), rs.pred = numeric(5))

br.pred.9$rs.pred[1] <- mean(BRIERE(cf.briere[[1]], cf.briere[[2]], cf.briere[[3]], tempdata$X16.dtr9))
br.pred.9$rs.pred[2] <- mean(BRIERE(cf.briere[[1]], cf.briere[[2]], cf.briere[[3]], tempdata$X20.dtr9))
br.pred.9$rs.pred[3] <- mean(BRIERE(cf.briere[[1]], cf.briere[[2]], cf.briere[[3]], tempdata$X24.dtr9))
br.pred.9$rs.pred[4] <- mean(BRIERE(cf.briere[[1]], cf.briere[[2]], cf.briere[[3]], tempdata$X28.dtr9))
br.pred.9$rs.pred[5] <- mean(BRIERE(cf.briere[[1]], cf.briere[[2]], cf.briere[[3]], tempdata$X32.dtr9))

br.pred.12$rs.pred[1] <- mean(BRIERE(cf.briere[[1]], cf.briere[[2]], cf.briere[[3]], tempdata$X16.dtr12))
br.pred.12$rs.pred[2] <- mean(BRIERE(cf.briere[[1]], cf.briere[[2]], cf.briere[[3]], tempdata$X20.dtr12))
br.pred.12$rs.pred[3] <- mean(BRIERE(cf.briere[[1]], cf.briere[[2]], cf.briere[[3]], tempdata$X24.dtr12))
br.pred.12$rs.pred[4] <- mean(BRIERE(cf.briere[[1]], cf.briere[[2]], cf.briere[[3]], tempdata$X28.dtr12))
br.pred.12$rs.pred[5] <- mean(BRIERE(cf.briere[[1]], cf.briere[[2]], cf.briere[[3]], tempdata$X32.dtr12))


# Calculate estimated fecundity for given temperatures for quad
f.pred.9 <- data.frame(meantemp = c(16, 20, 24, 28, 32), rs.pred = numeric(5))
f.pred.12 <- data.frame(meantemp = c(16, 20, 24, 28, 32), rs.pred = numeric(5))

f.pred.9$rs.pred[1] <- mean(QUAD(f.quad[[1]], f.quad[[2]], f.quad[[3]], tempdata$X16.dtr9))
f.pred.9$rs.pred[2] <- mean(QUAD(f.quad[[1]], f.quad[[2]], f.quad[[3]], tempdata$X20.dtr9))
f.pred.9$rs.pred[3] <- mean(QUAD(f.quad[[1]], f.quad[[2]], f.quad[[3]], tempdata$X24.dtr9))
f.pred.9$rs.pred[4] <- mean(QUAD(f.quad[[1]], f.quad[[2]], f.quad[[3]], tempdata$X28.dtr9))
f.pred.9$rs.pred[5] <- mean(QUAD(f.quad[[1]], f.quad[[2]], f.quad[[3]], tempdata$X32.dtr9))

f.pred.12$rs.pred[1] <- mean(QUAD(f.quad[[1]], f.quad[[2]], f.quad[[3]], tempdata$X16.dtr12))
f.pred.12$rs.pred[2] <- mean(QUAD(f.quad[[1]], f.quad[[2]], f.quad[[3]], tempdata$X20.dtr12))
f.pred.12$rs.pred[3] <- mean(QUAD(f.quad[[1]], f.quad[[2]], f.quad[[3]], tempdata$X24.dtr12))
f.pred.12$rs.pred[4] <- mean(QUAD(f.quad[[1]], f.quad[[2]], f.quad[[3]], tempdata$X28.dtr12))
f.pred.12$rs.pred[5] <- mean(QUAD(f.quad[[1]], f.quad[[2]], f.quad[[3]], tempdata$X32.dtr12))




par(mfrow = c(2,3))
plot(constant.data$Treatment, constant.data$lifespan, xlim= c(0,50), ylim = c(0,45), ylab = "Lifespan (days)", xlab = "(Temperature (C)")
lines(seq(0,50,0.1),predict(lf.quad,data.frame(Treatment= seq(0,50,0.1))))
text(x = 40, y = 40, labels = "Quad AIC = 80")

plot(constant.data$Treatment, constant.data$lifespan, xlim= c(0,50), ylim = c(0,45))
lines(seq(0,50,0.1),predict(lf.briere,data.frame(Treatment= seq(0,50,0.1))))
text(x = 40, y = 40, labels = "Briere AIC = 98")

plot(constant.data$Treatment, constant.data$lifespan, xlim= c(0,50), ylim = c(0,45))
lines(seq(0,50,0.1),predict(lf.norberg,data.frame(Treatment= seq(0,50,0.1))))
text(x = 40, y = 40, labels = "Norberg AIC = 75")

plot(constant.data$Treatment, constant.data$lifespan, xlim= c(0,50), ylim = c(0,45))
lines(seq(0,50,0.1),predict(lf.modSS,data.frame(TempK= seq(273.15,273.15+50,0.1))))
text(x = 40, y = 40, labels = "SS AIC = 71")

plot(constant.data$Treatment, constant.data$lifespan, xlim= c(0,50), ylim = c(0,45))
lines(CLBriere.points ~ T)
text(x = 40, y = 40, labels = "CL Briere AIC = ???")




plot(constant.data$Treatment, constant.data$lifespan, xlim= c(0,50), ylim = c(0,45))
points(dtr9.data$Treatment, dtr9.data$lifespan, col = "blue")
points(dtr12.data$Treatment, dtr12.data$lifespan, col = "dodgerblue")

plot(constant.data$Treatment, constant.data$bite.rate, xlim= c(0,50), ylim = c(0,0.6))
points(dtr9.data$Treatment, dtr9.data$bite.rate, col = "blue")
points(dtr12.data$Treatment, dtr12.data$bite.rate, col = "dodgerblue")

plot(constant.data$Treatment, constant.data$lifetime.eggs, xlim= c(0,50), ylim = c(0,545))
points(dtr9.data$Treatment, dtr9.data$lifetime.eggs, col = "blue")
points(dtr12.data$Treatment, dtr12.data$lifetime.eggs, col = "dodgerblue")

#lf Norberg
plot(constant.data$Treatment, constant.data$lifespan, xlim= c(0,50), ylim = c(0,45), col = "white")
lines(seq(0,50,0.1),predict(lf.norberg,data.frame(Treatment= seq(0,50,0.1))), lwd = 2)
points(dtr9.data$Treatment, dtr9.data$lifespan, col = "blue")
points(dtr12.data$Treatment, dtr12.data$lifespan, col = "dodgerblue")

points(rs.pred.9$rs.pred ~ rs.pred.9$meantemp, pch = 19, col = "blue")
points(rs.pred.12$rs.pred ~ rs.pred.12$meantemp, pch = 19, col = "dodgerblue")

#lf Quad
plot(constant.data$Treatment, constant.data$lifespan, xlim= c(0,50), ylim = c(0,45), col = "white", ylab = "Lifespan (days)", xlab = "(Temperature (C)")
lines(seq(0,50,0.1),predict(lf.quad,data.frame(Treatment= seq(0,50,0.1))), lwd = 2)
points(dtr9.data$Treatment, dtr9.data$lifespan, col = "blue")
points(dtr12.data$Treatment, dtr12.data$lifespan, col = "dodgerblue")

points(rs.pred.9q$rs.pred ~ rs.pred.9q$meantemp, pch = 19, col = "blue")
points(rs.pred.12q$rs.pred ~ rs.pred.12q$meantemp, pch = 19, col = "dodgerblue")

#br briere
plot(constant.data$Treatment, constant.data$bite.rate, xlim= c(0,50), ylim = c(0,0.6), col = "black", ylab = "Biting rate (1/day)", xlab = "(Temperature (C)")
plot(constant.data$Treatment, constant.data$bite.rate, xlim= c(0,50), ylim = c(0,0.6), col = "white", ylab = "Biting rate (1/day)", xlab = "(Temperature (C)")
lines(seq(0,50,0.1),predict(br.briere,data.frame(Treatment= seq(0,50,0.1))), lwd = 2)
points(dtr9.data$Treatment, dtr9.data$bite.rate, col = "blue")
points(dtr12.data$Treatment, dtr12.data$bite.rate, col = "dodgerblue")

points(br.pred.9$rs.pred ~  br.pred.9$meantemp, pch = 19, col = "blue")
points(br.pred.12$rs.pred ~ br.pred.12$meantemp, pch = 19, col = "dodgerblue")

#lifetime eggs Quad
plot(constant.data$Treatment, constant.data$lifetime.eggs, xlim= c(0,50), ylim = c(0,545), col = "black", ylab = "Lifetime eggs", xlab = "(Temperature (C)")
plot(constant.data$Treatment, constant.data$lifetime.eggs, xlim= c(0,50), ylim = c(0,545), col = "white", ylab = "Lifetime eggs", xlab = "(Temperature (C)")
lines(seq(0,50,0.1),predict(f.quad,data.frame(Treatment= seq(0,50,0.1))), lwd = 2)
points(dtr9.data$Treatment, dtr9.data$lifetime.eggs, col = "blue")
points(dtr12.data$Treatment, dtr12.data$lifetime.eggs, col = "dodgerblue")

points(f.pred.9$rs.pred ~ f.pred.9$meantemp, pch = 19, col = "blue")
points(f.pred.12$rs.pred ~ f.pred.12$meantemp, pch = 19, col = "dodgerblue")
