## Marta Shocket, Stanford University
## Updated by Erin Mordecai, August 10, 2018
## Updated by Kerri Miazgowicz, 2019
## Updated by KM April 13, 2020
##
## Purpose: Perform sensitivity analysis to determine which mosquito and parasite traits drive the
##          thermal response of R0.
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) Deriviative and helper functions
##           3) Calculate trait means across temp
##           4) Sensitivity analysis #1 - partial derivitatives
##           4) Sensitivity analysis #2 - holding single parameters constant
##           5) Uncertainty analysis


##########
###### 1. Set up workspace, load packages, get data, etc.
##########

# Set working directory
mainDir = "C:/Users/Kerri/Desktop/Chapter1 Submission"
setwd(mainDir)

##Load prior fits with informative priors where possible

############ Load traits fits
# Fits from Miazgowicz and Shapiro data - Informative priors
load("saved posteriors inf/lifespan_quadratic_inf.Rdata")
lifespan.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
#for senstivity analysis #1 I need these files
lifespan.out <- model.out$BUGSoutput$sims.list


load("saved posteriors inf/bite.rate_briere_inf.Rdata")
bite.rate.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
#for senstivity analysis #1 I need these files
bite.rate.out <- model.out$BUGSoutput$sims.list


#Note: due to note having an appropiate trait to fit informative priors to we use the uniformative fits
load("saved posteriors/lifeeggs_briere_withT_means_uniform.Rdata")
lifetime.eggs.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
lifetime.eggs.out <-  model.out$BUGSoutput$sims.list


#Note: due to the high amount of uncertainty on the inf fit; we are using the uni for this trait
load("saved posteriors/est.EFD_briere_withT_means_uniform.Rdata")
EFD.1stGC.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
EFD.1stGC.out <-  model.out$BUGSoutput$sims.list

load("saved posteriors inf/est.bite_briere_inf.Rdata")
est.bite.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
est.bite.out <-  model.out$BUGSoutput$sims.list

load("saved posteriors inf/mu_exponential_quadratic_inf.Rdata")
mu.exp.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred  # actually fits for 1/mu
mu.exp.out <-  model.out$BUGSoutput$sims.list

load("saved posteriors inf/shapiro_PDR_briere_inf.Rdata")
PDR.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred  # actually fits for 1/EIP50
PDR.out <-  model.out$BUGSoutput$sims.list

load("saved posteriors inf/shapiro_bc_quadratic_inf.Rdata")
bc.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
bc.out <-  model.out$BUGSoutput$sims.list

#no informative priors can be added to this one
load("saved posteriors/gamma_quad_withT_all_uniform.Rdata")
gamma.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
gamma.out <-  model.out$BUGSoutput$sims.list

#Note: to prevent overlap in preds variable names I add a "K" after pea and MDR
load("saved posteriors inf/peaK_quadratic_inf.Rdata")
pEA.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
pEA.out <-  model.out$BUGSoutput$sims.list

load("saved posteriors inf/MDRK_briere_inf.Rdata")
MDR.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
MDR.out <-  model.out$BUGSoutput$sims.list


# Temperature levels and # MCMC steps
load("saved posteriors/temps.Rdata")
N.Temp.xs <- length(Temp.xs)
nMCMC <- 7500

# Creating a small constant to keep denominators from being zero
ec <- 0.000001

##########
###### 2. Derivative and helper Functions
##########

# Arguments:
#   t           vector of temp gradient
#   q, Tm, T0   thermal response function coefficient posterior distributions

# Function for derivative of Briere thermal response
d_briere = function(t, T0, Tm, q) {
  
  b <- c()
  
  for (i in 1:length(t)) {
    if (t[i]>T0 && t[i]<Tm) {b[i] <- (q*(-5*(t[i]^2) + 3*t[i]*T0 + 4*t[i]*Tm - 2*T0*Tm)/(2*sqrt(Tm-t[i])))}
    else {b[i] <- 0}
  }
  
  b # return output
  
}

# Function for derivative of quadratic thermal response
d_quad = function(t, T0, Tm, q){
  
  b <- c()
  
  for (i in 1:length(t)){
    if (t[i]>T0 && t[i]<Tm) {b[i] <- -1*q*(2*t[i] - T0 - Tm)}
    else {b[i] <- 0}
  }
  
  b # return output

}

### Two R0 functions:
# 1.An.stephensi estimated model: M depends on eggs per female per day
R0.est = function(a, bc, lf, PDR, EFD, pEA, MDR){
  M = EFD * pEA * MDR * (lf+ec.lf)^2
  (a^2 * bc * exp(-(1/(lf+ec.lf))*(1/(PDR+ec.pdr))) * M * (lf+ec.lf) )^0.5
}

# 2.An. stephensi lifetime model: M depends on lifetime fecundity
R0.full = function(a, bc, lf, y, B, pEA, MDR){
  M = B * pEA * MDR * (lf+ec.lf)
  (a^2 * bc * y * M * (lf+ec.lf) )^0.5
}

# constant to keep lifespan from being numerically zero
# assume minimum survival time is half an hour
ec.lf = 1/48

# constant to keep PDR from being numerically zero
# assume the maximum EIP is 100 days
ec.pdr = 1/100


# Function to calculate mean & quantiles
calcPostQuants = function(input, grad.xs) {
  
  # Get length of gradient
  N.grad.xs <- length(grad.xs)
  
  # Create output dataframe
  output.df <- data.frame("mean" = numeric(N.Temp.xs), "median" = numeric(N.Temp.xs), "lowerCI" = numeric(N.Temp.xs), "upperCI" = numeric(N.Temp.xs), temp = grad.xs)
  
  # Calculate mean & quantiles
  for(i in 1:N.grad.xs){
    output.df$mean[i] <- mean(input[ ,i])
    output.df$lowerCI[i] <- quantile(input[ ,i], 0.025)
    output.df$upperCI[i] <- quantile(input[ ,i], 0.975)
    output.df$median[i] <- quantile(input[ ,i], 0.5)
  }
  
  output.df # return output
  
}



##########
###### 3. Calculate trait means across temp gradient
##########

bc.m = colMeans(bc.preds)
bite.rate.m = colMeans(bite.rate.preds)
EFD.1stGC.m = colMeans(EFD.1stGC.preds)
PDR.m = colMeans(PDR.preds)
est.bite.m = colMeans(est.bite.preds)
lf.m = colMeans(lifespan.preds)
lifetime.eggs.m = colMeans(lifetime.eggs.preds)
MDR.m = colMeans(MDR.preds)
mu.exp.m = colMeans(mu.exp.preds)
pEA.m = colMeans(pEA.preds)
gamma.m = colMeans(gamma.preds) #Shapiro-based; only for the 'MSP lifetime model'


##########
###### 4. Sensitivity Analysis #1
##########

### Full lifetime model
# Create matrices to hold sensitivity results
#dR0.da <- dR0.dbc <- dR0.dlf <- dR0.dPDR <- dR0.dB <- dR0.dpEA <- dR0.dMDR <- dR0.dT <- matrix(NA, nMCMC, N.Temp.xs)
dR0.da <- dR0.dbc <- dR0.dlf <- dR0.dgamma <- dR0.dB <- dR0.dpEA <- dR0.dMDR <- dR0.dT <- matrix(NA, nMCMC, N.Temp.xs)


# Calculate dy/dt and dR0/dy for each MCMC step across the temp gradient
for(i in 1:nMCMC){ # loop through MCMC steps
  # put in a step counter
  if(i%%500==0) cat("iteration ", i, "\n")
  
  # Calculate derivative of all traits w/r/t temp (dy/dt) across temp gradient (for a single MCMC step)
  da.dT <- d_briere(Temp.xs, 
                    bite.rate.out[[1]][i],
                    bite.rate.out[[2]][i],
                    bite.rate.out[[3]][i])
  dbc.dT <- d_quad(Temp.xs, 
                   bc.out[[1]][i],
                   bc.out[[2]][i],
                   bc.out[[3]][i])
  dlf.dT <- d_quad(Temp.xs, 
                   lifespan.out[[1]][i],
                   lifespan.out[[2]][i],
                   lifespan.out[[3]][i])
 dgamma.dT <- d_quad(Temp.xs,
                      gamma.out[[1]][i],
                      gamma.out[[2]][i],
                      gamma.out[[3]][i])
  dMDR.dT <- d_briere(Temp.xs,
                      MDR.out[[1]][i],
                      MDR.out[[2]][i],
                      MDR.out[[3]][i])
  dB.dT <- d_briere(Temp.xs,
                    lifetime.eggs.out[[1]][i],
                    lifetime.eggs.out[[2]][i],
                    lifetime.eggs.out[[3]][i])
  dpEA.dT <- d_quad(Temp.xs,
                    pEA.out[[1]][i],
                    pEA.out[[2]][i],
                    pEA.out[[3]][i])

  # Calculate sensitivity (dR0/dy * dy/dt) across temp gradient (for a single MCMC step)
  # See Mathematica notebook for derivative calculations
  dR0.da[i, ] <- R0.full(bite.rate.preds[i, ], bc.m, lf.m, gamma.m, lifetime.eggs.m, pEA.m, MDR.m)/(bite.rate.preds[i, ]+ec) * da.dT
  dR0.dbc[i, ] <- 1/2 * (R0.full(bite.rate.m, bc.preds[i, ], lf.m, gamma.m, lifetime.eggs.m, pEA.m, MDR.m)/(bc.preds[i, ]+ec) * dbc.dT)
  dR0.dlf[i, ] <- R0.full(bite.rate.m, bc.m, lifespan.preds[i, ], gamma.m, lifetime.eggs.m, pEA.m, MDR.m) /(lifespan.preds[i, ]+ec) * dlf.dT #Lifespan in twice; but gamma is used instead of e^[-1/{lf*PDR}]
  dR0.dgamma[i, ] <- 1/2 * (R0.full(bite.rate.m, bc.m, lf.m, gamma.preds[i, ], lifetime.eggs.m,pEA.m, MDR.m)/(gamma.preds[i, ]+ec) * dgamma.dT)
  dR0.dB[i, ] <- 1/2 * (R0.full(bite.rate.m, bc.m, lf.m, gamma.m, lifetime.eggs.preds[i, ], pEA.m, MDR.m)/(lifetime.eggs.preds[i, ]+ec) * dB.dT)
  dR0.dpEA[i, ] <- 1/2 * (R0.full(bite.rate.m, bc.m, lf.m, gamma.m, lifetime.eggs.m, pEA.preds[i, ], MDR.m)/(pEA.preds[i, ]+ec) * dpEA.dT)
  dR0.dMDR[i, ] <- 1/2 * (R0.full(bite.rate.m, bc.m, lf.m, gamma.m, lifetime.eggs.m, pEA.m, MDR.preds[i, ])/(MDR.preds[i, ]+ec) * dMDR.dT)
  dR0.dT[i, ] <-  dR0.da[i, ] + dR0.dbc[i, ] + dR0.dlf[i, ] + dR0.dgamma[i, ] + dR0.dB[i, ] + dR0.dpEA[i, ] + dR0.dMDR[i, ]

} # end MCMC loop

# Calculate R0 from trait means scaling
R0.med <- R0.full(bite.rate.m, bc.m, lf.m, gamma.m, lifetime.eggs.m, pEA.m, MDR.m)

dR0.R0da = calcPostQuants(dR0.da, Temp.xs)
dR0.R0dbc = calcPostQuants(dR0.dbc, Temp.xs)
dR0.R0dlf = calcPostQuants(dR0.dlf, Temp.xs)
dR0.R0dgamma = calcPostQuants(dR0.dgamma, Temp.xs)
dR0.R0dB = calcPostQuants(dR0.dB, Temp.xs)
dR0.R0dpEA = calcPostQuants(dR0.dpEA, Temp.xs)
dR0.R0dMDR = calcPostQuants(dR0.dMDR, Temp.xs)
dR0.R0dT = calcPostQuants(dR0.dT, Temp.xs)

# Plot posterior median sensitivities
#################
pdf("plots/sensitivity1_lifetime_model_inf.pdf")
plot(Temp.xs, dR0.R0dT$median, type = "l", lwd = 2, col = 1, xlab = expression(paste("Temperature (",degree,"C)")), xlim = c(14, 34), ylab = "Relative sensitivity")
abline(h = 0, col = "gray80")
lines(Temp.xs, dR0.R0da$median, lwd = 2, col = 2)
lines(Temp.xs, dR0.R0dbc$median, lwd = 2, col = 3)
lines(Temp.xs, dR0.R0dlf$median, lwd = 2, col = 4)
lines(Temp.xs, dR0.R0dgamma$median, lwd = 2, col = 5)
lines(Temp.xs, dR0.R0dB$median, lwd = 2, col = 6)
lines(Temp.xs, dR0.R0dpEA$median, lwd = 2, col = 7)
lines(Temp.xs, dR0.R0dMDR$median, lwd = 2, col = 8)
legend('topright', legend = c(expression(R[0]), "a", "bc", "lf", expression(gamma), "B", expression(p[EA]), "MDR"), col = c(1:8), lty = 1, lwd = 2, bty = 'n')
dev.off()
################





### Estimated trait model
# Create matrices to hold sensitivity results
dR0.da <- dR0.dbc <- dR0.dlf <- dR0.dPDR <- dR0.dEFD <- dR0.dpEA <- dR0.dMDR <- dR0.dT <- matrix(NA, nMCMC, N.Temp.xs)

# Calculate dy/dt and dR0/dy for each MCMC step across the temp gradient
for(i in 1:nMCMC){ # loop through MCMC steps
  # put in a step counter
  if(i%%500==0) cat("iteration ", i, "\n")
  
  # Calculate derivative of all traits w/r/t temp (dy/dt) across temp gradient (for a single MCMC step)
  da.dT <- d_briere(Temp.xs,
                    est.bite.out[[1]][i],
                    est.bite.out[[2]][i],
                    est.bite.out[[3]][i])
  dbc.dT <- d_quad(Temp.xs,
                   bc.out[[1]][i],
                   bc.out[[2]][i],
                   bc.out[[3]][i])
  dlf.dT <- d_quad(Temp.xs,
                   mu.exp.out[[1]][i],
                   mu.exp.out[[2]][i],
                   mu.exp.out[[3]][i])
  dPDR.dT <- d_briere(Temp.xs,
                      PDR.out[[1]][i],
                      PDR.out[[2]][i],
                      PDR.out[[3]][i])
  dMDR.dT <- d_briere(Temp.xs,
                      MDR.out[[1]][i],
                      MDR.out[[2]][i],
                      MDR.out[[3]][i])
  dEFD.dT <- d_briere(Temp.xs,
                      EFD.1stGC.out[[1]][i],
                      EFD.1stGC.out[[2]][i],
                      EFD.1stGC.out[[3]][i])
  dpEA.dT <- d_quad(Temp.xs,
                    pEA.out[[1]][i],
                    pEA.out[[2]][i],
                    pEA.out[[3]][i])
  
  # Calculate sensitivity (dR0/dy * dy/dt) across temp gradient (for a single MCMC step)
  # See Mathematica notebook for derivative calculations
  dR0.da[i, ] <- R0.est(est.bite.preds[i, ], bc.m, mu.exp.m, PDR.m, EFD.1stGC.m, pEA.m, MDR.m)/(est.bite.preds[i, ]+ec) * da.dT
  dR0.dbc[i, ] <- 1/2 * (R0.est(est.bite.m, bc.preds[i, ], mu.exp.m, PDR.m, EFD.1stGC.m, pEA.m, MDR.m)/(bc.preds[i, ]+ec) * dbc.dT)
  dR0.dlf[i, ] <- 1/2 * (R0.est(est.bite.m, bc.m, mu.exp.preds[i, ], PDR.m, EFD.1stGC.m, pEA.m, MDR.m) *
                           (1+3*PDR.m*mu.exp.preds[i, ]) / ((mu.exp.preds[i, ] + ec)^2 * (PDR.m + ec) ) * dlf.dT)
  dR0.dPDR[i, ] <- 1/2 * (R0.est(est.bite.m, bc.m, mu.exp.m, PDR.preds[i, ], EFD.1stGC.m, pEA.m, MDR.m)/((mu.exp.m + ec)*(PDR.preds[i, ] + ec)^2) * dPDR.dT)
  dR0.dEFD[i, ] <- 1/2 * (R0.est(est.bite.m, bc.m, mu.exp.m, PDR.m, EFD.1stGC.preds[i, ], pEA.m, MDR.m)/(EFD.1stGC.preds[i, ]+ec) * dEFD.dT)
  dR0.dpEA[i, ] <- 1/2 * (R0.est(est.bite.m, bc.m, mu.exp.m, PDR.m, EFD.1stGC.m, pEA.preds[i, ], MDR.m)/(pEA.preds[i, ]+ec) * dpEA.dT)
  dR0.dMDR[i, ] <- 1/2 * (R0.est(est.bite.m, bc.m, mu.exp.m, PDR.m, EFD.1stGC.m, pEA.m, MDR.preds[i, ])/(MDR.preds[i, ]+ec) * dMDR.dT)
  dR0.dT[i, ] <-  dR0.da[i, ] + dR0.dbc[i, ] + dR0.dlf[i, ] + dR0.dPDR[i, ] + dR0.dEFD[i, ] + dR0.dpEA[i, ] + dR0.dMDR[i, ]
  
} # end MCMC loop

# Calculate R0 from trait means scaling
R0.med <- R0.est(est.bite.m, bc.m, mu.exp.m, PDR.m, EFD.1stGC.m, pEA.m, MDR.m)

dR0.R0da2 = calcPostQuants(dR0.da, Temp.xs)
dR0.R0dbc2 = calcPostQuants(dR0.dbc, Temp.xs)
dR0.R0dlf2 = calcPostQuants(dR0.dlf, Temp.xs)
dR0.R0dPDR2 = calcPostQuants(dR0.dPDR, Temp.xs)
dR0.R0dEFD2 = calcPostQuants(dR0.dEFD, Temp.xs)
dR0.R0dpEA2 = calcPostQuants(dR0.dpEA, Temp.xs)
dR0.R0dMDR2 = calcPostQuants(dR0.dMDR, Temp.xs)
dR0.R0dT2 = calcPostQuants(dR0.dT, Temp.xs)

# Plot posterior median sensitivities
#################
pdf("plots/sensitivity1_Anestimated_model_inf.pdf")
plot(Temp.xs, dR0.R0dT2$median, type = "l", lwd = 2, col = 1, xlab = expression(paste("Temperature (",degree,"C)")), xlim = c(14, 34), ylab = "Relative sensitivity")
abline(h = 0, col = "gray80")
lines(Temp.xs, dR0.R0da2$median, lwd = 2, col = 2)
lines(Temp.xs, dR0.R0dbc2$median, lwd = 2, col = 3)
lines(Temp.xs, dR0.R0dlf2$median, lwd = 2, col = 4)
lines(Temp.xs, dR0.R0dPDR2$median, lwd = 2, col = 5)
lines(Temp.xs, dR0.R0dEFD2$median, lwd = 2, col = 6)
lines(Temp.xs, dR0.R0dpEA2$median, lwd = 2, col = 7)
lines(Temp.xs, dR0.R0dMDR2$median, lwd = 2, col = 8)
legend('topright', legend = c(expression(R[0]), "a", "bc", "lf", "PDR", "EFD", expression(p[EA]), "MDR"), col = c(1:8), lty = 1, lwd = 2, bty = 'n')
dev.off()
################

##########
###### 4. Sensitivity Analysis #2 - holding single parameters constant
##########

# Calculate R0 (lifetime formulation, informative priors except gamma) holding each parameter constant
R0.sens2.a = R0.full(1, bc.preds, lifespan.preds, gamma.preds, lifetime.eggs.preds, pEA.preds, MDR.preds)
R0.sens2.bc = R0.full(bite.rate.preds, 1, lifespan.preds, gamma.preds, lifetime.eggs.preds, pEA.preds, MDR.preds)
R0.sens2.lf = R0.full(bite.rate.preds, bc.preds, 1, gamma.preds, lifetime.eggs.preds, pEA.preds, MDR.preds)
#R0.sens2.PDR = R0.full(bite.rate.preds, bc.preds, lifespan.preds, 1, lifetime.eggs.preds, pEA.preds, MDR.preds)
R0.sens2.gamma = R0.full(bite.rate.preds, bc.preds, lifespan.preds, 1, lifetime.eggs.preds, pEA.preds, MDR.preds)
R0.sens2.fec = R0.full(bite.rate.preds, bc.preds, lifespan.preds, gamma.preds, 1, pEA.preds, MDR.preds)
R0.sens2.pEA = R0.full(bite.rate.preds, bc.preds, lifespan.preds, gamma.preds, lifetime.eggs.preds, 1, MDR.preds)
R0.sens2.MDR = R0.full(bite.rate.preds, bc.preds, lifespan.preds, gamma.preds, lifetime.eggs.preds, pEA.preds, 1)
R0.sens2 = R0.full(bite.rate.preds, bc.preds, lifespan.preds, gamma.preds, lifetime.eggs.preds, pEA.preds, MDR.preds)

# Calculate R0 (estimated formulation, informative priors) holding each parameter constant
R0.sens3.a = R0.est(1, bc.preds, mu.exp.preds, PDR.preds, EFD.1stGC.preds, pEA.preds, MDR.preds)
R0.sens3.bc = R0.est(est.bite.preds, 1, mu.exp.preds, PDR.preds, EFD.1stGC.preds, pEA.preds, MDR.preds)
R0.sens3.lf = R0.est(est.bite.preds, bc.preds, 1, PDR.preds, EFD.1stGC.preds, pEA.preds, MDR.preds)
R0.sens3.PDR = R0.est(est.bite.preds, bc.preds, mu.exp.preds, 1, EFD.1stGC.preds, pEA.preds, MDR.preds)
R0.sens3.fec = R0.est(est.bite.preds, bc.preds, mu.exp.preds, PDR.preds, 1, pEA.preds, MDR.preds)
R0.sens3.pEA = R0.est(est.bite.preds, bc.preds, mu.exp.preds, PDR.preds, EFD.1stGC.preds, 1, MDR.preds)
R0.sens3.MDR = R0.est(est.bite.preds, bc.preds, mu.exp.preds, PDR.preds, EFD.1stGC.preds, pEA.preds, 1)
R0.sens3 = R0.est(est.bite.preds, bc.preds, mu.exp.preds, PDR.preds, EFD.1stGC.preds, pEA.preds, MDR.preds)


# Get posterior quantiles for plotting
a.sens2.out <- calcPostQuants(R0.sens2.a, Temp.xs)
bc.sens2.out <- calcPostQuants(R0.sens2.bc, Temp.xs)
lf.sens2.out <- calcPostQuants(R0.sens2.lf, Temp.xs)
#PDR.sens2.out <- calcPostQuants(R0.sens2.PDR, Temp.xs)
gamma.sens2.out <- calcPostQuants(R0.sens2.gamma, Temp.xs)
fec.sens2.out <- calcPostQuants(R0.sens2.fec, Temp.xs)
pEA.sens2.out <- calcPostQuants(R0.sens2.pEA, Temp.xs)
MDR.sens2.out <- calcPostQuants(R0.sens2.MDR, Temp.xs)
R0.sens2.out <- calcPostQuants(R0.sens2, Temp.xs)

a.sens3.out <- calcPostQuants(R0.sens3.a, Temp.xs)
bc.sens3.out <- calcPostQuants(R0.sens3.bc, Temp.xs)
lf.sens3.out <- calcPostQuants(R0.sens3.lf, Temp.xs)
PDR.sens3.out <- calcPostQuants(R0.sens3.PDR, Temp.xs)
fec.sens3.out <- calcPostQuants(R0.sens3.fec, Temp.xs)
pEA.sens3.out <- calcPostQuants(R0.sens3.pEA, Temp.xs)
MDR.sens3.out <- calcPostQuants(R0.sens3.MDR, Temp.xs)
R0.sens3.out <- calcPostQuants(R0.sens3, Temp.xs)

##################
pdf("plots/sensitivity2_An_lifetime_model_inf.pdf", width = 10, height = 10)

##### Comparing R0s with single traits held constant - scaled, no CIs
plot(R0.sens2.out$mean / max(R0.sens2.out$mean) ~ Temp.xs, type = "l", main = "", ylim = c(0, 1.1), xlim = c(15,35),
     xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 0.9, cex.lab = 1.1, ylab = "", yaxt = "n", col = 1, lwd = 2)
lines(a.sens2.out$mean / max(a.sens2.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 2)
lines(bc.sens2.out$mean / max(bc.sens2.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 3)
lines(lf.sens2.out$mean / max(lf.sens2.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 4)
lines(gamma.sens2.out$mean / max(gamma.sens2.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 5)
lines(fec.sens2.out$mean / max(fec.sens2.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 6)
lines(pEA.sens2.out$mean / max(pEA.sens2.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 7)
lines(MDR.sens2.out$mean / max(MDR.sens2.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 8)
mtext(text = expression(paste('Relative R'[0])), side = 2, cex = 1.1, line = 1, las = 0)
# legend("topleft", legend = "A", bty = "n", adj = 1.5)
legend(x = 14.5, y = 1.025, col = c("black"), lwd = 3, lty = 1, legend = c("Lifetime model"), bty = "n", cex = 1)
text(x = 17, y = 0.85, "Constant parameter:", cex = 1)
legend(x = 14.5, y = 0.85, col = c(2:8), lwd = 1.5,
       lty = 1, legend = c("a", "bc", "lf", "PDR", "B", "pEA", "MDR"), bty = "n", cex = 1)

dev.off()
#####################

# Compare to the model built from estimated traits
################
pdf("plots/sensitivity2_Anestimated_model_inf.pdf", width = 10, height = 10)
##### Comparing R0s with single traits held constant - scaled, no CIs
plot(R0.sens3.out$mean / max(R0.sens3.out$mean) ~ Temp.xs, type = "l", main = "", ylim = c(0, 1.1), xlim = c(15,35),
     xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 0.9, cex.lab = 1.1, ylab = "", yaxt = "n", col = 1, lwd = 2)
lines(a.sens3.out$mean / max(a.sens3.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 2)
lines(bc.sens3.out$mean / max(bc.sens3.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 3)
lines(lf.sens3.out$mean / max(lf.sens3.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 4)
lines(PDR.sens3.out$mean / max(PDR.sens3.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 5)
lines(fec.sens3.out$mean / max(fec.sens3.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 6)
lines(pEA.sens3.out$mean / max(pEA.sens3.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 7)
lines(MDR.sens3.out$mean / max(MDR.sens3.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 8)
mtext(text = expression(paste('Relative R'[0])), side = 2, cex = 1.1, line = 1, las = 0)
# legend("topleft", legend = "A", bty = "n", adj = 1.5)
legend(x = 14.5, y = 1.025, col = c("black"), lwd = 3, lty = 1, legend = c("Estimated model"), bty = "n", cex = 1)
text(x = 17, y = 0.85, "Constant parameter:", cex = 1)
legend(x = 14.5, y = 0.85, col = c(2:8), lwd = 1.5,
       lty = 1, legend = c("a", "bc", "lf", "PDR", "B", "pEA", "MDR"), bty = "n", cex = 1)
dev.off()
########################



##########
###### 5. Uncertainty Analysis
##########

### Lifetime model
# Build matrices to hold results
R0.a<-R0.bc<-R0.fec<-R0.pEA<-R0.MDR<-R0.lf<-R0.gamma<-R0.lifetime<-matrix(NA,nrow(bite.rate.preds), ncol(bite.rate.preds))

## For uncertainty analysis: calculate posterior samples for R0 with
## all but a single component fixed the posterior mean.

for (j in 1:nrow(bite.rate.preds)){
  
  R0.a[j,] = R0.full(bite.rate.preds[j,], bc.m, lf.m, gamma.m, lifetime.eggs.m, pEA.m, MDR.m)
  R0.bc[j,] = R0.full(bite.rate.m, bc.preds[j,], lf.m, gamma.m, lifetime.eggs.m, pEA.m, MDR.m)
  R0.fec[j,] = R0.full(bite.rate.m, bc.m, lf.m, gamma.m, lifetime.eggs.preds[j,], pEA.m, MDR.m)
  R0.pEA[j,] = R0.full(bite.rate.m, bc.m, lf.m, gamma.m, lifetime.eggs.m, pEA.preds[j,], MDR.m)
  R0.MDR[j,] = R0.full(bite.rate.m, bc.m, lf.m, gamma.m, lifetime.eggs.m, pEA.m, MDR.preds[j,])
  R0.lf[j,] = R0.full(bite.rate.m, bc.m, lifespan.preds[j,], gamma.m, lifetime.eggs.m, pEA.m, MDR.m)
  #R0.PDR[j,] = R0.full(bite.rate.m, bc.m, lf.m, PDR.preds[j,], lifetime.eggs.m, pEA.m, MDR.m)
  R0.gamma[j,] = R0.full(bite.rate.m, bc.m, lf.m, gamma.preds[j,], lifetime.eggs.m, pEA.m, MDR.m)
  R0.lifetime[j,] = R0.full(bite.rate.preds[j,], bc.preds[j,], lifespan.preds[j,], gamma.preds[j,], lifetime.eggs.preds[j,], pEA.preds[j,], MDR.preds[j,])
  
}

## Calculate the distance within the inner 95% quantile for R0 overall
## (R0.q) and for the posterior of R0 with each component held fixed.

R0.q<-  apply(R0.lifetime, 2, FUN=quantile, probs=0.925, na.rm=F)- apply(R0.lifetime, 2, FUN=quantile, probs=0.025, na.rm=F)

a.q<-  apply(R0.a, 2, FUN=quantile, probs=0.925)- apply(R0.a, 2, FUN=quantile, probs=0.025)
bc.q<- apply(R0.bc, 2, FUN=quantile, probs=0.925)- apply(R0.bc, 2, FUN=quantile, probs=0.025)
fec.q<- apply(R0.fec, 2, FUN=quantile, probs=0.925)- apply(R0.fec, 2, FUN=quantile, probs=0.025)
pEA.q<-apply(R0.pEA, 2, FUN=quantile, probs=0.925)- apply(R0.pEA, 2, FUN=quantile, probs=0.025)
MDR.q<-  apply(R0.MDR, 2, FUN=quantile, probs=0.925)- apply(R0.MDR, 2, FUN=quantile, probs=0.025)
lf.q <-  apply(R0.lf, 2, FUN=quantile, probs=0.925)- apply(R0.lf, 2, FUN=quantile, probs=0.025)
gamma.q<- apply(R0.gamma, 2, FUN=quantile, probs=0.925)- apply(R0.gamma, 2, FUN=quantile, probs=0.025)

## Next plot relative width of quantiles 

# Creating a small constant to keep denominators from being zero.
ec<-0.1 

####################
pdf("plots/uncertainty_MSP_lifetime_model_inf.pdf")
plot(Temp.xs, a.q/(R0.q+ec), col=2, type="l", lwd=2, ylim = c(0,1), xlim = c(17, 34),
     xlab=expression(paste("Temperature (",degree,"C)")), ylab="Width of quantiles")
lines(Temp.xs, bc.q/(R0.q+ec), col=3, lwd=2)
lines(Temp.xs, lf.q/(R0.q+ec), col=4, lwd=2)
lines(Temp.xs, gamma.q/(R0.q+ec), col=5, lwd=2)
lines(Temp.xs, fec.q/(R0.q+ec), col=6, lwd=2)
lines(Temp.xs, pEA.q/(R0.q+ec), col=7, lwd=2)
lines(Temp.xs, MDR.q/(R0.q+ec), col=8, lwd=2)

# Adding a legend to the plot.
legend(18, 0.85, legend = c("a", "bc", "lf", "y", "B", expression(p[EA]), "MDR"), col = c(2:8), lty = 1, lwd = 2, bty = 'n')
dev.off()
###################

#### Repeat for the estimated model
# Build matrices to hold results
R0.a<-R0.bc<-R0.fec<-R0.pEA<-R0.MDR<-R0.lf<-R0.PDR<-R0.lifetime<-matrix(NA,nrow(est.bite.preds), ncol(est.bite.preds))

## For uncertainty analysis: calculate posterior samples for R0 with
## all but a single component fixed the posterior mean.

for (j in 1:nrow(est.bite.preds)){
  
  R0.a[j,] = R0.est(est.bite.preds[j,], bc.m, mu.exp.m, PDR.m, EFD.1stGC.m, pEA.m, MDR.m)
  R0.bc[j,] = R0.est(est.bite.m, bc.preds[j,], mu.exp.m, PDR.m, EFD.1stGC.m, pEA.m, MDR.m)
  R0.fec[j,] = R0.est(est.bite.m, bc.m, mu.exp.m, PDR.m, EFD.1stGC.preds[j,], pEA.m, MDR.m)
  R0.pEA[j,] = R0.est(est.bite.m, bc.m, mu.exp.m, PDR.m, EFD.1stGC.m, pEA.preds[j,], MDR.m)
  R0.MDR[j,] = R0.est(est.bite.m, bc.m, mu.exp.m, PDR.m, EFD.1stGC.m, pEA.m, MDR.preds[j,])
  R0.lf[j,] = R0.est(est.bite.m, bc.m, mu.exp.preds[j,], PDR.m, EFD.1stGC.m, pEA.m, MDR.m)
  R0.PDR[j,] = R0.est(est.bite.m, bc.m, mu.exp.m, PDR.preds[j,], EFD.1stGC.m, pEA.m, MDR.m)
  R0.lifetime[j,] = R0.est(est.bite.preds[j,], bc.preds[j,], mu.exp.preds[j,], PDR.preds[j,], EFD.1stGC.preds[j,], pEA.preds[j,], MDR.preds[j,])
  
}

## Calculate the distance within the inner 95% quantile for R0 overall
## (R0.q) and for the posterior of R0 with each component held fixed.

R0.q2<-  apply(R0.lifetime, 2, FUN=quantile, probs=0.925, na.rm=F)- apply(R0.lifetime, 2, FUN=quantile, probs=0.025, na.rm=F)

a.q2<-  apply(R0.a, 2, FUN=quantile, probs=0.925)- apply(R0.a, 2, FUN=quantile, probs=0.025)
bc.q2<- apply(R0.bc, 2, FUN=quantile, probs=0.925)- apply(R0.bc, 2, FUN=quantile, probs=0.025)
fec.q2<- apply(R0.fec, 2, FUN=quantile, probs=0.925)- apply(R0.fec, 2, FUN=quantile, probs=0.025)
pEA.q2<-apply(R0.pEA, 2, FUN=quantile, probs=0.925)- apply(R0.pEA, 2, FUN=quantile, probs=0.025)
MDR.q2<-  apply(R0.MDR, 2, FUN=quantile, probs=0.925)- apply(R0.MDR, 2, FUN=quantile, probs=0.025)
lf.q2 <-  apply(R0.lf, 2, FUN=quantile, probs=0.925)- apply(R0.lf, 2, FUN=quantile, probs=0.025)
PDR.q2<- apply(R0.PDR, 2, FUN=quantile, probs=0.925)- apply(R0.PDR, 2, FUN=quantile, probs=0.025)

## Next plot relative width of quantiles 

# Creating a small constant to keep denominators from being zero.
ec<-0.1 

####################
pdf("plots/uncertainty_MSP_estimated_model_inf.pdf")
plot(Temp.xs, a.q2/(R0.q2+ec), col=2, type="l", lwd=2, ylim = c(0,1), xlim = c(17, 34),
     xlab=expression(paste("Temperature (",degree,"C)")), ylab="Width of quantiles")
lines(Temp.xs, bc.q2/(R0.q2+ec), col=3, lwd=2)
lines(Temp.xs, lf.q2/(R0.q2+ec), col=4, lwd=2)
lines(Temp.xs, PDR.q2/(R0.q2+ec), col=5, lwd=2)
lines(Temp.xs, fec.q2/(R0.q2+ec), col=6, lwd=2)
lines(Temp.xs, pEA.q2/(R0.q2+ec), col=7, lwd=2)
lines(Temp.xs, MDR.q2/(R0.q2+ec), col=8, lwd=2)

# Adding a legend to the plot.
legend(18, 0.75, legend = c("a", "bc", "lf", "PDR", "EFD", expression(p[EA]), "MDR"), col = c(2:8), lty = 1, lwd = 2, bty = 'n')
dev.off()
###################


######################
######Create a final plot which displays four panels into a single figure
###################

# Alternative function that uses shaded polygons instead of dashed lines -> changed settings to be compatible with SI_FIGURES 2&3
plot.R0.poly = function(temp = Temp.xs, df, modname = "", add = "false", cols = "black", ltys = 1){
  if (add == "false") {
    par(mar=c(5,5,5,1))
    plot(temp, df$upperCI/max(df$upperCI, na.rm = T), lty = 1, xlab = expression(paste("Temperature (",degree,"C)")), ylab = expression("Relative R"[0]), xlim = c(15,35), cex.axis = 1.4, cex.lab = 1.4, main = modname, cex.main = 1.4, yaxt = 'n', col = "white")
    polygon(c(temp, rev(temp)), c(df$upperCI/max(df$upperCI, na.rm = T), rev(df$lowerCI/max(df$upperCI, na.rm = T))), col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  } else {
    polygon(c(temp, rev(temp)), c(df$upperCI/max(df$upperCI, na.rm = T), rev(df$lowerCI/max(df$upperCI, na.rm = T))), col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  }
  
  lines(temp, df$mean/max(df$upperCI, na.rm = T), lty = ltys, lwd = 2, col = cols)
  #points(temp, df$mean/max(df$upperCI, na.rm = T), col =cols, lwd =1) #Line to see why there is not a difference
}

#Note: R0.6.out required

########Supplemental Figure 2###########
pdf("plots/Anstephensi_lifetime_SIFIG2_sensitivity_uncertainty_revised414.pdf", width = 8, height = 8)
par(mfrow =c(2,2))

plot.R0.poly(df = R0.6.out, add = "false", cols = "black")
text1 <- c(expression(italic("An. stephensi")), expression("lifetime")) 
legend("topleft", legend=c(text1), lty=c(1,0,0,0),lwd=2, col = c("black"), y.intersp = 1,  bty="n") 
mtext("A", side = 3, at = 10, cex = 1.5, line = 2)

plot(Temp.xs, dR0.R0dT$median, type = "l", lwd = 2, col = 1, xlab = expression(paste("Temperature (",degree,"C)")), xlim = c(15, 34), ylim = c(min(dR0.R0dT$median), max(dR0.R0dT$median)),cex.axis = 1.4, cex.lab = 1.4, ylab = "Relative sensitivity")
abline(h = 0, col = "gray80")
lines(Temp.xs, dR0.R0da$median, lwd = 2, col = 2)
lines(Temp.xs, dR0.R0dbc$median, lwd = 2, col = 3)
lines(Temp.xs, dR0.R0dlf$median, lwd = 2, col = 4)
lines(Temp.xs, dR0.R0dgamma$median, lwd = 2, col = 5)
lines(Temp.xs, dR0.R0dB$median, lwd = 2, col = 6)
lines(Temp.xs, dR0.R0dpEA$median, lwd = 2, col = 7)
lines(Temp.xs, dR0.R0dMDR$median, lwd = 2, col = 8)
legend('bottomleft', ncol =2, legend = c(expression(R[0]), "a", "bc", "lf", expression(gamma), "B", expression(p[EA]), "MDR"), col = c(1:8), lty = 1, lwd = 2, bty = 'n')
mtext("B", side = 3, at = 10, cex = 1.5, line = 2)


plot(R0.sens2.out$mean / max(R0.sens2.out$mean) ~ Temp.xs, type = "l", main = "", ylim = c(0, 1.1), xlim = c(15,35),
     xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, ylab = expression(paste("Relative R"[0])), col = 1, lwd = 2)

lines(a.sens2.out$mean / max(a.sens2.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 2)
lines(bc.sens2.out$mean / max(bc.sens2.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 3)
lines(lf.sens2.out$mean / max(lf.sens2.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 4)
lines(gamma.sens2.out$mean / max(gamma.sens2.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 5)
lines(fec.sens2.out$mean / max(fec.sens2.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 6)
lines(pEA.sens2.out$mean / max(pEA.sens2.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 7)
lines(MDR.sens2.out$mean / max(MDR.sens2.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 8)
lines(R0.sens2.out$mean/ max(R0.sens2.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 1)
legend("topright", col = c("black"), lwd = 3, lty = 1, legend = c(expression(R[0])), bty = "n", cex = 1)
text(x = 19, y = 1.1, "Constant parameter:", cex = 1)
legend(x = 14, y = 1.1, col = c(2:8), lwd = 1.5,
       lty = 1, legend = c( "a", "bc", "lf",  expression(gamma), "B", "pEA", "MDR"), bty = "n", cex = 1)
mtext("C", side = 3, at = 10, cex = 1.5, line = 2)

plot(Temp.xs, a.q/(R0.q+ec), col=2, type="l", lwd=2, ylim = c(0,1), xlim = c(17, 34),cex.axis = 1.4, cex.lab = 1.4,
     xlab=expression(paste("Temperature (",degree,"C)")), ylab="Width of quantiles")
lines(Temp.xs, bc.q/(R0.q+ec), col=3, lwd=2)
lines(Temp.xs, lf.q/(R0.q+ec), col=4, lwd=2)
lines(Temp.xs, gamma.q/(R0.q+ec), col=5, lwd=2)
lines(Temp.xs, fec.q/(R0.q+ec), col=6, lwd=2)
lines(Temp.xs, pEA.q/(R0.q+ec), col=7, lwd=2)
lines(Temp.xs, MDR.q/(R0.q+ec), col=8, lwd=2)
# Adding a legend to the plot.
legend(18,.8,ncol =3, legend = c("a", "bc", "lf",  expression(gamma), "B", expression(p[EA]), "MDR"), col = c(2:8), lty = 1, lwd = 2, bty = 'n')
mtext("D", side = 3, at = 12, cex = 1.5, line = 2)

dev.off()


#Note R0.7.out required

#####################################
########Supplemental Figure 3###########
pdf("plots/Anstephensi_SIFIG3__estimated_sensitivity_uncertainty_revised414.pdf", width = 8, height = 8)
par(mfrow =c(2,2))

#Plot just the An. stephensi estimated model
#######################

plot.R0.poly(df = R0.7.out, add = "false", cols = "dodgerblue")
text1 <- c(expression(italic("An. stephensi")), expression("estimated")) 
legend("topleft", legend=c(text1), lty=c(1,0,0,0),lwd=2, col = c("dodgerblue"), y.intersp = 1,  bty="n") 
mtext("A", side = 3, at = 10, cex = 1.5, line = 2)

plot(Temp.xs, dR0.R0dT2$median, type = "l", lwd = 2, col = 1, xlab = expression(paste("Temperature (",degree,"C)")),cex.axis = 1.4, cex.lab = 1.4, xlim = c(15, 34), ylab = "Relative sensitivity")
abline(h = 0, col = "gray80")
lines(Temp.xs, dR0.R0da2$median, lwd = 2, col = 2)
lines(Temp.xs, dR0.R0dbc2$median, lwd = 2, col = 3)
lines(Temp.xs, dR0.R0dlf2$median, lwd = 2, col = 4)
lines(Temp.xs, dR0.R0dPDR2$median, lwd = 2, col = 5)
lines(Temp.xs, dR0.R0dEFD2$median, lwd = 2, col = 6)
lines(Temp.xs, dR0.R0dpEA2$median, lwd = 2, col = 7)
lines(Temp.xs, dR0.R0dMDR2$median, lwd = 2, col = 8)
legend('bottomleft', ncol =2,legend = c(expression(R[0]), "a*", "bc", "lf*", "PDR", "EFD*", expression(p[EA]), "MDR"), col = c(1:8), lty = 1, lwd = 2, bty = 'n')
mtext("B", side = 3, at = 10, cex = 1.5, line = 2)


plot(R0.sens3.out$mean / max(R0.sens3.out$mean) ~ Temp.xs, type = "l", main = "", ylim = c(0, 1.05), xlim = c(15,35),
     xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, ylab = expression(paste("Relative R"[0])), col = 1, lwd = 2)
lines(a.sens3.out$mean / max(a.sens3.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 2)
lines(bc.sens3.out$mean / max(bc.sens3.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 3)
lines(lf.sens3.out$mean / max(lf.sens3.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 4)
lines(PDR.sens3.out$mean / max(PDR.sens3.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 5)
lines(fec.sens3.out$mean / max(fec.sens3.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 6)
lines(pEA.sens3.out$mean / max(pEA.sens3.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 7)
lines(MDR.sens3.out$mean / max(MDR.sens3.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 8)
lines(R0.sens3.out$mean/ max(R0.sens3.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = 1)
# legend("topleft", legend = "A", bty = "n", adj = 1.5)
legend("topright", col = c("black"), lwd = 3, lty = 1, legend = c(expression(R[0])), bty = "n", cex = 1)
text(x = 19, y = 1, "Constant parameter:", cex = 1)
legend(x = 14, y = 1, col = c(2:8), lwd = 1.5,
       lty = 1, legend = c("a*", "bc", "lf*", "PDR", "EFD*", "pEA", "MDR"), bty = "n", cex = 1)
mtext("C", side = 3, at = 10, cex = 1.5, line = 2)

plot(Temp.xs, a.q2/(R0.q2+ec), col=2, type="l", lwd=2, ylim = c(0,1), xlim = c(17, 34), cex.axis = 1.4, cex.lab = 1.4,
     xlab=expression(paste("Temperature (",degree,"C)")), ylab="Width of quantiles")
lines(Temp.xs, bc.q2/(R0.q2+ec), col=3, lwd=2)
lines(Temp.xs, lf.q2/(R0.q2+ec), col=4, lwd=2)
lines(Temp.xs, PDR.q2/(R0.q2+ec), col=5, lwd=2)
lines(Temp.xs, fec.q2/(R0.q2+ec), col=6, lwd=2)
lines(Temp.xs, pEA.q2/(R0.q2+ec), col=7, lwd=2)
lines(Temp.xs, MDR.q2/(R0.q2+ec), col=8, lwd=2)
# Adding a legend to the plot.
legend("topleft",ncol =3, legend = c("a*", "bc", "lf*", "PDR", "EFD*", expression(p[EA]), "MDR"), col = c(2:8), lty = 1, lwd = 2, bty = 'n')
mtext("D", side = 3, at = 12, cex = 1.5, line = 2)
dev.off()

