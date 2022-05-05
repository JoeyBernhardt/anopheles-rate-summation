##Kerri Miazgowicz May 07,2020

## Purpose: Calculate temperature-dependent R0 from traits fit Use Bayesian Inference (JAGS)
##
## Contents: 1) Set-up,load packages, get data, etc.
##           3) Functions defining R0 equations
##           4) Calculate R0 posteriors and quantiles
##           5) Function to calculate distributions of T0, Tmax, and peak R0
##           6) Calculate distributions of T0, Tmax, and peak R0


##########
###### 1. Set up workspace, load packages, get data, etc.
##########

###set up working directory
#Set working directory
mainDir = "C:/Users/Kerri/Desktop/Chapter2 InProgress"
setwd(mainDir)

library('coda') # Calculate HPD intervals
library('dplyr') # Slice and dice data frames

#Load fits; use the un-summarised version of the fits

#load fits
load("saved posteriors/temps.Rdata")

#Constant temperature fits
load("saved posteriors/constant_lf_c.quad_T.uniform.Rdata")
# Calculate summary statistics of trait trajectories for plotting
lf.c.quad_T.summary <- trait.trajs.quad(lf.c.quad_T, Temp.xs, summary = FALSE)
load("saved posteriors/constant_br.c.briere_T.uniform.Rdata")
br.c.briere_T.summary <- trait.trajs.briere(br.c.briere_T, Temp.xs, summary = FALSE)
load("saved posteriors/constant_le.c.briere_T.uniform.Rdata")
le.c.briere_T.summary <- trait.trajs.briere(le.c.briere_T, Temp.xs, summary = FALSE)
load("saved posteriors/constant_gamma_c.quad_T.uniform.Rdata")
g.c.quad_T.summary<- trait.trajs.quad(gamma.c.quad_T,Temp.xs, summary = FALSE)


###literature based traits
load("saved posteriors/bc_quad_withT_uniform.Rdata")
bc.c.quad_T.summary <- trait.trajs.quad(bc.c.quad_T,Temp.xs,summary = FALSE)
load("saved posteriors/MDR_briere_withT_uniform.Rdata")
MDR.c.briere_T.summary <- trait.trajs.briere(MDR.c.briere_T, Temp.xs, summary = FALSE)
load("saved posteriors/PDR_briere_withT_uniform.Rdata")
PDR.c.briere_T.summary <- trait.trajs.briere(PDR.c.briere_T, Temp.xs,summary= FALSE)
load("saved posteriors/pea_quad_withT_uniform.Rdata")
pea.c.quad_T.summary <- trait.trajs.quad(pea.c.quad_T, Temp.xs, summary = FALSE)



##########
###### 2. calcPostQuants - Function to calculate quantiles of derived quantity posteriors (for plotting)
##########

# Arguments:
#   input       data frame with posterior distributions of trait predicted over a temperature gradient (i.e., 'trait.preds')
#   temp.list   the temperature gradient itself (i.e., 'Temp.xs')

calcPostQuants = function(input, temp.list) {
  
  # Get length of gradient
  N.grad.xs <- length(temp.list)
  
  # Create output dataframe
  output.df <- data.frame("mean" = numeric(N.grad.xs), "median" = numeric(N.grad.xs), "lowerCI" = numeric(N.grad.xs), "upperCI" = numeric(N.grad.xs))
  
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
###### 3. Functions defining R0 equations
##########

# Arguments: dataframes with posterior distributions of each trait predicted over a gradient (i.e., 'trait.preds')

# constant to keep lifespan from being numerically zero
# assume minimum survival time is half an hour
ec.lf = 1/48

# constant to keep PDR from being numerically zero
# assume the maximum EIP is 100 days
ec.pdr = 1/100

# 2. M depends on lifetime fecundity; gamma TPC is substituted for exp^(-1/(lf*PDR))expression
R0.full = function(a, bc, lf, y, B, pEA, MDR){
  M = B * pEA * MDR * (lf+ec.lf)
  (a^2 * bc * y * M * (lf+ec.lf) )^0.5
}


##########
###### 4. Function to calculate distributions of T0, Tmax, and peak R0
##########

# Arguments:
#   input       data frame with posterior distributions of R0 (or other quantity) predicted over a temperature gradient
#   temp.list   the temperature gradient itself (i.e., 'Temp.xs')

calcThreshPeakDists = function(input, temp.list) {
  
  # Create output dataframe
  output.df <- data.frame("peak" = numeric(nrow(input)), "T0" = numeric(nrow(input)), "Tmax" = numeric(nrow(input)))
  
  for (i in 1:nrow(input)) { # loop through each row of the input (MCMC step)
    
    output.df$peak[i] <- temp.list[which.max(input[i, ])] # Calculate index of R0 peak & store corresponding temperature
    
    # Makes lists of T0 and Tmax
    index.list <- which(input[i, ] > 0) # Create vector of list of indices where relative R0 > 0
    length.index.list <- length(index.list)
    output.df$T0[i] <- temp.list[index.list[1] - 1] # Store T0 (index prior to first value in index.list)
    output.df$Tmax[i] <- temp.list[index.list[length.index.list] + 1] # Store Tmax (index after to last value in index.list)
  }
  
  output.df # return
  
}

##########
###### 5. Calculate R0 posteriors and quantiles
##########
#R0.full = function(a, bc, lf, y, B, pEA, MDR)

# An. lifetime full constant model: Miazgowicz gomp survival; lifetime a,B,lf; Paaijmans et al. 2013 MDR and pea; Shapiro 2017 bc and EIP50
R0.constant = R0.full(br.c.briere_T.summary, bc.c.quad_T.summary, lf.c.quad_T.summary, g.c.quad_T.summary, le.c.briere_T.summary, pea.c.quad_T.summary, MDR.c.briere_T.summary)
#Test if functions work after I transpose the matrix
R0.constant =  t(R0.constant) #yup that fixed it.
R0.c.uni.out = calcPostQuants(R0.constant, Temp.xs)
R0.c.TPdists = calcThreshPeakDists(R0.constant, Temp.xs)

#Export the output so that rate summation can be used over it
#write.csv(R0.c.uni.out, "data/R0_constanst_PostQuants.csv")
#write.csv(R0.c.TPdists, "data/R0_c_TPdists.csv")
library(MASS)
#write.matrix(R0.constant, "data/R0_constant_preds.matrix")
##########
###### 6. Plot the resulting models
##########

# Alternative function that uses shaded polygons instead of dashed lines
plot.R0.poly = function(temp = Temp.xs, df, modname = "", add = "false", cols = "black"){
  if (add == "false") {
    par(mar=c(5,5,5,1))
    plot(temp, df$upperCI/max(df$upperCI, na.rm = T), lty = 1, xlim = c(0,40),xlab = expression(paste("Temperature (",degree,"C)")), ylab = expression(S[T]), cex.axis = 1.4, cex.lab = 1.4, main = modname, cex.main = 1.5, yaxt = 'n', col = "white")
    polygon(c(temp, rev(temp)), c(df$upperCI/max(df$upperCI, na.rm = T), rev(df$lowerCI/max(df$upperCI, na.rm = T))), col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  } else {
    polygon(c(temp, rev(temp)), c(df$upperCI/max(df$upperCI, na.rm = T), rev(df$lowerCI/max(df$upperCI, na.rm = T))), col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  }
  
  lines(temp, df$mean/max(df$upperCI, na.rm = T), lty = 1, lwd = 2, col = cols)
}

# Plot all the model means and credible intervals in different panels
#####################
pdf("plots/R0_constant_only.pdf", width = 4, height = 20)
par(mfrow = c(3, 1))
plot.R0.poly(df = R0.c.uni.out, modname = "Constant relative R0", cols = "black")
par(mfrow = c(1,1))
dev.off()
######################

#########################################
## After 04_RateSummation_ContiuousLP_withError.R script is run calculate R0 from the RS estimates

# An. lifetime full constant model: Miazgowicz gomp survival; lifetime a,B,lf; Paaijmans et al. 2013 MDR and pea; Shapiro 2017 bc and EIP50
R0.RS.dtr9 = R0.full(br.dtr9.sims.PostQuants, bc.dtr9.sims.PostQuants, lf.dtr9.sims.PostQuants, g.dtr9.sims.PostQuants, le.dtr9.sims.PostQuants, pea.dtr9.sims.PostQuants, mdr.dtr9.sims.PostQuants)
R0.RS.dtr9.allposteriors <- R0.full(br.dtr9.matrix.sims, bc.dtr9.matrix.sims, lf.dtr9.matrix.sims, g.dtr9.matrix.sims, le.dtr9.matrix.sims,pea.dtr9.matrix.sims,mdr.dtr9.matrix.sims)

#Test if functions work after I transpose the matrix
R0.RS.dtr9 =  t(dplyr::select(R0.RS.dtr9,-temps)) #yup that fixed it.-->Needed to drop the hour column to work properly
R0.RS.dtr9.out = calcPostQuants(R0.RS.dtr9, Temp.xs) #works as is, column order doesn't matter
R0.RS.dtr9.allpost.out <- calcPostQuants(t(R0.RS.dtr9.allposteriors), Temp.xs)

R0.RS.dtr9.TPdists = calcThreshPeakDists(R0.RS.dtr9, Temp.xs)
R0.RS.dtr9.allpost.TPdists = calcThreshPeakDists(R0.RS.dtr9.allposteriors, Temp.xs)

#Try plotting to see if it looks close at all?
plot.R0.poly(df = R0.RS.dtr9.out, modname = "RS dtr9 R0", cols = "blue") #Yes this is working correctly

####Repeat for a dtr of 12
R0.RS.dtr12 = R0.full(br.dtr12.sims.PostQuants, br.dtr12.sims.PostQuants, lf.dtr12.sims.PostQuants, g.dtr12.sims.PostQuants, le.dtr12.sims.PostQuants, pea.dtr12.sims.PostQuants, mdr.dtr12.sims.PostQuants)
R0.RS.dtr12.allposteriors <- R0.full(br.dtr12.matrix.sims, bc.dtr12.matrix.sims, lf.dtr12.matrix.sims, g.dtr12.matrix.sims, le.dtr12.matrix.sims,pea.dtr12.matrix.sims,mdr.dtr12.matrix.sims)

R0.RS.dtr12 =  t(dplyr::select(R0.RS.dtr12,-temps))
R0.RS.dtr12.out = calcPostQuants(R0.RS.dtr12, Temp.xs) #works as is, column order doesn't matter
R0.RS.dtr12.allpost.out <- calcPostQuants(t(R0.RS.dtr12.allposteriors), Temp.xs)

R0.RS.dtr12.TPdists = calcThreshPeakDists(R0.RS.dtr12, Temp.xs)
R0.RS.dtr12.allpost.TPdists = calcThreshPeakDists(R0.RS.dtr12.allposteriors, Temp.xs)


plot.R0.poly(df = R0.RS.dtr12.out, modname = "RS dtr12 R0", cols = "purple") #Yes this is working correctly
plot.R0.poly(df = R0.RS.dtr9.allpost.out, modname = "test", cols = "dodgerblue")
#Export the TPdists - rename the rownames
row.names(R0.RS.dtr9.TPdists)<- c("mean", "median","lowerCI","upperCI")
row.names(R0.RS.dtr9.TPdists)<- c("mean", "median","lowerCI", "upperCI")
#write.csv(R0.RS.dtr9.TPdists, "data/R0_RS_dtr9_TPdists.csv")
#write.csv(R0.RS.dtr12.TPdists, "data/R0_RS_dtr12_TPdists.csv")


#Customize this plotting function for constant temperature R0
plot.R0.poly3 = function(temp = Temp.xs, df, modname = "", add = "false", cols = "black"){
  if (add == "false") {
    par(mar=c(5,5,5,1))
    plot(temp, df$upperCI/max(df$upperCI, na.rm = T), lty = 1, xlim = c(0,40),xlab = expression(paste("Temperature (",degree,"C)")), ylab = expression("Suitability (S)"),ylim = c(0,max(df$mean/max(df$mean))*1.2), cex.axis = 1.4, cex.lab = 1.4, main = modname, cex.main = 1.5, yaxt = 'n', col = "white")
    polygon(c(temp, rev(temp)), c(df$upperCI/max(df$upperCI, na.rm = T), rev(df$lowerCI/max(df$upperCI, na.rm = T))), col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  } else {
    polygon(c(temp, rev(temp)), c(df$upperCI/max(df$upperCI, na.rm = T), rev(df$lowerCI/max(df$upperCI, na.rm = T))), col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  }
  
  lines(temp, df$mean/max(df$upperCI, na.rm = T), lty = 1, lwd = 2, col = cols)
}

#Go to other files for plotting RS estimates over constant temp TPC
pdf("plots/S_RS_over_constantTPC.pdf", width = 8, height = 10)
par(mfrow = c(4,3))
plot.R0.poly3(Temp.xs, df = R0.c.uni.out, cols = "black")
lines(Temp.xs, R0.RS.dtr9.out$median/max(R0.RS.dtr9.out$median), col='darkcyan', lty = 5, lwd = 2)
polygon(c(Temp.xs, rev(Temp.xs)), c(R0.RS.dtr9.out$upperCI/max(R0.RS.dtr9.out$median), rev(R0.RS.dtr9.out$lowerCI/max(R0.RS.dtr9.out$median))), col = adjustcolor("darkcyan", alpha.f = 0.2), border = NA)

lines(Temp.xs, R0.RS.dtr12.out$median/max(R0.RS.dtr12.out$median), col='#4a3fa8', lty = 5, lwd = 2)
polygon(c(Temp.xs, rev(Temp.xs)), c(R0.RS.dtr12.out$upperCI/max(R0.RS.dtr12.out$median), rev(R0.RS.dtr12.out$lowerCI/max(R0.RS.dtr12.out$median))), col = adjustcolor("#4a3fa8", alpha.f = 0.2), border = NA)
mtext("h", side = 3, at = -5, cex = 1.8, line = 2)
par(mfrow = c(1,1))
dev.off()

####Main Figure 2d- thermal suitability
# Plot thermal performance curves for observed data over complete RS estimates
pdf("plots/Main_Fig2_S_overlay_RSwithError.pdf", width = 8, height = 8)
par(mfrow = c(3,2))
plot.R0.poly3(Temp.xs, df = R0.dtr9.uni.out, cols = "darkcyan")
plot.R0.poly3(Temp.xs, df = R0.dtr12.uni.out, cols = "#4a3fa8", add ="true")

#RS estimates associated with each S -dtr9
lines(Temp.xs, R0.RS.dtr9.out$median/max(R0.RS.dtr9.out$median), col='darkcyan', lty = 2, lwd = 2)
lines(Temp.xs, R0.RS.dtr9.out$upperCI/max(R0.RS.dtr9.out$median), lty = 3, lwd = 1, col = adjustcolor("darkcyan", alpha.f = 0.5))
lines(Temp.xs, R0.RS.dtr9.out$lowerCI/max(R0.RS.dtr9.out$median), lty = 3, lwd = 1, col = adjustcolor("darkcyan", alpha.f = 0.5))
polygon(c(Temp.xs, rev(Temp.xs)), c(R0.RS.dtr9.out$upperCI/max(R0.RS.dtr9.out$median), rev(R0.RS.dtr9.out$lowerCI/max(R0.RS.dtr9.out$median))), col = adjustcolor("darkcyan", alpha.f = 0.2), border = NA)

#RS estimates associated with each S-dtr12
lines(Temp.xs, R0.RS.dtr12.out$median/max(R0.RS.dtr12.out$median), col='#4a3fa8', lty = 2, lwd = 2)
lines(Temp.xs, R0.RS.dtr12.out$upperCI/max(R0.RS.dtr12.out$median), lty = 3, lwd = 1, col = adjustcolor("#4a3fa8", alpha.f = 0.5))
lines(Temp.xs, R0.RS.dtr12.out$lowerCI/max(R0.RS.dtr12.out$median), lty = 3, lwd = 1, col = adjustcolor("#4a3fa8", alpha.f = 0.5))
polygon(c(Temp.xs, rev(Temp.xs)), c(R0.RS.dtr12.out$upperCI/max(R0.RS.dtr12.out$median), rev(R0.RS.dtr12.out$lowerCI/max(R0.RS.dtr12.out$median))), col = adjustcolor("#4a3fa8", alpha.f = 0.2), border = NA)

mtext("d", side = 3, at = -5, cex = 1.8, line = 2)
par(mfrow = c(1,1))
dev.off()

