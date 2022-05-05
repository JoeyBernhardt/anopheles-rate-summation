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


##temp fluctuation dtr 9
load("saved posteriors/dtr_lf_f9.quad_T.uniform.Rdata")
lf.f9.quad_T.summary <- trait.trajs.quad(lf.f9.quad_T, Temp.xs, summary = FALSE)
load("saved posteriors/dtr_br.f9.briere_T.uniform.Rdata")
br.f9.briere_T.summary <- trait.trajs.briere(br.f9.briere_T, Temp.xs, summary = FALSE)
load("saved posteriors/dtr_le.f9.briere_T.uniform.Rdata")
le.f9.briere_T.summary <- trait.trajs.briere(le.f9.briere_T, Temp.xs, summary = FALSE)
load("saved posteriors/dtr_gamma_f9.quad_T.uniform.Rdata")
g.f9.quad_T.summary<- trait.trajs.quad(gamma.f9.quad_T,Temp.xs, summary = FALSE)

#To save space remove all the objects I do not need
remove(lf.f9.quad_T)
remove(br.f9.briere_T)
remove(le.f9.briere_T)
remove(gamma.f9.quad_T)

##temp fluctuation dtr 12
load("saved posteriors/dtr_lf_f12.quad_T.uniform.Rdata")
lf.f12.quad_T.summary <- trait.trajs.quad(lf.f12.quad_T, Temp.xs, summary = FALSE)
load("saved posteriors/dtr_br.f12.briere_T.uniform.Rdata")
br.f12.briere_T.summary <- trait.trajs.briere(br.f12.briere_T,Temp.xs, summary = FALSE)
load("saved posteriors/dtr_le.f12.briere_T.uniform.Rdata")
le.f12.briere_T.summary <- trait.trajs.briere(le.f12.briere_T, Temp.xs, summary = FALSE)
load("saved posteriors/dtr_gamma_f12.quad_T.uniform.Rdata")
g.f12.quad_T.summary<- trait.trajs.quad(gamma.f12.quad_T,Temp.xs, summary = FALSE)

#To save space remove all the objects I do not need
remove(lf.f12.quad_T)
remove(br.f12.briere_T)
remove(le.f12.briere_T)
remove(gamma.f12.quad_T)

###literature based traits
load("saved posteriors/bc_quad_withT_uniform.Rdata")
bc.c.quad_T.summary <- trait.trajs.quad(bc.c.quad_T,Temp.xs,summary = FALSE)
load("saved posteriors/MDR_briere_withT_uniform.Rdata")
MDR.c.briere_T.summary <- trait.trajs.briere(MDR.c.briere_T, Temp.xs, summary = FALSE)
load("saved posteriors/PDR_briere_withT_uniform.Rdata")
PDR.c.briere_T.summary <- trait.trajs.briere(PDR.c.briere_T, Temp.xs,summary= FALSE)
load("saved posteriors/pea_quad_withT_uniform.Rdata")
pea.c.quad_T.summary <- trait.trajs.quad(pea.c.quad_T, Temp.xs, summary = FALSE)

remove(bc.c.quad_T)
remove(MDR.c.briere_T)
remove(PDR.c.briere_T)
remove(pea.c.quad_T)

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

# An. lifetime full dtr9 model: Miazgowicz gomp survival; lifetime a,B,lf; Paaijmans et al. 2013 MDR and pea; Shapiro 2017 bc and EIP50
R0.dtr9 = R0.full(br.f9.briere_T.summary, bc.c.quad_T.summary, lf.f9.quad_T.summary, g.f9.quad_T.summary, le.f9.briere_T.summary, pea.c.quad_T.summary, MDR.c.briere_T.summary)
#Test if functions work after I transpose the matrix
R0.dtr9 =  t(R0.dtr9) #yup that fixed it.
R0.dtr9.uni.out = calcPostQuants(R0.dtr9, Temp.xs)
R0.dtr9.TPdists = calcThreshPeakDists(R0.dtr9, Temp.xs)

#Export the output so that rate summation can be used over it
#write.csv(R0.dtr9.uni.out, "data/R0_dtr9_PostQuants.csv")
#write.csv(R0.dtr9.TPdists, "data/R0_dtr9_TPdists.csv")

# An. lifetime full dtr12 model: Miazgowicz gomp survival; lifetime a,B,lf; Paaijmans et al. 2013 MDR and pea; Shapiro 2017 bc and EIP50
R0.dtr12 = R0.full(br.f12.briere_T.summary, bc.c.quad_T.summary, lf.f12.quad_T.summary, g.f12.quad_T.summary, le.f12.briere_T.summary, pea.c.quad_T.summary, MDR.c.briere_T.summary)
#Test if functions work after I transpose the matrix
R0.dtr12 =  t(R0.dtr12) #yup that fixed it.
R0.dtr12.uni.out = calcPostQuants(R0.dtr12, Temp.xs)
R0.dtr12.TPdists = calcThreshPeakDists(R0.dtr12, Temp.xs)

#Export the output so that rate summation can be used over it
#write.csv(R0.dtr12.uni.out, "data/R0_dtr12_PostQuants.csv")
#write.csv(R0.dtr12.TPdists, "data/R0_dtr12_TPdists.csv")

##########
###### 6. Plot the resulting models
##########

# Alternative function that uses shaded polygons instead of dashed lines
plot.R0.poly = function(temp = Temp.xs, df, modname = "", add = "false", cols = "black"){
  if (add == "false") {
    par(mar=c(5,5,5,1))
    plot(temp, df$upperCI/max(df$upperCI, na.rm = T), lty = 1, xlab = expression(paste("Temperature (",degree,"C)")), ylab = expression(R[0]), cex.axis = 1.2, cex.lab = 1.2, main = modname, cex.main = 1.5, yaxt = 'n', col = "white")
    polygon(c(temp, rev(temp)), c(df$upperCI/max(df$upperCI, na.rm = T), rev(df$lowerCI/max(df$upperCI, na.rm = T))), col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  } else {
    polygon(c(temp, rev(temp)), c(df$upperCI/max(df$upperCI, na.rm = T), rev(df$lowerCI/max(df$upperCI, na.rm = T))), col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  }
  
  lines(temp, df$mean/max(df$upperCI, na.rm = T), lty = 1, lwd = 2, col = cols)
}

# Plot all the model means and credible intervals in different panels
#####################
pdf("plots/R0_dtr_each.pdf", width = 4, height = 20)
par(mfrow = c(3, 1))
plot.R0.poly(df = R0.dtr9.uni.out, modname = "DTR9 relative R0", cols = "darkcyan")
plot.R0.poly(df = R0.dtr12.uni.out, modname = "DTR12 relative R0", cols = "#4a3fa8")

par(mfrow = c(1,1))
dev.off()
######################

