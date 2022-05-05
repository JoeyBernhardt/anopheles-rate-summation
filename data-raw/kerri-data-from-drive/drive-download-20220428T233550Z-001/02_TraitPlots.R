### Code written by Kerri Miazgowicz on 4/22

##Standardize the plotting of the thermal responses

###set up working directory
#Set working directory
mainDir = "C:/Users/Kerri/Desktop/Chapter2 InProgress"
setwd(mainDir)

######  Load raw trait data
cdata <- data.frame(read.csv("data/constant.individual.trait.csv"))
c.gamma.data <- data.frame(read.csv("data/constant.gamma.values.csv"))

f9data <- data.frame(read.csv("data/fluc9.individual.trait.csv"))
f9.gamma.data <- data.frame(read.csv("data/dtr9.gamma.values.csv"))

f12data <- data.frame(read.csv("data/fluc12.individual.trait.csv"))
f12.gamma.data <- data.frame(read.csv("data/dtr12.gamma.values.csv"))

pea.MDR.data <- data.frame(read.csv("data/Krijn_Raw_Data.csv"))
bc.PDR.data <- data.frame(read.csv("data/forErin_ShapiroData.csv"))

#######Generate Block specific data values for fitting the curve over 
library(plyr)
library(dplyr)
library(tidyr)

c.data.block.temp <- data.frame(cdata %>%
  dplyr::group_by(Treatment,Block)%>%
  dplyr::summarise(lifespan = mean(lifespan),
                   bite.rate = mean(bite.rate),
                   lifetime.eggs = mean(lifetime.eggs))) %>%
  ungroup()

plot(lifespan ~ Treatment, data = c.data.block.temp) 
plot(bite.rate ~ Treatment, data = c.data.block.temp) 
plot(lifetime.eggs ~ Treatment, data = c.data.block.temp) 

################dtr 9#################
f9.data.block.temp <- data.frame(f9data %>%
  dplyr::group_by(Treatment,Block)%>%
  dplyr::summarise(lifespan = mean(lifespan),
                   bite.rate = mean(bite.rate),
                   lifetime.eggs = mean(lifetime.eggs))) %>%
  ungroup()

plot(lifespan ~ Treatment, data = f9.data.block.temp) 
plot(bite.rate ~ Treatment, data = f9.data.block.temp) 
plot(lifetime.eggs ~ Treatment, data = f9.data.block.temp) 
###############dtr 12##################3
f12.data.block.temp <- data.frame(f12data %>%
  dplyr::group_by(Treatment,Block)%>%
  dplyr::summarise(lifespan = mean(lifespan),
                   bite.rate = mean(bite.rate),
                   lifetime.eggs = mean(lifetime.eggs))) %>%
  ungroup()

plot(lifespan ~ Treatment, data = f12.data.block.temp) 
plot(bite.rate ~ Treatment, data = f12.data.block.temp) 
plot(lifetime.eggs ~ Treatment, data = f12.data.block.temp) 
plot(gamma~temp, data = c.gamma.data)
plot(gamma~temp, data = f9.gamma.data)
plot(gamma~temp, data = f12.gamma.data)


#load fits
load("saved posteriors/temps.Rdata")

#Constant temperature fits
load("saved posteriors/constant_lf_c.quad_T.uniform.Rdata")
# Calculate summary statistics of trait trajectories for plotting
lf.c.quad_T.summary <- trait.trajs.quad(lf.c.quad_T, Temp.xs, summary = TRUE)
load("saved posteriors/constant_br.c.briere_T.uniform.Rdata")
br.c.briere_T.summary <- trait.trajs.briere(br.c.briere_T, Temp.xs, summary = TRUE)
load("saved posteriors/constant_le.c.briere_T.uniform.Rdata")
le.c.briere_T.summary <- trait.trajs.briere(le.c.briere_T, Temp.xs, summary = TRUE)
load("saved posteriors/constant_gamma_c.quad_T.uniform.Rdata")
g.c.quad_T.summary<- trait.trajs.quad(gamma.c.quad_T,Temp.xs, summary = TRUE)

##temp fluctuation dtr 9
load("saved posteriors/dtr_lf_f9.quad_T.uniform.Rdata")
lf.f9.quad_T.summary <- trait.trajs.quad(lf.f9.quad_T, Temp.xs, summary = TRUE)
load("saved posteriors/dtr_br.f9.briere_T.uniform.Rdata")
br.f9.briere_T.summary <- trait.trajs.briere(br.f9.briere_T, Temp.xs, summary = TRUE)
load("saved posteriors/dtr_le.f9.briere_T.uniform.Rdata")
le.f9.briere_T.summary <- trait.trajs.briere(le.f9.briere_T, Temp.xs, summary = TRUE)
load("saved posteriors/dtr_gamma_f9.quad_T.uniform.Rdata")
g.f9.quad_T.summary<- trait.trajs.quad(gamma.f9.quad_T,Temp.xs, summary = TRUE)



##temp fluctuation dtr 12
load("saved posteriors/dtr_lf_f12.quad_T.uniform.Rdata")
lf.f12.quad_T.summary <- trait.trajs.quad(lf.f12.quad_T, Temp.xs, summary = TRUE)
load("saved posteriors/dtr_br.f12.briere_T.uniform.Rdata")
br.f12.briere_T.summary <- trait.trajs.briere(br.f12.briere_T,Temp.xs, summary = TRUE)
load("saved posteriors/dtr_le.f12.briere_T.uniform.Rdata")
le.f12.briere_T.summary <- trait.trajs.briere(le.f12.briere_T, Temp.xs, summary = TRUE)
load("saved posteriors/dtr_gamma_f12.quad_T.uniform.Rdata")
g.f12.quad_T.summary<- trait.trajs.quad(gamma.f12.quad_T,Temp.xs, summary = TRUE)


###literature based traits
load("saved posteriors/bc_quad_withT_uniform.Rdata")
bc.c.quad_T.summary <- trait.trajs.quad(bc.c.quad_T,Temp.xs,summary = TRUE)
load("saved posteriors/MDR_briere_withT_uniform.Rdata")
MDR.c.briere_T.summary <- trait.trajs.briere(MDR.c.briere_T, Temp.xs, summary = TRUE)
load("saved posteriors/PDR_briere_withT_uniform.Rdata")
PDR.c.briere_T.summary <- trait.trajs.briere(PDR.c.briere_T, Temp.xs,summary= TRUE)
load("saved posteriors/pea_quad_withT_uniform.Rdata")
pea.c.quad_T.summary <- trait.trajs.quad(pea.c.quad_T, Temp.xs, summary = TRUE)


#load datapoints
#############
####### 2. Plot thermal performance curves
#############
###### Plot TPCs with the raw data


trait.plots = function(traitname, summary, df, labs, pch = 1, cols = c("black"), add ="false"){
  # Specify the temperature sequences and trait to plot
  par(mar=c(5,6,5,1))
  temp = df$Treatment
  traitplot = (df[ , which(colnames(df)==traitname)])

  # plot the data
  if(add== "false"){
    plot(traitplot ~ temp, xlim = c(0, 45), ylim = c(0, max(traitplot, na.rm = T)*1.02), 
         ylab = labs, xlab = expression(paste("Temperature (",degree,"C)")), cex.lab = 1.4, cex.axis = 1.4, pch = pch)
  }else{
    points(traitplot ~ temp, data = df, pch = pch, col = cols)
  }
  # plot the TPC mean and 95% credible interval
  lines(summary$lower ~ summary$temp, lty = 2, col = cols)
  lines(summary$upper ~ summary$temp, lty = 2, col = cols)
  lines(summary$mean ~ summary$temp, lwd =2, col = cols)
}



##############################################
##First I need to redefine the mapping function
###########################################
pea.MDR.data$Treatment <- pea.MDR.data$temp
bc.PDR.data$Treatment <- bc.PDR.data$temp
bc.PDR.data$PDR <- 1/bc.PDR.data$EIP50
c.gamma.data$Treatment <- c.gamma.data$temp
f9.gamma.data$Treatment <- f9.gamma.data$temp
f12.gamma.data$Treatment <- f12.gamma.data$temp

# Plot thermal performance curves with data they were fit over
pdf("plots/trait_fits_final_means.pdf", width = 8, height = 10)
par(mfrow = c(4,3))
trait.plots("lifespan",lf.c.quad_T.summary , c.data.block.temp, "Lifespan(days)", cols = "black", pch = 1)
mtext("A", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots("bite.rate",br.c.briere_T.summary, c.data.block.temp, "Bite rate (days^-1)", cols = "black", pch = 1)
mtext("B", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots("lifetime.eggs",le.c.briere_T.summary, c.data.block.temp, "Fecundity (eggs)", cols = "black", pch = 1)
mtext("C", side = 3, at = -5, cex = 1.8, line = 2)

trait.plots("lifespan",lf.f9.quad_T.summary , f9.data.block.temp, "Lifespan(days)", cols = "darkcyan", pch = 1)
mtext("D", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots("bite.rate",br.f9.briere_T.summary , f9.data.block.temp, "Bite rate(days^-1)", cols = "darkcyan", pch = 1)
mtext("E", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots("lifetime.eggs",le.f9.briere_T.summary, f9.data.block.temp, "Fecundity (eggs)", cols = "darkcyan", pch = 1)
mtext("F", side = 3, at = -5, cex = 1.8, line = 2)

trait.plots("lifespan",lf.f12.quad_T.summary , f12.data.block.temp, "Lifespan(days)", cols = "#4a3fa8", pch = 1)
mtext("G", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots("bite.rate",br.f12.briere_T.summary , f12.data.block.temp, "Bite rate(days^1)", cols = "#4a3fa8", pch = 1)
mtext("H", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots("lifetime.eggs",le.f12.briere_T.summary, f12.data.block.temp, "Fecundity (eggs)", cols = "#4a3fa8", pch = 1)
mtext("I", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots("Pea", pea.c.quad_T.summary,pea.MDR.data,"Prob. e2a survival", cols = "#AD343E", pch=1 )
mtext("J", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots("MDR", MDR.c.briere_T.summary, pea.MDR.data,"Mosq. development rate (days-1)", cols = "#AD343E", pch=1 )
mtext("K", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots("PDR", PDR.c.briere_T.summary, bc.PDR.data,"Path. development rate(days-1)", cols = "#AD343E", pch=1 )
mtext("L", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots("bc", bc.c.quad_T.summary,bc.PDR.data,"Vector competence", cols = "#AD343E", pch=1 )
mtext("M", side = 3, at = -5, cex = 1.8, line = 2)

trait.plots("gamma",g.c.quad_T.summary, c.gamma.data, "Gamma", cols = "black", pch = 1)
mtext("N", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots("gamma",g.f9.quad_T.summary, f9.gamma.data, "Gamma", cols = "darkcyan", pch = 1)
mtext("O", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots("gamma",g.f12.quad_T.summary, f12.gamma.data, "Gamma", cols = "#4a3fa8", pch = 1)
mtext("P", side = 3, at = -5, cex = 1.8, line = 2)

par(mfrow = c(1,1))
dev.off()

#####################Alternative plotting function
#credible intervals as shaded lines instead of dashed lines

trait.plots.poly = function(traitname, summary,df, add = "false", cols = "black", pch = 1, labs = ""){
  temp = df$Treatment
  traitplot = (df[ , which(colnames(df)==traitname)])
  if (add == "false") {
    par(mar=c(5,5,5,1))
    plot(traitplot ~ temp, lty = 1,xlim = c(0,45), xlab = expression(paste("Temperature (",degree,"C)")), ylab = labs, cex.axis = 1.4, cex.lab = 1.4, ylim = c(0, max(traitplot, na.rm = T)*1.02), cex.main = 1.5, col = cols, pch = pch)
    polygon(c(summary$temp, rev(summary$temp)), c(summary$upper, rev(summary$lower)), col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  } else {
    polygon(c(summary$temp, rev(summary$temp)), c(summary$upper, rev(summary$lower)), col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  }
  points(traitplot~temp, col = cols, pch = pch)
  lines(summary$temp, summary$mean, lty = 1, lwd = 2, col =  adjustcolor(cols, alpha.f = 1))
}


# Plot thermal performance curves with data they were fit over
pdf("plots/trait_fits_poly_means.pdf", width = 8, height = 10)
par(mfrow = c(4,3))
trait.plots.poly("lifespan",lf.c.quad_T.summary , c.data.block.temp, "Lifespan (days)", cols = "black", pch = 1, add = "false")
trait.plots.poly("bite.rate",br.c.briere_T.summary, c.data.block.temp, "Bite rate (days^-1)", cols = "black", pch = 1, add = "false")
trait.plots.poly("lifetime.eggs",le.c.briere_T.summary, c.data.block.temp, "Fecundity (eggs)", cols = "black", pch = 1, add = "false")

trait.plots.poly("lifespan",lf.f9.quad_T.summary , f9.data.block.temp, "Lifespan (days)", cols = "darkcyan", pch = 1, add ="false")
trait.plots.poly("bite.rate",br.f9.briere_T.summary , f9.data.block.temp, "Bite rate(days^-1)", cols = "darkcyan", pch = 1, add = "false")
trait.plots.poly("lifetime.eggs",le.f9.briere_T.summary, f9.data.block.temp, "Fecundity (eggs)", cols = "darkcyan", pch = 1, add = "false")

trait.plots.poly("lifespan",lf.f12.quad_T.summary , f12.data.block.temp, "Lifespan (days)", cols = "#4a3fa8", pch = 1, add = "false")
trait.plots.poly("bite.rate",br.f12.briere_T.summary , f12.data.block.temp, "Bite rate(days^1)", cols = "#4a3fa8", pch = 1, add = "false")
trait.plots.poly("lifetime.eggs",le.f12.briere_T.summary, f12.data.block.temp, "Fecundity (eggs)", cols = "#4a3fa8", pch = 1, add = "false")

trait.plots.poly("Pea", pea.c.quad_T.summary, pea.MDR.data,"Prob. e2a survival", cols = "#AD343E", pch=1, add = "false" )
trait.plots.poly("MDR", MDR.c.briere_T.summary, pea.MDR.data,"Mosq. dev. rate (days-1)", cols = "#AD343E", pch=1, add = "false" )
trait.plots.poly("PDR", PDR.c.briere_T.summary, bc.PDR.data,"Path. dev. rate(days-1)", cols = "#AD343E", pch=1, add = "false" )
trait.plots.poly("bc", bc.c.quad_T.summary,bc.PDR.data,"Vector competence", cols = "#AD343E", pch=1, add = "false" )

trait.plots.poly("gamma",g.c.quad_T.summary, c.gamma.data, "Gamma", cols = "black", pch = 1, add = "false")
trait.plots.poly("gamma",g.f9.quad_T.summary, f9.gamma.data, "Gamma", cols = "darkcyan", pch = 1, add = "false")
trait.plots.poly("gamma",g.f12.quad_T.summary, f12.gamma.data, "Gamma", cols = "#4a3fa8", pch = 1, add = "false")

par(mfrow = c(1,1))
dev.off()

####Overlay by trait figure
# Plot thermal performance curves with data they were fit over
pdf("plots/trait_fits_overlay_poly_means.pdf", width = 4, height = 8)
par(mfrow = c(3,1))
trait.plots.poly("lifespan",lf.c.quad_T.summary , c.data.block.temp, "Lifespan (days)", cols = "black", pch = 16, add = "false")
trait.plots.poly("lifespan",lf.f9.quad_T.summary , f9.data.block.temp, "Lifespan (days)", cols = "darkcyan", pch = 16, add ="true")
trait.plots.poly("lifespan",lf.f12.quad_T.summary , f12.data.block.temp, "Lifespan (days)", cols = "#4a3fa8", pch = 16, add = "true")
trait.plots.poly("bite.rate",br.c.briere_T.summary, c.data.block.temp, expression(paste("Bite rate (days"^"-1", ")")), cols = "black", pch = 16, add = "false")
trait.plots.poly("bite.rate",br.f9.briere_T.summary , f9.data.block.temp, expression(paste("Bite rate (days"^"-1", ")")), cols = "darkcyan", pch = 16, add = "true")
trait.plots.poly("bite.rate",br.f12.briere_T.summary , f12.data.block.temp, expression(paste("Bite rate (days"^"-1", ")")), cols = "#4a3fa8", pch = 16, add = "true")
trait.plots.poly("lifetime.eggs",le.c.briere_T.summary, c.data.block.temp, "Fecundity (eggs)", cols = "black", pch = 16, add = "false")
trait.plots.poly("lifetime.eggs",le.f9.briere_T.summary, f9.data.block.temp, "Fecundity (eggs)", cols = "darkcyan", pch = 16, add = "true")
trait.plots.poly("lifetime.eggs",le.f12.briere_T.summary, f12.data.block.temp, "Fecundity (eggs)", cols = "#4a3fa8", pch = 16, add = "true")
trait.plots.poly("gamma",g.c.quad_T.summary, c.gamma.data, "Gamma", cols = "black", pch = 16, add = "false")
trait.plots.poly("gamma",g.f9.quad_T.summary, f9.gamma.data, "Gamma", cols = "darkcyan", pch = 16, add = "true")
trait.plots.poly("gamma",g.f12.quad_T.summary, f12.gamma.data, "Gamma", cols = "#4a3fa8", pch = 16, add = "true")
par(mfrow = c(1,1))
dev.off()

###########################################################################
###########################################################################
### Supplemental Figure 2 - all thermal responses

# Plot thermal performance curves with data they were fit over
pdf("plots/SIFig2_trait_fits_poly_means.pdf", width = 10.6, height = 10)
par(mfrow = c(4,4))
trait.plots.poly("lifespan",lf.c.quad_T.summary , c.data.block.temp, "Lifespan (days)", cols = "black", pch = 1, add = "false")
mtext("a", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly("bite.rate",br.c.briere_T.summary, c.data.block.temp, "Bite rate (days-1)", cols = "black", pch = 1, add = "false")
mtext("b", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly("lifetime.eggs",le.c.briere_T.summary, c.data.block.temp, "Fecundity (eggs)", cols = "black", pch = 1, add = "false")
mtext("c", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly("gamma",g.c.quad_T.summary, c.gamma.data, "Gamma", cols = "black", pch = 1, add = "false")
mtext("d", side = 3, at = -5, cex = 1.8, line = 2)

trait.plots.poly("lifespan",lf.f9.quad_T.summary , f9.data.block.temp, "Lifespan (days)", cols = "darkcyan", pch = 1, add ="false")
mtext("e", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly("bite.rate",br.f9.briere_T.summary , f9.data.block.temp, "Bite rate(days-1)", cols = "darkcyan", pch = 1, add = "false")
mtext("f", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly("lifetime.eggs",le.f9.briere_T.summary, f9.data.block.temp, "Fecundity (eggs)", cols = "darkcyan", pch = 1, add = "false")
mtext("g", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly("gamma",g.f9.quad_T.summary, f9.gamma.data, "Gamma", cols = "darkcyan", pch = 1, add = "false")
mtext("h", side = 3, at = -5, cex = 1.8, line = 2)

trait.plots.poly("lifespan",lf.f12.quad_T.summary , f12.data.block.temp, "Lifespan (days)", cols = "#4a3fa8", pch = 1, add = "false")
mtext("i", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly("bite.rate",br.f12.briere_T.summary , f12.data.block.temp, "Bite rate(days-1)", cols = "#4a3fa8", pch = 1, add = "false")
mtext("j", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly("lifetime.eggs",le.f12.briere_T.summary, f12.data.block.temp, "Fecundity (eggs)", cols = "#4a3fa8", pch = 1, add = "false")
mtext("k", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly("gamma",g.f12.quad_T.summary, f12.gamma.data, "Gamma", cols = "#4a3fa8", pch = 1, add = "false")
mtext("l", side = 3, at = -5, cex = 1.8, line = 2)

trait.plots.poly("Pea", pea.c.quad_T.summary, pea.MDR.data,"Prob. e2a survival", cols = "#AD343E", pch=1, add = "false" )
mtext("m", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly("MDR", MDR.c.briere_T.summary, pea.MDR.data,"Mosq. dev. rate (days-1)", cols = "#AD343E", pch=1, add = "false" )
mtext("n", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly("PDR", PDR.c.briere_T.summary, bc.PDR.data,"Path. dev. rate(days-1)", cols = "#AD343E", pch=1, add = "false" )
mtext("o", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly("bc", bc.c.quad_T.summary,bc.PDR.data,"Vector competence", cols = "#AD343E", pch=1, add = "false" )
mtext("p", side = 3, at = -5, cex = 1.8, line = 2)
par(mfrow = c(1,1))
dev.off()

####Main Figure 1, requires the code for 03 and 03b was run; and then the beginning of this code rerun
# Plot thermal performance curves with data they were fit over
pdf("plots/Main_Fig1_fits_overlay_poly_means.pdf", width = 8, height = 8)
par(mfrow = c(3,2))
trait.plots.poly("lifespan",lf.c.quad_T.summary , c.data.block.temp, "Lifespan (days)", cols = "black", pch = 16, add = "false")
trait.plots.poly("lifespan",lf.f9.quad_T.summary , f9.data.block.temp, "Lifespan (days)", cols = "darkcyan", pch = 16, add ="true")
trait.plots.poly("lifespan",lf.f12.quad_T.summary , f12.data.block.temp, "Lifespan (days)", cols = "#4a3fa8", pch = 16, add = "true")
mtext("a", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly("bite.rate",br.c.briere_T.summary, c.data.block.temp, expression(paste("Bite rate (days"^"-1", ")")), cols = "black", pch = 16, add = "false")
trait.plots.poly("bite.rate",br.f9.briere_T.summary , f9.data.block.temp, expression(paste("Bite rate (days"^"-1", ")")), cols = "darkcyan", pch = 16, add = "true")
trait.plots.poly("bite.rate",br.f12.briere_T.summary , f12.data.block.temp, expression(paste("Bite rate (days"^"-1", ")")), cols = "#4a3fa8", pch = 16, add = "true")
mtext("b", side = 3, at = -5, cex = 1.8, line = 2)
trait.plots.poly("lifetime.eggs",le.c.briere_T.summary, c.data.block.temp, "Lifetime egg production", cols = "black", pch = 16, add = "false")
trait.plots.poly("lifetime.eggs",le.f9.briere_T.summary, f9.data.block.temp, "Lifetime egg production", cols = "darkcyan", pch = 16, add = "true")
trait.plots.poly("lifetime.eggs",le.f12.briere_T.summary, f12.data.block.temp, "Lifetime egg production", cols = "#4a3fa8", pch = 16, add = "true")
mtext("c", side = 3, at = -5, cex = 1.8, line = 2)
plot.R0.poly(df = R0.c.uni.out, cols = "black", add = "false")
plot.R0.poly(df = R0.dtr9.uni.out,  cols = "darkcyan", add ="true")
plot.R0.poly(df = R0.dtr12.uni.out, cols = "#4a3fa8", add ="true")
mtext("d", side = 3, at = -5, cex = 1.8, line = 2)
par(mfrow = c(1,1))
dev.off()

####Supplemental Figure X: Density plot, requires the code for 03 and 03b was run; and then the beginning of this code rerun
# Plot thermal performance curves with data they were fit over
##########
###### 4. Function to calculate distributions of T0, Tmax, and peak R0
##########
# Arguments:
#   input       data frame with posterior distributions of R0 (or other quantity) predicted over a temperature gradient
#           input structure [:,:]
#   temp.list   the temperature gradient itself (i.e., 'Temp.xs')

calcThreshPeakDists = function(input, temp.list) {
  
  # Create output dataframe
  output.df <- data.frame("peak" = numeric(nrow(input)), "T0" = numeric(nrow(input)), "Tmax" = numeric(nrow(input)))
  
  for (i in 1:nrow(input)) { # loop through each row of the input (MCMC step)
    
    output.df$peak[i] <- temp.list[which.max(input[i, ])] # Calculate index of R0 peak & store corresponding temperature
    
    # Makes lists of T0 and Tmax
    index.list <- which(input[i, ] > 0) # Create vector of list of indices where R0 > 1
    length.index.list <- length(index.list)
    output.df$T0[i] <- temp.list[index.list[1] - 1] # Store T0 (index prior to first value in index.list)
    output.df$Tmax[i] <- temp.list[index.list[length.index.list] + 1] # Store Tmax (index after to last value in index.list)
  }
  return(output.df) # return
}

#Make bayesian models ready for use with TPdits functions
#Constant temperature fits
lf.c.quad_T.preds <- trait.trajs.quad(lf.c.quad_T, Temp.xs, summary = FALSE)
lf.c.TPdists <- calcThreshPeakDists(t(lf.c.quad_T.preds), Temp.xs)
br.c.briere_T.preds <- trait.trajs.briere(br.c.briere_T, Temp.xs, summary = FALSE)
br.c.TPdists<- calcThreshPeakDists(t(br.c.briere_T.preds), Temp.xs)
le.c.briere_T.preds <- trait.trajs.briere(le.c.briere_T, Temp.xs, summary = FALSE)
le.c.TPdists <- calcThreshPeakDists(t(le.c.briere_T.preds), Temp.xs)
g.c.quad_T.preds<- trait.trajs.quad(gamma.c.quad_T,Temp.xs, summary = FALSE)
g.c.TPdists <- calcThreshPeakDists(t(g.c.quad_T.preds), Temp.xs)

##temp fluctuation dtr 9
lf.f9.quad_T.preds <- trait.trajs.quad(lf.f9.quad_T, Temp.xs, summary = FALSE)
lf.f9.TPdists <- calcThreshPeakDists(t(lf.f9.quad_T.preds), Temp.xs)
br.f9.briere_T.preds <- trait.trajs.briere(br.f9.briere_T, Temp.xs, summary = FALSE)
br.f9.TPdists <- calcThreshPeakDists(t(br.f9.briere_T.preds), Temp.xs)
le.f9.briere_T.preds <- trait.trajs.briere(le.f9.briere_T, Temp.xs, summary = FALSE)
le.f9.TPdists <- calcThreshPeakDists(t(le.f9.briere_T.preds), Temp.xs)
g.f9.quad_T.preds <- trait.trajs.quad(gamma.f9.quad_T,Temp.xs, summary = FALSE)
g.f9.TPdists <- calcThreshPeakDists(t(g.f9.quad_T.preds), Temp.xs)

##temp fluctuation dtr 12
lf.f12.quad_T.preds <- trait.trajs.quad(lf.f12.quad_T, Temp.xs, summary = FALSE)
lf.f12.TPdists <- calcThreshPeakDists(t(lf.f12.quad_T.preds), Temp.xs)
br.f12.briere_T.preds <- trait.trajs.briere(br.f12.briere_T,Temp.xs, summary = FALSE)
br.f12.TPdists <- calcThreshPeakDists(t(br.f12.briere_T.preds), Temp.xs)
le.f12.briere_T.preds <- trait.trajs.briere(le.f12.briere_T, Temp.xs, summary = FALSE)
le.f12.TPdists <- calcThreshPeakDists(t(le.f12.briere_T.preds), Temp.xs)
g.f12.quad_T.preds <- trait.trajs.quad(gamma.f12.quad_T,Temp.xs, summary = FALSE)
g.f12.TPdists <- calcThreshPeakDists(t(g.f12.quad_T.preds), Temp.xs)

break.n <- seq(0,50,0.4)
#Tmin - lifespan
p1 <- hist(lf.c.TPdists$T0, breaks = break.n,freq= F)
p2 <- hist(lf.f9.TPdists$T0, breaks = break.n,freq= F)
p3 <- hist(lf.f12.TPdists$T0, breaks = break.n,freq= F)
#Topt - lifespan
p4 <- hist(lf.c.TPdists$peak, breaks = break.n,freq= F)
p5 <- hist(lf.f9.TPdists$peak, breaks = break.n,freq= F)
p6 <- hist(lf.f12.TPdists$peak, breaks = break.n,freq= F)
#Tmax - lifespan
p7 <- hist(lf.c.TPdists$Tmax, breaks = break.n,freq= F)
p8 <- hist(lf.f9.TPdists$Tmax, breaks = break.n,freq= F)
p9 <- hist(lf.f12.TPdists$Tmax, breaks = break.n,freq= F)
#Tmin - bite rate
p10 <- hist(br.c.TPdists$T0, breaks = break.n,freq= F)
p11 <- hist(br.f9.TPdists$T0, breaks = break.n,freq= F)
p12 <- hist(br.f12.TPdists$T0, breaks = break.n,freq= F)
#Topt - bite rate
p13 <- hist(br.c.TPdists$peak, breaks = break.n,freq= F)
p14 <- hist(br.f9.TPdists$peak, breaks = break.n,freq= F)
p15 <- hist(br.f12.TPdists$peak, breaks = break.n,freq= F)
#Tmax - bite rate
p16 <- hist(br.c.TPdists$Tmax, breaks = break.n,freq= F)
p17 <- hist(br.f9.TPdists$Tmax, breaks = break.n,freq= F)
p18 <- hist(br.f12.TPdists$Tmax, breaks = break.n,freq= F)
#Tmin - lifetime eggs
p20 <- hist(le.c.TPdists$T0, breaks = break.n,freq= F)
p21 <- hist(le.f9.TPdists$T0, breaks = break.n,freq= F)
p22 <- hist(le.f12.TPdists$T0, breaks = break.n,freq= F)
#Topt - lifetime eggs
p23 <- hist(le.c.TPdists$peak, breaks = break.n,freq= F)
p24 <- hist(le.f9.TPdists$peak, breaks = break.n,freq= F)
p25 <- hist(le.f12.TPdists$peak, breaks = break.n,freq= F)
#Tmax - lifetime eggs
p26 <- hist(le.c.TPdists$Tmax, breaks = break.n,freq= F)
p27 <- hist(le.f9.TPdists$Tmax, breaks = break.n,freq= F)
p28 <- hist(le.f12.TPdists$Tmax, breaks = break.n,freq= F)

#Tmin - ST
p30 <- hist(R0.c.TPdists$T0, breaks = break.n, freq = F)
p31 <- hist(R0.dtr9.TPdists$T0 , breaks = break.n, freq = F)
p32 <- hist(R0.dtr12.TPdists$T0, breaks = break.n, freq = F)
#Topt - ST
p33 <- hist(R0.c.TPdists$peak, breaks = break.n, freq = F)
p34 <- hist(R0.dtr9.TPdists$peak , breaks = break.n, freq = F)
p35 <- hist(R0.dtr12.TPdists$peak, breaks = break.n, freq = F)
#Tmax - ST
p36 <- hist(R0.c.TPdists$Tmax, breaks = break.n, freq = F)
p37 <- hist(R0.dtr9.TPdists$Tmax, breaks = break.n, freq = F)
p38 <- hist(R0.dtr12.TPdists$Tmax, breaks = break.n, freq = F)


###############################Ok make a seperate plot for each trait
pdf("plots/Thresholds_density_lf_a_B_S_experimental.pdf", width = 12, height = 12) #width = 8, height = 16) #alter dims to fit onto a single page// 4 by 8/// 2 by 4// 1 x 2//
par(mar=c(5,6,5,1))
par(mfrow = c(4,3))
#lifespan
plot(p1, col =scales::alpha('grey43',0.4),freq= F,border=scales::alpha('black',0.6), xlim= c(0,15),ylim = c(0,0.3),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, main = expression(italic(paste("T"[min])- lf)))
plot(p2, col= scales::alpha('darkcyan',0.4),freq= F,border=scales::alpha('darkcyan',0.6), ylab = "Density",add = T)
plot(p3, col= scales::alpha('#4a3fa8',0.4),freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
legend('topright', legend = c("DTR 0C", "DTR 9C", "DTR 12C"), lty = 1, col = c("black", "darkcyan","#4a3fa8"), bty = 'n')
mtext("a", side = 3, at = -3, cex = 1.8, line = 2)
plot(p4, col =scales::alpha('grey43',0.4),freq= F,border=scales::alpha('black',0.6), xlim= c(12,24),ylim = c(0,0.5),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, main = expression(italic(paste("T"[opt])-lf)))
plot(p5, col= scales::alpha('darkcyan',0.4),freq= F,border=scales::alpha('darkcyan',0.6), ylab = "Density",add = T)
plot(p6, col= scales::alpha('#4a3fa8',0.4),freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
#legend('topright', legend = c("DTR 0C", "DTR 9C", "DTR 12C"), lty = 1, col = c("black", "darkcyan","#4a3fa8"), bty = 'n')
plot(p7, col =scales::alpha('grey43',0.4),freq= F,border=scales::alpha('black',0.6), xlim= c(24,40),ylim = c(0,0.4),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, main = expression(italic(paste("T"[max])-lf)))
plot(p8, col= scales::alpha('darkcyan',0.4),freq= F,border=scales::alpha('darkcyan',0.6), ylab = "Density",add = T)
plot(p9, col= scales::alpha('#4a3fa8',0.4),freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
#legend('topright', legend = c("DTR 0C", "DTR 9C", "DTR 12C"), lty = 1, col = c("black", "darkcyan","#4a3fa8"), bty = 'n')

#bite rate
plot(p10, col =scales::alpha('grey43',0.4),freq= F,border=scales::alpha('black',0.6), xlim= c(0,8),ylim = c(0,.8),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, main = expression(italic(paste("T"[min])-a)))
plot(p11, col= scales::alpha('darkcyan',0.4),freq= F,border=scales::alpha('darkcyan',0.6), ylab = "Density",add = T)
plot(p12, col= scales::alpha('#4a3fa8',0.4),freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
legend('topright', legend = c("DTR 0C", "DTR 9C", "DTR 12C"), lty = 1, col = c("black", "darkcyan","#4a3fa8"), bty = 'n')
mtext("b", side = 3, at = -1.5, cex = 1.8, line = 2)
plot(p13, col =scales::alpha('grey43',0.4),freq= F,border=scales::alpha('black',0.6), xlim= c(26,36),ylim = c(0,0.6),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, main = expression(italic(paste("T"[opt])-a)))
plot(p14, col= scales::alpha('darkcyan', 0.4), freq= F,border=scales::alpha('darkcyan',0.6), ylab = "Density",add = T)
plot(p15, col= scales::alpha('#4a3fa8',0.4),freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
#legend('topright', legend = c("DTR 0C", "DTR 9C", "DTR 12C"), lty = 1, col = c("black", "darkcyan","#4a3fa8"), bty = 'n')
plot(p16, col =scales::alpha('grey43',0.4),freq= F,border=scales::alpha('black',0.6), xlim= c(32,46),ylim = c(0,0.4),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, main = expression(italic(paste("T"[max])-a)))
plot(p17, col= scales::alpha('darkcyan',0.4),freq= F,border=scales::alpha('darkcyan',0.6), ylab = "Density",add = T)
plot(p18, col= scales::alpha('#4a3fa8',0.4),freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
#legend('topright', legend = c("DTR 0C", "DTR 9C", "DTR 12C"), lty = 1, col = c("black", "darkcyan","#4a3fa8"), bty = 'n')

#lifetime eggs
plot(p20, col =scales::alpha('grey43',0.4),freq= F,border=scales::alpha('black',0.6), xlim= c(0,20),ylim = c(0,.12),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, main = expression(italic(paste("T"[min])-B)))
plot(p21, col= scales::alpha('darkcyan',0.4),freq= F,border=scales::alpha('darkcyan',0.6), ylab = "Density",add = T)
plot(p22, col= scales::alpha('#4a3fa8',0.4),freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
legend('topleft', legend = c("DTR 0C", "DTR 9C", "DTR 12C"), lty = 1, col = c("black", "darkcyan","#4a3fa8"), bty = 'n')
mtext("c", side = 3, at = -4, cex = 1.8, line = 2)
plot(p23, col =scales::alpha('grey43',0.4),freq= F,border=scales::alpha('black',0.6), xlim= c(20,34),ylim = c(0,0.5),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, main = expression(italic(paste("T"[opt])-B)))
plot(p24, col= scales::alpha('darkcyan', 0.4), freq= F,border=scales::alpha('darkcyan',0.6), ylab = "Density",add = T)
plot(p25, col= scales::alpha('#4a3fa8',0.4),freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
#legend('topright', legend = c("DTR 0C", "DTR 9C", "DTR 12C"), lty = 1, col = c("black", "darkcyan","#4a3fa8"), bty = 'n')
plot(p26, col =scales::alpha('grey43',0.4),freq= F,border=scales::alpha('black',0.6), xlim= c(25,42),ylim = c(0,0.8),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, main = expression(italic(paste("T"[max])-B)))
plot(p27, col= scales::alpha('darkcyan',0.4),freq= F,border=scales::alpha('darkcyan',0.6), ylab = "Density",add = T)
plot(p28, col= scales::alpha('#4a3fa8',0.4),freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
#legend('topright', legend = c("DTR 0C", "DTR 9C", "DTR 12C"), lty = 1, col = c("black", "darkcyan","#4a3fa8"), bty = 'n')

#thermal suitability 
plot(p30, col =scales::alpha('grey43',0.4),freq= F,border=scales::alpha('black',0.6), xlim= c(12,22),ylim = c(0,1),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, main = expression(italic(paste("T"[min])-S)))
plot(p31, col= scales::alpha('darkcyan',0.4),freq= F,border=scales::alpha('darkcyan',0.6), ylab = "Density",add = T)
plot(p32, col= scales::alpha('#4a3fa8',0.4),freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
legend('topright', legend = c("DTR 0C", "DTR 9C", "DTR 12C"), lty = 1, col = c("black", "darkcyan","#4a3fa8"), bty = 'n')
mtext("d", side = 3, at = 10, cex = 1.8, line = 2)
plot(p33, col =scales::alpha('grey43',0.4),freq= F,border=scales::alpha('black',0.6), xlim= c(20,30),ylim = c(0,1),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, main = expression(italic(paste("T"[opt])-S)))
plot(p34, col= scales::alpha('darkcyan', 0.4), freq= F,border=scales::alpha('darkcyan',0.6), ylab = "Density",add = T)
plot(p35, col= scales::alpha('#4a3fa8',0.4),freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
#legend('topright', legend = c("DTR 0C", "DTR 9C", "DTR 12C"), lty = 1, col = c("black", "darkcyan","#4a3fa8"), bty = 'n')
plot(p36, col =scales::alpha('grey43',0.4),freq= F,border=scales::alpha('black',0.6), xlim= c(24,38),ylim = c(0,0.6),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, main = expression(italic(paste("T"[max])-S)))
plot(p37, col= scales::alpha('darkcyan',0.4),freq= F,border=scales::alpha('darkcyan',0.6), ylab = "Density",add = T)
plot(p38, col= scales::alpha('#4a3fa8',0.4),freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
#legend('topright', legend = c("DTR 0C", "DTR 9C", "DTR 12C"), lty = 1, col = c("black", "darkcyan","#4a3fa8"), bty = 'n')

par(mfrow = c(1,1))
dev.off()

##############################################################################
##############Calculate and export the summary of the TPdists for each trait
# Output used in main table 2
#########################################################################
# Function to summarize the median and 95% CI
summary.fn = function(df){
  med = median(df, na.rm = T)
  int = c(HPDinterval(mcmc(df)))
  out = c(med, int)
  names(out) = c("median", "lowerCI", "upperCI")
  out
}

#Apply summary.fn to each TPdists df
lf.c.interval <- apply(lf.c.TPdists, MARGIN =2, summary.fn)
lf.dtr9.interval <- apply(lf.f9.TPdists, MARGIN =2, summary.fn)
lf.dtr12.interval <- apply(lf.f12.TPdists, MARGIN =2, summary.fn)
br.c.interval <- apply(br.c.TPdists, MARGIN =2, summary.fn)
br.dtr9.interval <- apply(br.f9.TPdists, MARGIN =2, summary.fn)
br.dtr12.interval <- apply(br.f12.TPdists, MARGIN =2, summary.fn)
le.c.interval <- apply(le.c.TPdists, MARGIN =2, summary.fn)
le.dtr9.interval <- apply(le.f9.TPdists, MARGIN =2, summary.fn)
le.dtr12.interval <- apply(le.f12.TPdists, MARGIN =2, summary.fn)
s.c.interval <- apply(R0.c.TPdists, MARGIN =2, summary.fn)
s.dtr9.interval <- apply(R0.dtr9.TPdists, MARGIN =2, summary.fn)
s.dtr12.interval <- apply(R0.dtr12.TPdists, MARGIN =2, summary.fn)

experimental.trait.TPdists <- list(lf.c.interval= lf.c.interval,
                                   lf.dtr9.interval =lf.dtr9.interval,
                                   lf.dtr12.interval=lf.dtr12.interval,
                                   br.c.interval=br.c.interval,
                                   br.dtr9.interval= br.dtr9.interval,
                                   br.dtr12.interval =br.dtr12.interval,
                                   le.c.interval= le.c.interval,
                                   le.dtr9.interval=le.dtr9.interval,
                                   le.dtr12.interval= le.dtr12.interval,
                                   s.c.interval=s.c.interval,
                                   s.dtr9.interval = s.dtr9.interval,
                                   s.dtr12.interval = s.dtr12.interval)
#sink("exp_trait_thresholds_lf_a_B_s.txt") #file name
#print(experimental.trait.TPdists) #object to be exported as text
#sink() #exporting to source file location
