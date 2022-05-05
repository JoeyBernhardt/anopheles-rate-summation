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



#load RS estimates
all.RSestimates <- read.csv("data/traits.integrated.csv")
all.RSestimates <- read.csv("data/traits.integrated.CTmaxcutoff.csv")
#Make sure the 03_Calc_R0_constant.R code file was ran and the files are still accessible. (need the matrix object for R0.constant)

#Calculate R0.estimates (using constant data for pea, mdr, and bc)
#R0.full = function(a, bc, lf, y, B, pEA, MDR)
R0.dtr9.estimates <- R0.full(subset(all.RSestimates$br, all.RSestimates$dtr == ".dtr9"), bc.c.quad_T.summary$mean,subset(all.RSestimates$lf, all.RSestimates$dtr == ".dtr9"), subset(all.RSestimates$g, all.RSestimates$dtr == ".dtr9"), subset(all.RSestimates$le, all.RSestimates$dtr == ".dtr9"), pea.c.quad_T.summary$mean, MDR.c.briere_T.summary$mean)
R0.dtr9.estimates<- data.frame(temp = Temp.xs, mean = R0.dtr9.estimates)

R0.dtr12.estimates <- R0.full(subset(all.RSestimates$br, all.RSestimates$dtr == ".dtr12"), bc.c.quad_T.summary$mean,subset(all.RSestimates$lf, all.RSestimates$dtr == ".dtr12"), subset(all.RSestimates$g, all.RSestimates$dtr == ".dtr12"), subset(all.RSestimates$le, all.RSestimates$dtr == ".dtr12"), pea.c.quad_T.summary$mean, MDR.c.briere_T.summary$mean)
R0.dtr12.estimates<- data.frame(temp = Temp.xs, mean = R0.dtr12.estimates)

#export results
#write.csv(R0.dtr9.estimates, "data/R0.dtr9_only_B_a_lf.csv")
#write.csv(R0.dtr12.estimates, "data/R0.dtr12_only_B_a_lf.csv")

###################################################
###############################################
#Calculate R0.estimates (using RS data for all traits)
R0.dtr9.all.estimates <- R0.full(subset(all.RSestimates$br, all.RSestimates$dtr == ".dtr9"), subset(all.RSestimates$bc, all.RSestimates$dtr == ".dtr9"),subset(all.RSestimates$lf, all.RSestimates$dtr == ".dtr9"), subset(all.RSestimates$g, all.RSestimates$dtr == ".dtr9"), subset(all.RSestimates$le, all.RSestimates$dtr == ".dtr9"), subset(all.RSestimates$pea, all.RSestimates$dtr == ".dtr9"), subset(all.RSestimates$mdr, all.RSestimates$dtr == ".dtr9"))
R0.dtr9.all.estimates <- data.frame(temp = Temp.xs, mean = R0.dtr9.all.estimates)

R0.dtr12.all.estimates <- R0.full(subset(all.RSestimates$br, all.RSestimates$dtr == ".dtr12"), subset(all.RSestimates$bc, all.RSestimates$dtr == ".dtr12"),subset(all.RSestimates$lf, all.RSestimates$dtr == ".dtr12"), subset(all.RSestimates$g, all.RSestimates$dtr == ".dtr12"), subset(all.RSestimates$le, all.RSestimates$dtr == ".dtr12"), subset(all.RSestimates$pea, all.RSestimates$dtr == ".dtr12"), subset(all.RSestimates$mdr, all.RSestimates$dtr == ".dtr12"))
R0.dtr12.all.estimates <- data.frame(temp = Temp.xs, mean = R0.dtr12.all.estimates)



#############
####### 2. Plot thermal performance curves
#############


##############################################
##First I need to redefine the mapping function
###########################################
pea.MDR.data$Treatment <- pea.MDR.data$temp
bc.PDR.data$Treatment <- bc.PDR.data$temp
bc.PDR.data$PDR <- 1/bc.PDR.data$EIP50
c.gamma.data$Treatment <- c.gamma.data$temp
f9.gamma.data$Treatment <- f9.gamma.data$temp
f12.gamma.data$Treatment <- f12.gamma.data$temp


#####################Alternative plotting function
#credible intervals as shaded lines instead of dashed lines

trait.plots.poly = function(traitname, summary,df, add = "false", cols = "black", pch = 1, labs = ""){
  temp = df$Treatment
  traitplot = (df[ , which(colnames(df)==traitname)])
  if (add == "false") {
    par(mar=c(5,5,5,1))
    plot(traitplot ~ temp, lty = 1,xlim = c(0,45), xlab = expression(paste("Temperature (",degree,"C)")), ylab = labs, cex.axis = 1.4, cex.lab = 1.4, ylim = c(0, max(traitplot, na.rm = T)*1.2), cex.main = 1.5, col = cols, pch = pch)
    polygon(c(summary$temp, rev(summary$temp)), c(summary$upper, rev(summary$lower)), col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  } else {
    polygon(c(summary$temp, rev(summary$temp)), c(summary$upper, rev(summary$lower)), col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  }
  points(traitplot~temp, col = cols, pch = pch)
  lines(summary$temp, summary$mean, lty = 1, lwd = 2, col =  adjustcolor(cols, alpha.f = 1))
}



# Alternative function that uses shaded polygons instead of dashed lines
#Makes the mean max at 1.
plot.R0.poly = function(temp = Temp.xs, df, modname = "", add = "false", cols = "black"){
  if (add == "false") {
    par(mar=c(5,5,5,1))
    plot(temp, df$upperCI/max(df$mean, na.rm = T), lty = 1, xlab = expression(paste("Temperature (",degree,"C)")), ylab = expression(S[T]), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, yaxt = 'n', col = "white")
    polygon(c(temp, rev(temp)), c(df$upperCI/max(df$mean, na.rm = T), rev(df$lowerCI/max(df$mean, na.rm = T))), col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  } else {
    polygon(c(temp, rev(temp)), c(df$upperCI/max(df$mean, na.rm = T), rev(df$lowerCI/max(df$mean, na.rm = T))), col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  }
  
  lines(temp, df$mean/max(df$mean, na.rm = T), lty = 1, lwd = 2, col = cols)
}


##########################
########plotting function to add RS estimates (assume 'add' always equals true)
trait.estimates = function(traitname, df, cols = "black", linetype = "dotted"){
  temp = df$mean.temp
  traitplot = (df[ , which(colnames(df)==traitname)])
  #points(traitplot~temp, col = cols, pch = 1)
  lines(traitplot~temp, lty = linetype, lwd = 2, col = cols)
}

##########################
########plotting function to add RS estimates (assume 'add' always equals true)
R0.estimates = function(traitname, df, cols = "black", linetype = "dotted"){
  temp = df$mean.temp
  traitplot = (df[ , which(colnames(df)==traitname)])
  #points(traitplot~temp, col = cols, pch = 1)
  lines(traitplot/max(traitplot, na.rm = TRUE)~temp, lty = linetype, lwd = 2, col = cols)
}

########plotting function to add RS estimates (assume 'add' always equals true)
####For the estimates when only the B, lf, gamma, and a fluc were added.
R0.estimates.partial = function(df, cols = "black", linetype = "dotted"){
  temp = df$temp
  traitplot = (df[ , which(colnames(df)=="mean")])
  #points(traitplot~temp, col = cols, pch = 1)
  lines(traitplot/max(traitplot, na.rm = TRUE)~temp, lty = linetype, lwd = 2, col = cols)
}

# Plot thermal performance curves with data they were fit over
# Plot size A
pdf("plots/RS_overlay_squareDims.pdf", width = 8, height = 10)
par(mfrow = c(4,3))

trait.plots.poly("lifespan",lf.f9.quad_T.summary , f9.data.block.temp, "Lifespan (days)", cols = "darkcyan", pch = 1, add ="false")
trait.plots.poly("lifespan",lf.f12.quad_T.summary , f12.data.block.temp, "Lifespan (days)", cols = "#4a3fa8", pch = 1, add = "true")
trait.estimates("lf", subset(all.RSestimates, dtr==".dtr9"), cols = "darkcyan", linetype = "dashed")
trait.estimates("lf", subset(all.RSestimates, dtr==".dtr12"), cols = "#4a3fa8", linetype = "dashed")

trait.plots.poly("bite.rate",br.f9.briere_T.summary , f9.data.block.temp, expression(paste("Bite rate (days"^"-1", ")")), cols = "darkcyan", pch = 1, add = "false")
trait.plots.poly("bite.rate",br.f12.briere_T.summary , f12.data.block.temp, expression(paste("Bite rate (days"^"-1", ")")), cols = "#4a3fa8", pch = 1, add = "true")
trait.estimates("br", subset(all.RSestimates, dtr==".dtr9"), cols = "darkcyan", linetype = "dashed")
trait.estimates("br", subset(all.RSestimates, dtr==".dtr12"), cols = "#4a3fa8", linetype = "dashed")

trait.plots.poly("lifetime.eggs",le.f9.briere_T.summary, f9.data.block.temp, "Fecundity (eggs)", cols = "darkcyan", pch = 1, add = "false")
trait.plots.poly("lifetime.eggs",le.f12.briere_T.summary, f12.data.block.temp, "Fecundity (eggs)", cols = "#4a3fa8", pch = 1, add = "true")
trait.estimates("le", subset(all.RSestimates, dtr==".dtr9"), cols = "darkcyan", linetype = "dashed")
trait.estimates("le", subset(all.RSestimates, dtr==".dtr12"), cols = "#4a3fa8", linetype = "dashed")

trait.plots.poly("gamma",g.f9.quad_T.summary, f9.gamma.data, "Gamma", cols = "darkcyan", pch = 1, add = "false")
trait.plots.poly("gamma",g.f12.quad_T.summary, f12.gamma.data, "Gamma", cols = "#4a3fa8", pch = 1, add = "true")
trait.estimates("g", subset(all.RSestimates, dtr==".dtr9"), cols = "darkcyan", linetype = "dashed")
trait.estimates("g", subset(all.RSestimates, dtr==".dtr12"), cols = "#4a3fa8", linetype = "dashed")

plot.R0.poly(df = R0.dtr9.uni.out, cols = "darkcyan")
plot.R0.poly(df = R0.dtr12.uni.out, cols = "#4a3fa8", add ="true")
R0.estimates.partial(df = R0.dtr9.estimates,  cols = "darkcyan", linetype = "dashed")
R0.estimates.partial(df = R0.dtr12.estimates,  cols = "#4a3fa8", linetype = "dashed")
par(mfrow = c(1,1))
dev.off()


####Overlay by trait figure
# Plot thermal performance curves with data they were fit over
pdf("plots/RS_overlay_rectDims.pdf", width = 4, height = 8)
par(mfrow = c(3,1))
trait.plots.poly("lifespan",lf.f9.quad_T.summary , f9.data.block.temp, "Lifespan (days)", cols = "darkcyan", pch = 16, add ="false")
trait.plots.poly("lifespan",lf.f12.quad_T.summary , f12.data.block.temp, "Lifespan (days)", cols = "#4a3fa8", pch = 16, add = "true")
trait.estimates("lf", subset(all.RSestimates, dtr==".dtr9"), cols = "darkcyan", linetype = "dashed")
trait.estimates("lf", subset(all.RSestimates, dtr==".dtr12"), cols = "#4a3fa8", linetype = "dashed")

trait.plots.poly("bite.rate",br.f9.briere_T.summary , f9.data.block.temp, expression(paste("Bite rate (days"^"-1", ")")), cols = "darkcyan", pch = 16, add = "false")
trait.plots.poly("bite.rate",br.f12.briere_T.summary , f12.data.block.temp, expression(paste("Bite rate (days"^"-1", ")")), cols = "#4a3fa8", pch = 16, add = "true")
trait.estimates("br", subset(all.RSestimates, dtr==".dtr9"), cols = "darkcyan", linetype = "dashed")
trait.estimates("br", subset(all.RSestimates, dtr==".dtr12"), cols = "#4a3fa8", linetype = "dashed")


trait.plots.poly("lifetime.eggs",le.f9.briere_T.summary, f9.data.block.temp, "Fecundity (eggs)", cols = "darkcyan", pch = 16, add = "false")
trait.plots.poly("lifetime.eggs",le.f12.briere_T.summary, f12.data.block.temp, "Fecundity (eggs)", cols = "#4a3fa8", pch = 16, add = "true")
trait.estimates("le", subset(all.RSestimates, dtr==".dtr9"), cols = "darkcyan", linetype = "dashed")
trait.estimates("le", subset(all.RSestimates, dtr==".dtr12"), cols = "#4a3fa8", linetype = "dashed")


trait.plots.poly("gamma",g.f9.quad_T.summary, f9.gamma.data, "Gamma", cols = "darkcyan", pch = 16, add = "false")
trait.plots.poly("gamma",g.f12.quad_T.summary, f12.gamma.data, "Gamma", cols = "#4a3fa8", pch = 16, add = "true")
trait.estimates("g", subset(all.RSestimates, dtr==".dtr9"), cols = "darkcyan", linetype = "dashed")
trait.estimates("g", subset(all.RSestimates, dtr==".dtr12"), cols = "#4a3fa8", linetype = "dashed")

plot.R0.poly(df = R0.dtr9.uni.out, cols = "darkcyan")
plot.R0.poly(df = R0.dtr12.uni.out, cols = "#4a3fa8", add ="true")
R0.estimates.partial(df = R0.dtr9.estimates,  cols = "darkcyan", linetype = "dashed")
R0.estimates.partial(df = R0.dtr12.estimates,  cols = "#4a3fa8", linetype = "dashed")

par(mfrow = c(1,1))
dev.off()


##############single figure overlay of the three ways to look at ST
#####1. constant temp R0
#####2. RS on all traits of R0
#####3. RS directly on R0

pdf("plots/RS_3methods_overlay_rectDims.pdf", width = 4, height = 8)
par(mfrow = c(3,1))
#####1. constant temp R0
plot.R0.poly(df = R0.c.uni.out, cols = "black")
#####2. RS on all traits of R0
R0.estimates.partial(df = R0.dtr9.all.estimates,  cols = "darkcyan", linetype = "dashed")
R0.estimates.partial(df = R0.dtr9.all.estimates,  cols = "#4a3fa8", linetype = "dashed")


##3. RS directly on R0
R0.estimates("R0", subset(all.RSestimates, dtr==".dtr9"), cols = "darkcyan", linetype = "dotted")
R0.estimates("R0", subset(all.RSestimates, dtr==".dtr12"), cols = "#4a3fa8", linetype = "dotted")
par(mfrow = c(1,1))
dev.off()


######################ST overlay
pdf("plots/ST_overlay_poly_means.pdf", width = 4, height = 8)
par(mfrow = c(3,1))
plot.R0.poly(df = R0.c.uni.out, cols = "black")
plot.R0.poly(df = R0.dtr9.uni.out, cols = "darkcyan", add = "true")
plot.R0.poly(df = R0.dtr12.uni.out, cols = "#4a3fa8", add ="true")
par(mfrow = c(1,1))
dev.off()
