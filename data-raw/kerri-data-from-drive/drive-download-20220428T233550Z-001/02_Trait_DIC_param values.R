###set up working directory
#Set working directory
mainDir = "C:/Users/Kerri/Desktop/Chapter2 InProgress"
setwd(mainDir)

#Constant temperature fits
load("saved posteriors/constant_lf_c.quad_T.uniform.Rdata")
lf.c.quad.DIC <- lf.c.quad_T$BUGSoutput$DIC
lf.c.quad.params <- lf.c.quad_T$BUGSoutput$summary[,c(1,3,5,7)] #takes only the mean and 95% interval
load("saved posteriors/constant_lf.c.briere_T.uniform.Rdata")
lf.c.briere.DIC <- lf.c.briere_T$BUGSoutput$DIC
lf.c.briere.params <- lf.c.briere_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/constant_br.c.briere_T.uniform.Rdata")
br.c.briere.DIC <- br.c.briere_T$BUGSoutput$DIC
br.c.briere.params <- br.c.briere_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/constant_br_c.quad_T.uniform.Rdata")
br.c.quad.DIC <- br.c.quad_T$BUGSoutput$DIC
br.c.quad.params <- br.c.quad_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/constant_le.c.briere_T.uniform.Rdata")
le.c.briere.DIC <- le.c.briere_T$BUGSoutput$DIC
le.c.briere.params <- le.c.briere_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/constant_le_c.quad_T.uniform.Rdata")
le.c.quad.DIC <- le.c.quad_T$BUGSoutput$DIC
le.c.quad.params <- le.c.quad_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/constant_gamma_c.quad_T.uniform.Rdata")
g.c.quad.DIC <- gamma.c.quad_T$BUGSoutput$DIC
g.c.quad.params <- gamma.c.quad_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/constant_gamma.c.briere_T.uniform.Rdata")
g.c.briere.DIC <- gamma.c.briere_T$BUGSoutput$DIC
g.c.briere.params <- gamma.c.briere_T$BUGSoutput$summary[,c(1,3,5,7)]

##temp fluctuation dtr 9
load("saved posteriors/dtr_lf_f9.quad_T.uniform.Rdata")
lf.f9.quad.DIC <- lf.f9.quad_T$BUGSoutput$DIC
lf.f9.quad.params <- lf.f9.quad_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/dtr_lf.f9.briere_T.uniform.Rdata")
lf.f9.briere.DIC <- lf.f9.briere_T$BUGSoutput$DIC
lf.f9.briere.params <- lf.f9.briere_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/dtr_br.f9.briere_T.uniform.Rdata")
br.f9.briere.DIC <- br.f9.briere_T$BUGSoutput$DIC
br.f9.briere.params <- br.f9.briere_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/dtr_br_f9.quad_T.uniform.Rdata")
br.f9.quad.DIC <- br.f9.quad_T$BUGSoutput$DIC
br.f9.quad.params <- br.f9.quad_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/dtr_le.f9.briere_T.uniform.Rdata")
le.f9.briere.DIC <- le.f9.briere_T$BUGSoutput$DIC
le.f9.briere.params <- le.f9.briere_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/dtr_le_f9.quad_T.uniform.Rdata")
le.f9.quad.DIC <- le.f9.quad_T$BUGSoutput$DIC
le.f9.quad.params <- le.f9.quad_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/dtr_gamma_f9.quad_T.uniform.Rdata")
g.f9.quad.DIC <- gamma.f9.quad_T$BUGSoutput$DIC
g.f9.quad.params <- gamma.f9.quad_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/dtr_gamma.f9.briere_T.uniform.Rdata")
g.f9.briere.DIC <- gamma.f9.briere_T$BUGSoutput$DIC
g.f9.briere.params <- gamma.f9.briere_T$BUGSoutput$summary[,c(1,3,5,7)]
#need to fit to the counter structure for these too

##temp fluctuation dtr 12
load("saved posteriors/dtr_lf_f12.quad_T.uniform.Rdata")
lf.f12.quad.DIC <- lf.f12.quad_T$BUGSoutput$DIC
lf.f12.quad.params <- lf.f12.quad_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/dtr_lf.f12.briere_T.uniform.Rdata")
lf.f12.briere.DIC <- lf.f12.briere_T$BUGSoutput$DIC
lf.f12.briere.params <- lf.f12.briere_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/dtr_br.f12.briere_T.uniform.Rdata")
br.f12.briere.DIC <- br.f12.briere_T$BUGSoutput$DIC
br.f12.briere.params <- br.f12.briere_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/dtr_br_f12.quad_T.uniform.Rdata")
br.f12.quad.DIC <- br.f12.quad_T$BUGSoutput$DIC
br.f12.quad.params <- br.f12.quad_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/dtr_le.f12.briere_T.uniform.Rdata")
le.f12.briere.DIC <- le.f12.briere_T$BUGSoutput$DIC
le.f12.briere.params <- le.f12.briere_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/dtr_le_f12.quad_T.uniform.Rdata")
le.f12.quad.DIC <- le.f12.quad_T$BUGSoutput$DIC
le.f12.quad.params <- le.f12.quad_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/dtr_gamma_f12.quad_T.uniform.Rdata")
g.f12.quad.DIC <- gamma.f12.quad_T$BUGSoutput$DIC
g.f12.quad.params <- gamma.f12.quad_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/dtr_gamma_f12.briere_T.uniform.Rdata")
g.f12.briere.DIC <- gamma.f12.briere_T$BUGSoutput$DIC
g.f12.briere.params <- gamma.f12.briere_T$BUGSoutput$summary[,c(1,3,5,7)]
#need to fit to the counter structure for these too

###literature based traits
load("saved posteriors/bc_quad_withT_uniform.Rdata")
bc.c.quad.DIC <- bc.c.quad_T$BUGSoutput$DIC
bc.c.quad.params <- bc.c.quad_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/bc_briere_withT_uniform.Rdata")
bc.c.briere.DIC <- bc.c.briere_T$BUGSoutput$DIC
bc.c.briere.params <- bc.c.briere_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/MDR_briere_withT_uniform.Rdata")
MDR.c.briere.DIC <- MDR.c.briere_T$BUGSoutput$DIC
MDR.c.briere.params <- MDR.c.briere_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/MDR_quad_withT_uniform.Rdata")
MDR.c.quad.DIC <- MDR.c.quad_T$BUGSoutput$DIC
MDR.c.quad.params <- MDR.c.quad_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/PDR_briere_withT_uniform.Rdata")
PDR.c.briere.DIC <- PDR.c.briere_T$BUGSoutput$DIC
PDR.c.briere.params <- PDR.c.briere_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/PDR_quad_withT_uniform.Rdata")
PDR.c.quad.DIC <- PDR.c.quad_T$BUGSoutput$DIC
PDR.c.quad.params <- PDR.c.quad_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/pea_quad_withT_uniform.Rdata")
pea.c.quad.DIC <- pea.c.quad_T$BUGSoutput$DIC
pea.c.quad.params <- pea.c.quad_T$BUGSoutput$summary[,c(1,3,5,7)]
load("saved posteriors/pea_briere_withT_uniform.Rdata")
pea.c.briere.DIC <- pea.c.briere_T$BUGSoutput$DIC
pea.c.briere.params <- pea.c.briere_T$BUGSoutput$summary[,c(1,3,5,7)]


#Place all of the model DICs and parameter values into a single file and export for easy reference
models.DIC_list <- list(lf.c.quad.DIC=lf.c.quad.DIC,
                        lf.c.briere.DIC=lf.c.briere.DIC,
                        br.c.briere.DIC =br.c.briere.DIC,
                        br.c.quad.DIC =br.c.quad.DIC,
                        le.c.briere.DIC=le.c.briere.DIC,
                        le.c.quad.DIC=le.c.quad.DIC,
                        g.c.quad.DIC=g.c.quad.DIC,
                        g.c.briere.DIC =g.c.briere.DIC,
                        bc.c.quad.DIC= bc.c.quad.DIC,
                        bc.c.briere.DIC=bc.c.briere.DIC,
                        MDR.c.briere.DIC=MDR.c.briere.DIC,
                        MDR.c.quad.DIC =MDR.c.quad.DIC,
                        PDR.c.briere.DIC=PDR.c.briere.DIC,
                        PDR.c.quad.DIC=PDR.c.quad.DIC,
                        pea.c.quad.DIC=pea.c.quad.DIC,
                        pea.c.briere.DIC=pea.c.briere.DIC,
                        lf.f9.quad.DIC=lf.f9.quad.DIC,
                        lf.f9.briere.DIC =lf.f9.briere.DIC,
                        br.f9.briere.DIC=br.f9.briere.DIC,
                        br.f9.quad.DIC =br.f9.quad.DIC,
                        le.f9.briere.DIC=le.f9.briere.DIC,
                        le.f9.quad.DIC=le.f9.quad.DIC,
                        g.f9.quad.DIC=g.f9.quad.DIC,
                        g.f9.briere.DIC=g.f9.briere.DIC,
                        lf.f12.quad.DIC=lf.f12.quad.DIC,
                        lf.f12.briere.DIC=lf.f12.briere.DIC,
                        br.f12.briere.DIC=br.f12.briere.DIC,
                        br.f12.quad.DIC=br.f12.quad.DIC,
                        le.f12.briere.DIC=le.f12.briere.DIC,
                        le.f12.quad.DIC=le.f12.quad.DIC,
                        g.f12.quad.DIC=g.f12.quad.DIC,
                        g.f12.briere.DIC=g.f12.briere.DIC)  
sink("tables/model.DIC_list.txt") #file name
print(models.DIC_list) #object to be exported as text
sink() #exporting to source file location

models.params_list <- list(lf.c.quad.params=lf.c.quad.params,
                           lf.c.briere.params=lf.c.briere.params,
                           br.c.briere.params =br.c.briere.params,
                           br.c.quad.params =br.c.quad.params,
                           le.c.briere.params=le.c.briere.params,
                           le.c.quad.params=le.c.quad.params,
                           g.c.quad.params=g.c.quad.params,
                           g.c.briere.params =g.c.briere.params,
                           bc.c.quad.params= bc.c.quad.params,
                           bc.c.briere.params=bc.c.briere.params,
                           MDR.c.briere.params=MDR.c.briere.params,
                           MDR.c.quad.params =MDR.c.quad.params,
                           PDR.c.briere.params=PDR.c.briere.params,
                           PDR.c.quad.params=PDR.c.quad.params,
                           pea.c.quad.params=pea.c.quad.params,
                           pea.c.briere.params=pea.c.briere.params,
                           lf.f9.quad.params=lf.f9.quad.params,
                           lf.f9.briere.params =lf.f9.briere.params,
                           br.f9.briere.params=br.f9.briere.params,
                           br.f9.quad.params =br.f9.quad.params,
                           le.f9.briere.params=le.f9.briere.params,
                           le.f9.quad.params=le.f9.quad.params,
                           g.f9.quad.params=g.f9.quad.params,
                           g.f9.briere.params=g.f9.briere.params,
                           lf.f12.quad.params=lf.f12.quad.params,
                           lf.f12.briere.params=lf.f12.briere.params,
                           br.f12.briere.params=br.f12.briere.params,
                           br.f12.quad.params=br.f12.quad.params,
                           le.f12.briere.params=le.f12.briere.params,
                           le.f12.quad.params=le.f12.quad.params,
                           g.f12.quad.params=g.f12.quad.params,
                           g.f12.briere.params=g.f12.briere.params) 

sink("tables/model.params_list.txt") #file name
print(models.params_list) #object to be exported as text
sink() #exporting to source file location

