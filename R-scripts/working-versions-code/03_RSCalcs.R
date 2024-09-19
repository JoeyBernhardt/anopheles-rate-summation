## Rate Summation Project
## Script for performing rate summation calculations
## Written by Marta Shocket in 2022

##	Table of contents:
##
##	1. Set-up workspace
##  2. Generate hourly temperatures with the Parton-Logan model
##  3. Perform rate summation calculations on traits
##  4. Process rate summation calculation output for plotting
##  5. Join RS output with corresponding empirical TPCs and plot each panel
##  6. Manuscript Figure 2



##########
###### 1. Set-up workspace
##########

##### Load Packages
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(patchwork)
library(progress)
#library(gridExtra)
library(RColorBrewer)
library(progress)

##### Load functions
source("R-scripts/working-versions-code/00_RSProjectFunctions.R")

##### Load fitted TPC models
load("saved-posteriors/constant_lifespan.Rdata")
load("saved-posteriors/constant_gamma.Rdata")
load("saved-posteriors/constant_eggs.Rdata")
load("saved-posteriors/constant_biterate.Rdata")

load("saved-posteriors/dtr9_lifespan.Rdata")
load("saved-posteriors/dtr9_gamma.Rdata")
load("saved-posteriors/dtr9_eggs.Rdata")
load("saved-posteriors/dtr9_biterate.Rdata")

load("saved-posteriors/dtr12_lifespan.Rdata")
load("saved-posteriors/dtr12_gamma.Rdata")
load("saved-posteriors/dtr12_eggs.Rdata")
load("saved-posteriors/dtr12_biterate.Rdata")

load("saved-posteriors/pea.Rdata")
load("saved-posteriors/MDR.Rdata")
load("saved-posteriors/EIP50.Rdata")
load("saved-posteriors/bc.Rdata")
load("saved-posteriors/temps.Rdata")

# load("saved-posteriors/constant_lifespan_negpred.Rdata")

##### Extract matrices of predicted trait values
# NOTE: These are getting passed to dplyr so they need to be data frames here
# for the sensitivity analysis the predictions need to be matrices *NOT* data frames
predictions_bite_rate_constant <- as.data.frame(model_bite_rate_constant$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
predictions_bite_rate_dtr9 <- as.data.frame(model_bite_rate_dtr9$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
predictions_bite_rate_dtr12 <- as.data.frame(model_bite_rate_dtr12$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)

predictions_lifespan_constant <- as.data.frame(model_lifespan_constant$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
predictions_lifespan_dtr9 <- as.data.frame(model_lifespan_dtr9$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
predictions_lifespan_dtr12 <- as.data.frame(model_lifespan_dtr12$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)

predictions_eggs_constant <- as.data.frame(model_eggs_constant$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
predictions_eggs_dtr9 <- as.data.frame(model_eggs_dtr9$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
predictions_eggs_dtr12 <- as.data.frame(model_eggs_dtr12$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)

predictions_gamma_constant <- as.data.frame(model_gamma_constant$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
predictions_gamma_dtr9 <- as.data.frame(model_gamma_dtr9$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
predictions_gamma_dtr12 <- as.data.frame(model_gamma_dtr12$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)

predictions_bc <- as.data.frame(model_bc$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
predictions_eip50 <- as.data.frame(model_eip50$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
predictions_pea <- as.data.frame(model_pea$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)
predictions_mdr <- as.data.frame(model_mdr$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)

predictions_lifespan_constant_negpred <- as.data.frame(model_lifespan_constant_negpred$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)

##########
###### 2. Generate hourly temperatures with the Parton-Logan model
##########

# Temperature gradient that matches TPC predictions
Temp.xs <- seq(0, 45, 0.1)

# Generate hourly temperature sequences across the temperature gradient
LPtemps_dtr9 <- LoganPartonCalc(dtr = 9, Temp.xs)
LPtemps_dtr12 <- LoganPartonCalc(dtr = 12, Temp.xs)

# Set negative temperature values to 0, since TPCs stop at 0 on the low end
LPtemps_dtr9[LPtemps_dtr9 < 0 ] <- 0
LPtemps_dtr12[LPtemps_dtr12 < 0 ] <- 0

# Set temperature values > 45 to 45, since TPCs stop at 45 on the hight end
LPtemps_dtr9[LPtemps_dtr9 > 45 ] <- 45
LPtemps_dtr12[LPtemps_dtr12 > 45 ] <- 45



##########
###### 3. Perform rate summation calculations on traits
##########

# These take about 1.5-2 minutes each to run

predictions_bite_rate_rs9 <- RSCalcTempGrad(predictions_bite_rate_constant, LPtemps_dtr9, Temp.xs)
predictions_bite_rate_rs12 <- RSCalcTempGrad(predictions_bite_rate_constant, LPtemps_dtr12, Temp.xs)

predictions_lifespan_rs9 <- RSCalcTempGrad(predictions_lifespan_constant, LPtemps_dtr9, Temp.xs)
predictions_lifespan_rs12 <- RSCalcTempGrad(predictions_lifespan_constant, LPtemps_dtr12, Temp.xs)

predictions_eggs_rs9 <- RSCalcTempGrad(predictions_eggs_constant, LPtemps_dtr9, Temp.xs)
predictions_eggs_rs12 <- RSCalcTempGrad(predictions_eggs_constant, LPtemps_dtr12, Temp.xs)

predictions_gamma_rs9 <- RSCalcTempGrad(predictions_gamma_constant, LPtemps_dtr9, Temp.xs)
predictions_gamma_rs12 <- RSCalcTempGrad(predictions_gamma_constant, LPtemps_dtr12, Temp.xs)

predictions_bc_rs9 <- RSCalcTempGrad(predictions_bc, LPtemps_dtr9, Temp.xs)
predictions_bc_rs12 <- RSCalcTempGrad(predictions_bc, LPtemps_dtr12, Temp.xs)

predictions_eip50_rs9 <- RSCalcTempGrad(predictions_eip50, LPtemps_dtr9, Temp.xs)
predictions_eip50_rs12 <- RSCalcTempGrad(predictions_eip50, LPtemps_dtr12, Temp.xs)

predictions_pea_rs9 <- RSCalcTempGrad(predictions_pea, LPtemps_dtr9, Temp.xs)
predictions_pea_rs12 <- RSCalcTempGrad(predictions_pea, LPtemps_dtr12, Temp.xs)

predictions_mdr_rs9 <- RSCalcTempGrad(predictions_mdr, LPtemps_dtr9, Temp.xs)
predictions_mdr_rs12 <- RSCalcTempGrad(predictions_mdr, LPtemps_dtr12, Temp.xs)

# predictions_lifespan_rs9_negpred <- RSCalcTempGrad(predictions_lifespan_constant, LPtemps_dtr9, Temp.xs)
# predictions_lifespan_rs12_negpred <- RSCalcTempGrad(predictions_lifespan_constant, LPtemps_dtr12, Temp.xs)

# Save output for sensitivity analysis
save(predictions_bite_rate_rs9, file = "data-processed/predictions_bite_rate_rs9.Rdata")
save(predictions_bite_rate_rs12, file = "data-processed/predictions_bite_rate_rs12.Rdata")

save(predictions_lifespan_rs9, file = "data-processed/predictions_lifespan_rs9.Rdata")
save(predictions_lifespan_rs12, file = "data-processed/predictions_lifespan_rs12.Rdata")

save(predictions_eggs_rs9, file = "data-processed/predictions_eggs_rs9.Rdata")
save(predictions_eggs_rs12, file = "data-processed/predictions_eggs_rs12.Rdata")

save(predictions_gamma_rs9, file = "data-processed/predictions_gamma_rs9.Rdata")
save(predictions_gamma_rs12, file = "data-processed/predictions_gamma_rs12.Rdata")

save(predictions_bc_rs9, file = "data-processed/predictions_bc_rs9.Rdata")
save(predictions_bc_rs12, file = "data-processed/predictions_bc_rs12.Rdata")

save(predictions_eip50_rs9, file = "data-processed/predictions_eip50_rs9.Rdata")
save(predictions_eip50_rs12, file = "data-processed/predictions_eip50_rs12.Rdata")

save(predictions_pea_rs9, file = "data-processed/predictions_pea_rs9.Rdata")
save(predictions_pea_rs12, file = "data-processed/predictions_pea_rs12.Rdata")

save(predictions_mdr_rs9, file = "data-processed/predictions_mdr_rs9.Rdata")
save(predictions_mdr_rs12, file = "data-processed/predictions_mdr_rs12.Rdata")



##########
###### 4. Process rate summation calculation output for plotting
##########

##### Load saved data frames of predicted trait values from RS calculations
load("data-processed/predictions_bite_rate_rs9.Rdata")
load("data-processed/predictions_bite_rate_rs12.Rdata")

load("data-processed/predictions_lifespan_rs9.Rdata")
load("data-processed/predictions_lifespan_rs12.Rdata")

load("data-processed/predictions_eggs_rs9.Rdata")
load("data-processed/predictions_eggs_rs12.Rdata")

load("data-processed/predictions_gamma_rs9.Rdata")
load("data-processed/predictions_gamma_rs12.Rdata")

load("data-processed/predictions_bc_rs9.Rdata")
load("data-processed/predictions_bc_rs12.Rdata")
 
load("data-processed/predictions_eip50_rs9.Rdata")
load("data-processed/predictions_eip50_rs12.Rdata")

load("data-processed/predictions_pea_rs9.Rdata")
load("data-processed/predictions_pea_rs12.Rdata")

load("data-processed/predictions_mdr_rs9.Rdata")
load("data-processed/predictions_mdr_rs12.Rdata")

#### Calculate TPC posterior summary data (means & quantiles)

predictions_bite_rate_rs9_summary <- calcPostQuants(predictions_bite_rate_rs9, "bite_rate_rs9", Temp.xs)
predictions_bite_rate_rs12_summary <- calcPostQuants(predictions_bite_rate_rs12, "bite_rate_rs12", Temp.xs)

predictions_lifespan_rs9_summary <- calcPostQuants(predictions_lifespan_rs9, "lifespan_rs9", Temp.xs)
predictions_lifespan_rs12_summary <- calcPostQuants(predictions_lifespan_rs12, "lifespan_rs12", Temp.xs)

predictions_eggs_rs9_summary <- calcPostQuants(predictions_eggs_rs9, "eggs_rs9", Temp.xs)
predictions_eggs_rs12_summary <- calcPostQuants(predictions_eggs_rs12, "eggs_rs12", Temp.xs)

predictions_gamma_rs9_summary <- calcPostQuants(predictions_gamma_rs9, "gamma_rs9", Temp.xs)
predictions_gamma_rs12_summary <- calcPostQuants(predictions_gamma_rs12, "gamma_rs12", Temp.xs)

predictions_bc_rs9_summary <- calcPostQuants(predictions_bc_rs9, "bc_rs9", Temp.xs)
predictions_bc_rs12_summary <- calcPostQuants(predictions_bc_rs12, "bc_rs12", Temp.xs)

predictions_eip50_rs9_summary <- calcPostQuants(predictions_eip50_rs9, "eip50_rs9", Temp.xs)
predictions_eip50_rs12_summary <- calcPostQuants(predictions_eip50_rs12, "eip50_rs12", Temp.xs)

predictions_pea_rs9_summary <- calcPostQuants(predictions_pea_rs9, "pea_rs9", Temp.xs)
predictions_pea_rs12_summary <- calcPostQuants(predictions_pea_rs12, "pea_rs12", Temp.xs)

predictions_mdr_rs9_summary <- calcPostQuants(predictions_mdr_rs9, "mdr_rs9", Temp.xs)
predictions_mdr_rs12_summary <- calcPostQuants(predictions_mdr_rs12, "mdr_rs12", Temp.xs)


### Get summary statistics of Tmin, Tmax, Topt, and Tbreadth

params_bite_rate_rs9_summary <- extractDerivedTPC(predictions_bite_rate_rs9, "bite_rate_rs9", Temp.xs)
params_bite_rate_rs12_summary <- extractDerivedTPC(predictions_bite_rate_rs12, "bite_rate_rs12", Temp.xs)

params_lifespan_rs9_summary <- extractDerivedTPC(predictions_lifespan_rs9, "lifespan_rs9", Temp.xs)
params_lifespan_rs12_summary <- extractDerivedTPC(predictions_lifespan_rs12, "lifespan_rs12", Temp.xs)

params_eggs_rs9_summary <- extractDerivedTPC(predictions_eggs_rs9, "eggs_rs9", Temp.xs)
params_eggs_rs12_summary <- extractDerivedTPC(predictions_eggs_rs12, "eggs_rs12", Temp.xs)

# NOTE: need to load fitted TPC summaries below for these calculations to work

# Calculate % reduction in maximum value at Topt for each fluctuation model vs constant temperature model
1 - max(predictions_bite_rate_dtr9_summary$median)/max(predictions_bite_rate_constant_summary$median) # 25.1%
1 - max(predictions_bite_rate_dtr12_summary$median)/max(predictions_bite_rate_constant_summary$median) # 23.5%
1 - max(predictions_bite_rate_rs9_summary$median)/max(predictions_bite_rate_constant_summary$median) # 5.1%
1 - max(predictions_bite_rate_rs12_summary$median)/max(predictions_bite_rate_constant_summary$median) # 8.9%

1 - max(predictions_lifespan_dtr9_summary$median)/max(predictions_lifespan_constant_summary$median) # 2.7%
1 - max(predictions_lifespan_dtr12_summary$median)/max(predictions_lifespan_constant_summary$median) # -10.7% - actually higher!
1 - max(predictions_lifespan_rs9_summary$median)/max(predictions_lifespan_constant_summary$median) # 3.0%
1 - max(predictions_lifespan_rs12_summary$median)/max(predictions_lifespan_constant_summary$median) # 5.2%

1 - max(predictions_eggs_dtr9_summary$median)/max(predictions_eggs_constant_summary$median) # 7.9%
1 - max(predictions_eggs_dtr12_summary$median)/max(predictions_eggs_constant_summary$median) # 14.8%
1 - max(predictions_eggs_rs9_summary$median)/max(predictions_eggs_constant_summary$median) # 6.5%
1 - max(predictions_eggs_rs12_summary$median)/max(predictions_eggs_constant_summary$median) # 11.5%

# Calculate % reduction in empirical fluctuating model vs rate summation fluctuating model
1 - max(predictions_bite_rate_dtr9_summary$median)/max(predictions_bite_rate_rs9_summary$median) # 21.1%
1 - max(predictions_bite_rate_dtr12_summary$median)/max(predictions_bite_rate_rs12_summary$median) # 16.1%

1 - max(predictions_lifespan_dtr9_summary$median)/max(predictions_lifespan_rs9_summary$median) # -0.2%
1 - max(predictions_lifespan_dtr12_summary$median)/max(predictions_lifespan_rs12_summary$median) # -16.8% - actually higher!

1 - max(predictions_eggs_dtr9_summary$median)/max(predictions_eggs_rs9_summary$median) # 1.5%
1 - max(predictions_eggs_dtr12_summary$median)/max(predictions_eggs_rs12_summary$median) # 3.7%

##########
###### 5. Join RS output with corresponding empirical TPCs
##########

### Biting Rate (a) 

### load empirical TPC summaries
predictions_bite_rate_constant_summary <- read_csv("data-processed/predictions_bite_rate_constant_summary.csv")
predictions_bite_rate_dtr9_summary <- read_csv("data-processed/predictions_bite_rate_dtr9_summary.csv")
predictions_bite_rate_dtr12_summary <- read_csv("data-processed/predictions_bite_rate_dtr12_summary.csv")

params_bite_rate_constant_summary <- read.csv("data-processed/params_bite_rate_constant_summary.csv")
params_bite_rate_dtr9_summary <- read.csv("data-processed/params_bite_rate_dtr9_summary.csv")
params_bite_rate_dtr12_summary <- read.csv("data-processed/params_bite_rate_dtr12_summary.csv")

### combine empirically measured and rate summation treatments (separately for dtr 9 and for dtr12) - plus constant for comparison
dtr9_rs9_bite_rate_predictions <- bind_rows(predictions_bite_rate_constant_summary, predictions_bite_rate_dtr9_summary, predictions_bite_rate_rs9_summary)
dtr12_rs12_bite_rate_predictions <- bind_rows(predictions_bite_rate_constant_summary, predictions_bite_rate_dtr12_summary, predictions_bite_rate_rs12_summary)

dtr9_rs9_bite_rate_params <- bind_rows(params_bite_rate_constant_summary, params_bite_rate_dtr9_summary, params_bite_rate_rs9_summary) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))
dtr12_rs12_bite_rate_params <- bind_rows(params_bite_rate_constant_summary, params_bite_rate_dtr12_summary, params_bite_rate_rs12_summary) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

### Lifespan (lf) 

### load empirical TPC summaries
predictions_lifespan_constant_summary <- read_csv("data-processed/predictions_lifespan_constant_summary.csv")
predictions_lifespan_dtr9_summary <- read_csv("data-processed/predictions_lifespan_dtr9_summary.csv")
predictions_lifespan_dtr12_summary <- read_csv("data-processed/predictions_lifespan_dtr12_summary.csv")

params_lifespan_constant_summary <- read.csv("data-processed/params_lifespan_constant_summary.csv")
params_lifespan_dtr9_summary <- read.csv("data-processed/params_lifespan_dtr9_summary.csv")
params_lifespan_dtr12_summary <- read.csv("data-processed/params_lifespan_dtr12_summary.csv")

### combine empirically measured and rate summation treatments (separately for dtr 9 and for dtr12) - plus constant for comparison
dtr9_rs9_lifespan_predictions <- bind_rows(predictions_lifespan_constant_summary, predictions_lifespan_dtr9_summary, predictions_lifespan_rs9_summary)
dtr12_rs12_lifespan_predictions <- bind_rows(predictions_lifespan_constant_summary, predictions_lifespan_dtr12_summary, predictions_lifespan_rs12_summary)

dtr9_rs9_lifespan_params <- bind_rows(params_lifespan_constant_summary, params_lifespan_dtr9_summary, params_lifespan_rs9_summary) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))
dtr12_rs12_lifespan_params <- bind_rows(params_lifespan_constant_summary, params_lifespan_dtr12_summary, params_lifespan_rs12_summary) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

### Lifetime eggs (B) 

### load empirical TPC summaries
predictions_eggs_constant_summary <- read_csv("data-processed/predictions_eggs_constant_summary.csv")
predictions_eggs_dtr9_summary <- read_csv("data-processed/predictions_eggs_dtr9_summary.csv")
predictions_eggs_dtr12_summary <- read_csv("data-processed/predictions_eggs_dtr12_summary.csv")

params_eggs_constant_summary <- read.csv("data-processed/params_eggs_constant_summary.csv")
params_eggs_dtr9_summary <- read.csv("data-processed/params_eggs_dtr9_summary.csv")
params_eggs_dtr12_summary <- read.csv("data-processed/params_eggs_dtr12_summary.csv")

### combine empirically measured and rate summation treatments (separately for dtr 9 and for dtr12) - plus constant for comparison
dtr9_rs9_eggs_predictions <- bind_rows(predictions_eggs_constant_summary, predictions_eggs_dtr9_summary, predictions_eggs_rs9_summary)
dtr12_rs12_eggs_predictions <- bind_rows(predictions_eggs_constant_summary, predictions_eggs_dtr12_summary, predictions_eggs_rs12_summary)

dtr9_rs9_eggs_params <- bind_rows(params_eggs_constant_summary, params_eggs_dtr9_summary, params_eggs_rs9_summary) %>%
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))
dtr12_rs12_eggs_params <- bind_rows(params_eggs_constant_summary, params_eggs_dtr12_summary, params_eggs_rs12_summary) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

##########
###### 6. Manuscript Figure 2
##########

### Biting Rate (a) panels

bite_rate_rs9_plot <- dtr9_rs9_bite_rate_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), size = 0.6) +
	geom_text(x = -13.5, y = 0.25, size = 7, angle = 90, label = "DTR 9") +
	geom_text(x = 25, y = 0.6, size = 7, label = "Bite rate (a)") +
	scale_color_manual(values = c(c_constant, c_emp09, c_rstraits09), labels = c("Constant", "Empirically fit", "Rate summation")) +
	scale_fill_manual(values = c(ct_constant, ct_emp09, ct_rstraits09), labels = c("Constant", "Empirically fit", "Rate summation")) +
	scale_linetype_manual(values = c(1, 1, 2), labels = c("Constant", "Empirically fit", "Rate summation")) +
	ylab(parse(text = "Bite~rate~(day^-1)")) + xlab("Temperature (°C)") +
	theme_classic() +
	theme(axis.title.x = element_blank(), legend.position = c(0.215, 0.80), legend.title=element_blank(), legend.text=element_text(size=11)) +
	theme(plot.margin = unit(c(2.5, 1, 1, 2.5), "lines"), axis.title.y = element_text(vjust = -2)) +
	coord_cartesian(clip = "off") +
	annotate("text", x = 0, y = 0.53, label = "A", size = 5)

params_plot_bite_rate_rs9 <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = dtr9_rs9_bite_rate_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_manual(values = c(c_constant, c_emp09, c_rstraits09), labels = c("Constant", "DTR 9", "RS 9")) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 0, label = "D", size = 5)

bite_rate_rs12_plot <- dtr12_rs12_bite_rate_predictions %>%
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), size = 0.6) +
	geom_text(x = -13.5, y = 0.25, size = 7, angle = 90, label = "DTR 12") +
	scale_color_manual(values = c(c_constant, c_emp12, c_rstraits12), labels = c("Constant", "Empirically fit", "Rate summation")) +
	scale_fill_manual(values = c(ct_constant, ct_emp12, ct_rstraits12), labels = c("Constant", "Empirically fit", "Rate summation")) +
	scale_linetype_manual(values = c(1, 1, 2), labels = c("Constant", "Empirically fit", "Rate summation")) +
	ylab(parse(text = "Bite~rate~(day^-1)")) + xlab("Temperature (°C)") +
	theme_classic() +
	theme(axis.title.x = element_blank(), legend.position = c(0.215, 0.80), legend.title=element_blank(), legend.text=element_text(size=11)) +
	theme(plot.margin = unit(c(1, 1, 1, 2.5), "lines"), axis.title.y = element_text(vjust = -2)) +
	coord_cartesian(clip = "off") +
	annotate("text", x = 0, y = 0.53, label = "G", size = 5)

params_plot_bite_rate_rs12 <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = dtr12_rs12_bite_rate_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_manual(values = c(c_constant, c_emp12, c_rstraits12), labels = c("Constant", "DTR 12", "RS 12")) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 0, label = "J", size = 5)

bite_rate9_plot <- bite_rate_rs9_plot / params_plot_bite_rate_rs9 + plot_layout(heights = c(3, 0.5))
bite_rate12_plot <- bite_rate_rs12_plot / params_plot_bite_rate_rs12 + plot_layout(heights = c(3, 0.5))

### Lifespan (lf) panels

lifespan_rs9_plot <- dtr9_rs9_lifespan_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), size = 0.6) +
	geom_text(x = 25, y = 44, size = 7, label = "Lifespan (lf)") +
	scale_color_manual(values = c(c_constant, c_emp09, c_rstraits09), labels = c("Constant", "Empirically fit", "Rate summation")) +
	scale_fill_manual(values = c(ct_constant, ct_emp09, ct_rstraits09), labels = c("Constant", "Empirically fit", "Rate summation")) +
	scale_linetype_manual(values = c(1, 1, 2), labels = c("Constant", "Empirically fitted", "Rate summation")) +
	ylab("Lifespan (days)") + xlab("Temperature (°C)") +
	theme_classic() +
	theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_text(vjust = -1)) +
	coord_cartesian(clip = "off") +
	theme(plot.margin = unit(c(2.5, 1, 1, 1), "lines")) +
	annotate("text", x = 0, y = 39, label = "B", size = 5)

params_plot_lifespan_rs9 <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = dtr9_rs9_lifespan_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_manual(values = c(c_constant, c_emp09, c_rstraits09), labels = c("Constant", "DTR 9", "RS 9")) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 0, label = "E", size = 5)

lifespan_rs12_plot <- dtr12_rs12_lifespan_predictions %>%
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), size = 0.6) +
	scale_color_manual(values = c(c_constant, c_emp12, c_rstraits12), labels = c("Constant", "Empirically fit", "Rate summation")) +
	scale_fill_manual(values = c(ct_constant, ct_emp12, ct_rstraits12), labels = c("Constant", "Empirically fit", "Rate summation")) +
	scale_linetype_manual(values = c(1, 1, 2), labels = c("Constant", "Empirically fitted", "Rate summation")) +
	ylab("Lifespan (days)") + xlab("Temperature (°C)") +
	theme_classic() +
	theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_text(vjust = -1)) + 
	annotate("text", x = 0, y = 44, label = "H", size = 5)

params_plot_lifespan_rs12 <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = dtr12_rs12_lifespan_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_manual(values = c(c_constant, c_emp12, c_rstraits12), labels = c("Constant", "DTR 12", "RS 12")) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 0, label = "K", size = 5)

lifespan9_plot <- lifespan_rs9_plot / params_plot_lifespan_rs9 + plot_layout(heights = c(3, 0.5))
lifespan12_plot <- lifespan_rs12_plot / params_plot_lifespan_rs12 + plot_layout(heights = c(3, 0.5))

### Lifetime eggs (B) panels

eggs_rs9_plot <- dtr9_rs9_eggs_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), size = 0.6) +
	geom_text(x = 25, y = 495, size = 7, label = "Lifetime eggs (B)") +
	scale_color_manual(values = c(c_constant, c_emp09, c_rstraits09), labels = c("Constant", "Empirically fit", "Rate summation")) +
	scale_fill_manual(values = c(ct_constant, ct_emp09, ct_rstraits09), labels = c("Constant", "Empirically fit", "Rate summation")) +
	scale_linetype_manual(values = c(1, 1, 2), labels = c("Constant", "Empirically fit", "Rate summation")) +
	ylab("Lifetime eggs") + xlab("Temperature (°C)") +
	theme_classic() +
	theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_text(vjust = -2)) +
	coord_cartesian(clip = "off") +
	theme(plot.margin = unit(c(2.5, 1, 1, 1), "lines")) +
	annotate("text", x = 0, y = 440, label = "C", size = 5)

params_plot_eggs_rs9 <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = dtr9_rs9_eggs_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_manual(values = c(c_constant, c_emp09, c_rstraits09), labels = c("Constant", "DTR 9", "RS 9")) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 0, label = "F", size = 5)

eggs_rs12_plot <- dtr12_rs12_eggs_predictions %>%
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), size = 0.6) +
	scale_color_manual(values = c(c_constant, c_emp12, c_rstraits12), labels = c("Constant", "Empirically fit", "Rate summation")) +
	scale_fill_manual(values = c(ct_constant, ct_emp12, ct_rstraits12), labels = c("Constant", "Empirically fit", "Rate summation")) +
	scale_linetype_manual(values = c(1, 1, 2), labels = c("Constant", "Empirically fit", "Rate summation")) +
	ylab("Lifetime eggs") + xlab("Temperature (°C)") +
	theme_classic() +
	theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_text(vjust = -2)) + 
	annotate("text", x = 0, y = 440, label = "I", size = 5)

params_plot_eggs_rs12 <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = dtr12_rs12_eggs_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_manual(values = c(c_constant, c_emp12, c_rstraits12), labels = c("Constant", "DTR 12", "RS 12")) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 0, label = "L", size = 5)

eggs9_plot <- eggs_rs9_plot / params_plot_eggs_rs9 + plot_layout(heights = c(3, 0.5))
eggs12_plot <- eggs_rs12_plot / params_plot_eggs_rs12 + plot_layout(heights = c(3, 0.5))

### All panels together

Fig2_plots <- wrap_plots(bite_rate9_plot, lifespan9_plot, eggs9_plot,
						 bite_rate12_plot, lifespan12_plot, eggs12_plot, ncol = 3)
ggsave('figures/Fig2_plots-Sept2024.pdf', Fig2_plots, width = 15, height = 12)

# Old figure without the parameter panels
Fig2_plots_no_params <- wrap_plots(bite_rate_rs9_plot, lifespan_rs9_plot, eggs_rs9_plot,
						 bite_rate_rs12_plot, lifespan_rs12_plot, eggs_rs12_plot, ncol = 3)
ggsave('figures/Fig2_plots-Jan4.pdf', Fig2_plots_no_params, width = 15, height = 10)

