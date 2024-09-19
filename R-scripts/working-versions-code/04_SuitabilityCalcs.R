## Rate Summation Project
## Script for calculating suitability (S(T) = relative R0)
## Written by Marta Shocket in 2022

##	Table of contents:
##
##	1. Set-up workspace
##  2. Define S(T)/R0 function
##	3. Calculate versions of Suitability S(T)
##  4. Process output
##  5. Manuscript Figure 3



##########
###### 1. Set-up workspace
##########

##### Load Packages
library(tidyverse)
library(progress)
library(gridExtra)
library(patchwork)

##### Load functions
source("R-scripts/working-versions-code/00_RSProjectFunctions.R")

# NOTE: Saved trait TPCs can be loaded in the top of file 03_RSCalcs.R


##########
###### 2. Define S(T)/R0 function
##########

# constant to keep lifespan from being numerically zero
# assume minimum survival time is half an hour
ec.lf = 1/48

# R0 formulation where: 1) mosquito density (M) depends on lifetime fecundity (B); 2) gamma (y) is substituted for exp^(-1/(lf*PDR))expression
R0eq = function(a, lf, B, y, bc, pEA, MDR) {
	M = B * pEA * MDR * (lf+ec.lf)
	R0 = (a^2 * bc * y * M * (lf+ec.lf) )^0.5
	return(R0)
}


##########
###### 5. Calculate versions of Suitability S(T)
##########

##### Summary of different versions of suitability calculations
#
# Five different versions of S(T) for specific pairwise comparisons:
# 	1) suit_const:			constant temperature TPCs
# 2) suit_empfluc:		empirically fit TPCs for a, lf, and B in fluctuating environments (other traits varying based on constant temperature TPCs)
# 3) suit_rs_traitsall:	rate summation on all trait TPCs
# 4) suit_rs_3traits:		rate summation on a, lf, and B (other traits varying based on constant temperature TPCs)
# 5) suit_rs_suit: 		rate summation on S(T) TPC
# 
# Pairwise comparisons and what they mean biologically or statistically:
# 1 vs 2 - Does fluctuating temperature (empirically) affect S(T)?
# 1 vs 3 - Does incorporating rate summation (on traits) affect predicted S(T)?
# 2 vs 4 - Does rate summation accurately predict the effect of fluctuating temperature on S(T)?
# 3 vs 5 - Does performing rate summation on traits vs. on S(T) TPC affect predicted S(T)?

# 1. All constant temperature TPCs
suit_const <- R0eq(predictions_bite_rate_constant, predictions_lifespan_constant, predictions_eggs_constant, predictions_gamma_constant, 
				   predictions_bc, predictions_pea, predictions_mdr)

# 2A. Empirically fit TPCs for a, lf, and B in dtr9, other traits constant temperature TPCs
suit_empfluc_dtr09 <- R0eq(predictions_bite_rate_dtr9, predictions_lifespan_dtr9, predictions_eggs_dtr9, predictions_gamma_dtr9,
						  predictions_bc, predictions_pea, predictions_mdr)

# 2B. Empirically fit TPCs for a, lf, and B in dtr12, other traits constant temperature TPCs
suit_empfluc_dtr12 <- R0eq(predictions_bite_rate_dtr12, predictions_lifespan_dtr12, predictions_eggs_dtr12, predictions_gamma_dtr12,
						   predictions_bc, predictions_pea, predictions_mdr)

# 3A. Rate summation on all traits assuming dtr9
suit_rs_traitsall_09 <- R0eq(predictions_bite_rate_rs9, predictions_lifespan_rs9, predictions_eggs_rs9, predictions_gamma_rs9,
							predictions_bc_rs9, predictions_pea_rs9, predictions_mdr_rs9)

# 3B. Rate summation on all traits assuming dtr12
suit_rs_traitsall_12 <- R0eq(predictions_bite_rate_rs12, predictions_lifespan_rs12, predictions_eggs_rs12, predictions_gamma_rs12,
							 predictions_bc_rs12, predictions_pea_rs12, predictions_mdr_rs12)

# 4A. Rate summation on a, lf, and B assuming dtr9, other traits constant temperature TPCs
suit_rs_3traits_09 <- R0eq(predictions_bite_rate_rs9, predictions_lifespan_rs9, predictions_eggs_rs9, predictions_gamma_rs9,
						  predictions_bc, predictions_pea, predictions_mdr)

# 4B. Rate summation on a, lf, and B assuming dtr12, other traits constant temperature TPCs
suit_rs_3traits_12 <- R0eq(predictions_bite_rate_rs12, predictions_lifespan_rs12, predictions_eggs_rs12, predictions_gamma_rs12,
						   predictions_bc, predictions_pea, predictions_mdr)

# 5A. Rate summation on suit_const assuming dtr9
suit_rs_suit_09 <- RSCalcTempGrad(suit_const, LPtemps_dtr9, Temp.xs)

# 5A. Rate summation on suit_const assuming dtr12
suit_rs_suit_12 <- RSCalcTempGrad(suit_const, LPtemps_dtr12, Temp.xs)


# ##### Change column names to temperature gradient values
# colnames(suit_const) <- colnames(suit_empfluc_dtr09) <- colnames(suit_empfluc_dtr12) <- Temp.xs
# colnames(suit_rs_traitsall_09) <- colnames(suit_rs_traitsall_12) <- colnames(suit_rs_3traits_09) <- colnames(suit_rs_3traits_12) <- Temp.xs
# colnames(suit_rs_suit_09) <- colnames(suit_rs_suit_12) <- Temp.xs


##########
###### 6. Process output
##########

# # Use a threshold for very small values of suitability (0.1) - magnitude of maximum S(T) is ~50-75
# suit_const[suit_const < 0.1] <- 0 
# suit_empfluc_dtr09[suit_empfluc_dtr09 < 0.1] <- 0 
# suit_empfluc_dtr12[suit_empfluc_dtr12 < 0.1] <- 0 
# suit_rs_traitsall_09[suit_rs_traitsall_09 < 0.1] <- 0 
# suit_rs_traitsall_12[suit_rs_traitsall_12 < 0.1] <- 0 
# suit_rs_3traits_09[suit_rs_3traits_09 < 0.1] <- 0 
# suit_rs_3traits_12[suit_rs_3traits_12 < 0.1] <- 0 
# suit_rs_suit_09[suit_rs_suit_09 < 0.1] <- 0 
# suit_rs_suit_12[suit_rs_suit_12 < 0.1] <- 0 


# Calculate quantiles for plotting
suit_const_summary <- calcPostQuants(suit_const, "suit_const", Temp.xs)

suit_empfluc_dtr09_summary <- calcPostQuants(suit_empfluc_dtr09, "suit_empfluc_dtr09", Temp.xs)
suit_empfluc_dtr12_summary <- calcPostQuants(suit_empfluc_dtr12, "suit_empfluc_dtr12", Temp.xs)

suit_rs_traitsall_09_summary <- calcPostQuants(suit_rs_traitsall_09, "suit_rs_traitsall_09", Temp.xs)
suit_rs_traitsall_12_summary <- calcPostQuants(suit_rs_traitsall_12, "suit_rs_traitsall_12", Temp.xs)

suit_rs_3traits_09_summary <- calcPostQuants(suit_rs_3traits_09, "suit_rs_3traits_09", Temp.xs)
suit_rs_3traits_12_summary <- calcPostQuants(suit_rs_3traits_12, "suit_rs_3traits_12", Temp.xs)

suit_rs_suit_09_summary <- calcPostQuants(suit_rs_suit_09, "suit_rs_suit_09", Temp.xs)
suit_rs_suit_12_summary <- calcPostQuants(suit_rs_suit_12, "suit_rs_suit_12", Temp.xs)

# Calculate maximum predicted R0 value for each model
max(suit_const_summary$median)
max(suit_empfluc_dtr09_summary$median)
max(suit_empfluc_dtr12_summary$median)
max(suit_rs_traitsall_09_summary$median)
max(suit_rs_traitsall_12_summary$median)
max(suit_rs_3traits_09_summary$median)
max(suit_rs_3traits_12_summary$median)
max(suit_rs_suit_09_summary$median)
max(suit_rs_suit_12_summary$median)

# Calculate % reduction in maximum R0 at Topt for each fluctuation model vs constant temperature model
1 - max(suit_empfluc_dtr09_summary$median)/max(suit_const_summary$median) # 32.0%
1 - max(suit_empfluc_dtr12_summary$median)/max(suit_const_summary$median) # 33.8%
1 - max(suit_rs_traitsall_09_summary$median)/max(suit_const_summary$median) # 18.1%
1 - max(suit_rs_traitsall_12_summary$median)/max(suit_const_summary$median) # 30.6%
1 - max(suit_rs_3traits_09_summary$median)/max(suit_const_summary$median) # 10.0%
1 - max(suit_rs_3traits_12_summary$median)/max(suit_const_summary$median) # 17.1%
1 - max(suit_rs_suit_09_summary$median)/max(suit_const_summary$median) # 19.9 %
1 - max(suit_rs_suit_12_summary$median)/max(suit_const_summary$median) # 32.0 %

# Scale quantiles to get relative R0/S(T)
#	  NOTE: scale everything to the median of the constant temperature model

suit_const_summary_rel <- calcRelSuit(suit_const_summary, suit_const_summary)

suit_empfluc_dtr09_summary_rel <- calcRelSuit(suit_empfluc_dtr09_summary, suit_const_summary)
suit_empfluc_dtr12_summary_rel <- calcRelSuit(suit_empfluc_dtr12_summary, suit_const_summary)

suit_rs_traitsall_09_summary_rel <- calcRelSuit(suit_rs_traitsall_09_summary, suit_const_summary)
suit_rs_traitsall_12_summary_rel <- calcRelSuit(suit_rs_traitsall_12_summary, suit_const_summary)

suit_rs_3traits_09_summary_rel <- calcRelSuit(suit_rs_3traits_09_summary, suit_const_summary)
suit_rs_3traits_12_summary_rel <- calcRelSuit(suit_rs_3traits_12_summary, suit_const_summary)

suit_rs_suit_09_summary_rel <- calcRelSuit(suit_rs_suit_09_summary, suit_const_summary)
suit_rs_suit_12_summary_rel <- calcRelSuit(suit_rs_suit_12_summary, suit_const_summary)


# Calculate the maximum value of the median S(T) from the constant temperature model
suitability_scale <- max(suit_const_summary$median)


# Get summary statistics of Tmin, Tmax, Topt, and Tbreadth 
params_suit_const <- extractDerivedTPC(suit_const, "suit_const", Temp.xs)

params_suit_empfluc_dtr09 <- extractDerivedTPC(suit_empfluc_dtr09, "suit_empfluc_dtr09", Temp.xs)
params_suit_empfluc_dtr12 <- extractDerivedTPC(suit_empfluc_dtr12, "suit_empfluc_dtr12", Temp.xs)

params_suit_rs_traitsall_09 <- extractDerivedTPC(suit_rs_traitsall_09, "suit_rs_traitsall_09", Temp.xs)
params_suit_rs_traitsall_12 <- extractDerivedTPC(suit_rs_traitsall_12, "suit_rs_traitsall_12", Temp.xs)

params_suit_rs_3traits_09 <- extractDerivedTPC(suit_rs_3traits_09, "suit_rs_3traits_09", Temp.xs)
params_suit_rs_3traits_12 <- extractDerivedTPC(suit_rs_3traits_12, "suit_rs_3traits_12", Temp.xs)

params_suit_rs_suit_09 <- extractDerivedTPC(suit_rs_suit_09, "suit_rs_suit_09", Temp.xs)
params_suit_rs_suit_12 <- extractDerivedTPC(suit_rs_suit_12, "suit_rs_suit_12", Temp.xs)



##########
###### 7. Figures
##########

##### Summary of different versions of suitability calculations
#
# Five different versions of S(T) for specific pairwise comparisons:
# 	1) suit_const:			constant temperature TPCs
#	2) suit_empfluc:		empirically fit TPCs for a, lf, and B in fluctuating environments (other traits varying based on constant temperature TPCs)
#	3) suit_rs_traitsall:	rate summation on all trait TPCs
#	4) suit_rs_3traits:		rate summation on a, lf, and B (other traits varying based on constant temperature TPCs)
#	5) suit_rs_suit: 		rate summation on S(T) TPC
#
# Pairwise comparisons and what they mean biologically or statistically:
#	1 vs 2 - Does fluctuating temperature (empirically) affect S(T)?
#	1 vs 3 - Does incorporating rate summation (on traits) affect predicted S(T)?
#	2 vs 4 - Does rate summation accurately predict the effect of fluctuating temperature on S(T)?
#	3 vs 5 - Does performing rate summation on traits vs. on S(T) TPC affect predicted S(T)?


################################ Combining treatments for different comparisons

#### Comparison #1: combine all empirical treatments
emp_suit_predictions <- bind_rows(suit_const_summary, suit_empfluc_dtr09_summary, suit_empfluc_dtr12_summary)
#emp_suit_predictions_rel <- bind_rows(suit_const_summary_rel, suit_empfluc_dtr09_summary_rel, suit_empfluc_dtr12_summary_rel)
emp_suit_params <- bind_rows(params_suit_const, params_suit_empfluc_dtr09, params_suit_empfluc_dtr12) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

#### Comparison #2: combine rate summation on 3 traits and empirical dtr treatments
rsaccuracy_suit_predictions <- bind_rows(suit_rs_3traits_09_summary, suit_rs_3traits_12_summary, suit_empfluc_dtr09_summary, suit_empfluc_dtr12_summary)
#rsaccuracy_suit_predictions_rel <- bind_rows(suit_rs_3traits_09_summary_rel, suit_rs_3traits_12_summary_rel, suit_empfluc_dtr09_summary_rel, suit_empfluc_dtr12_summary_rel)
rsaccuracy_suit_params <- bind_rows(params_suit_rs_3traits_09, params_suit_rs_3traits_12, params_suit_empfluc_dtr09, params_suit_empfluc_dtr12) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

#### Comparison #3: combine rate summation on all traits and rate summation on suitability
rslevel_suit_predictions <- bind_rows(suit_rs_traitsall_09_summary, suit_rs_traitsall_12_summary, suit_rs_suit_09_summary, suit_rs_suit_12_summary)
#rslevel_suit_predictions_rel <- bind_rows(suit_rs_traitsall_09_summary_rel, suit_rs_traitsall_12_summary_rel, suit_rs_suit_09_summary_rel, suit_rs_suit_12_summary_rel)
rslevel_suit_params <- bind_rows(params_suit_rs_traitsall_09, params_suit_rs_traitsall_12, params_suit_rs_suit_09, params_suit_rs_suit_12) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

# #### Comparison #4: combine rate summation on all traits and constant
rseffect_suit_predictions <- bind_rows(suit_const_summary, suit_rs_traitsall_09_summary, suit_rs_traitsall_12_summary)
#rseffect_suit_predictions_rel <- bind_rows(suit_const_summary_rel, suit_rs_traitsall_09_summary_rel, suit_rs_traitsall_12_summary_rel)
rseffect_suit_params <- bind_rows(params_suit_const, params_suit_rs_traitsall_09, params_suit_rs_traitsall_12) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

#### Combine all suitability calculations - for supplemental figure
all_suit_predictions <- bind_rows(suit_const_summary, suit_empfluc_dtr09_summary, suit_empfluc_dtr12_summary, suit_rs_3traits_09_summary, suit_rs_3traits_12_summary,
								  suit_rs_traitsall_09_summary, suit_rs_traitsall_12_summary, suit_rs_suit_09_summary, suit_rs_suit_12_summary)


#### Comparison #1: combine all empirical treatments

emp_suit_plot <- emp_suit_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), size = 0.6) +
	scale_color_manual(values = c(c_constant, c_emp09, c_emp12), labels = c("1: Constant", "2: Empirical DTR 9", "2: Empirical DTR 12")) +
	scale_fill_manual(values = c(ct_constant, ct_emp09, ct_emp12), labels = c("1: Constant", "2: Empirical DTR 9", "2: Empirical DTR 12")) +
	scale_linetype_manual(values = c(1, 1, 1), labels = c("1: Constant", "2: Empirical DTR 9", "2: Empirical DTR 12")) +
	ylab("Suitability for transmission - S(T)") + xlab("Temperature (°C)") + xlim(5, 43) +
	theme_classic() +
	theme(legend.position = c(0.2, 0.75), legend.title=element_blank(), axis.title.x = element_blank(),
		  axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
	annotate("text", x = 5, y = 67, label = "A", size = 5)

params_emp_suit_plot <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = emp_suit_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") + ylim(5, 43) +
	scale_color_manual(values = c(c_constant, c_emp09, c_emp12)) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 5, label = "E", size = 5)

emp_suit_plot_combined <- emp_suit_plot / params_emp_suit_plot + plot_layout(heights = c(3, 0.75))

#### Comparison #2: combine rate summation on 3 traits and empirical dtr treatments

rsaccuracy_suit_plot <- rsaccuracy_suit_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), size = 0.6) +
	scale_color_manual(values = c(c_emp09, c_emp12, c_rstraits09, c_rstraits12), labels = c("2: Empirical DTR 9", "2: Empirical DTR 12", "3: RS on 3 traits - DTR 9", "3: RS on 3 traits - DTR 12")) +
	scale_fill_manual(values = c(ct_emp09, ct_emp12, ct_rstraits09, ct_rstraits12), labels = c("2: Empirical DTR 9", "2: Empirical DTR 12", "3: RS on 3 traits - DTR 9", "3: RS on 3 traits - DTR 12")) +
	scale_linetype_manual(values = c(1, 1, 2, 2), labels = c("2: Empirical DTR 9", "2: Empirical DTR 12", "3: RS on 3 traits - DTR 9", "3: RS on 3 traits - DTR 12")) +
	ylab("Suitability for transmission - S(T)") + xlab("Temperature (°C)") + xlim(5, 43) +
	theme_classic() +
	theme(legend.position = c(0.2, 0.75), legend.title=element_blank(), axis.title.x = element_blank(),
		  axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
	annotate("text", x = 5, y = 60, label = "B", size = 5)

params_rsaccuracy_suit_plot <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = rsaccuracy_suit_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") + ylim(5, 43) +
	scale_color_manual(values = c(c_emp09, c_emp12, c_rstraits09, c_rstraits12)) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 5, label = "F", size = 5)

rsaccuracy_suit_plot_combined <- rsaccuracy_suit_plot / params_rsaccuracy_suit_plot + plot_layout(heights = c(3, 0.75))

#### Comparison #3: combine rate summation on all traits and rate summation on suitability

rslevel_suit_plot <- rslevel_suit_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), size = 0.6) +
	scale_color_manual(values = c(c_rssuit09, c_rssuit12, c_rstraits09, c_rstraits12), labels = c("5: RS on S(T) - DTR 9", "5: RS on S(T) - DTR 12", "4: RS on all traits - DTR 9", "4: RS on all traits - DTR 12")) +
	scale_fill_manual(values = c(ct_rssuit09, ct_rssuit12, ct_rstraits09, ct_rstraits12), labels = c("5: RS on S(T) - DTR 9", "5: RS on S(T) - DTR 12", "4: RS on all traits - DTR 9", "4: RS on all traits - DTR 12")) +
	scale_linetype_manual(values = c(3, 3, 2, 2), labels = c("5: RS on S(T) - DTR 9", "5: RS on S(T) - DTR 12", "4: RS on all traits - DTR 9", "4: RS on all traits - DTR 12")) +
	ylab("Suitability for transmission - S(T)") + xlab("Temperature (°C)") + xlim(5, 43) +
	theme_classic() +
	theme(legend.position = c(0.2, 0.75), legend.title=element_blank(), axis.title.x = element_blank(),
		  axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
	annotate("text", x = 5, y = 54, label = "C", size = 5)

params_rslevel_suit_plot <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = rslevel_suit_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") + ylim(5, 43) +
	scale_color_manual(values = c(c_rssuit09, c_rssuit12, c_rstraits09, c_rstraits12)) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 5, label = "G", size = 5)

rslevel_suit_plot_combined <- rslevel_suit_plot / params_rslevel_suit_plot + plot_layout(heights = c(3, 0.75))

#### Panel D: Lines for all models, without CIs

rseffect_suit_plot <- rseffect_suit_predictions %>%
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), size = 0.6) +
	scale_color_manual(values = c(c_constant, c_rstraits09, c_rstraits12), labels = c("1: Constant", "4: RS all traits - DTR 9", "4: RS all traits - DTR 12")) +
	scale_fill_manual(values = c(ct_constant, ct_rstraits09, ct_rstraits12), labels = c("1: Constant", "4: RS all traits - DTR 9", "4: RS all traits - DTR 12")) +
	scale_linetype_manual(values = c(1, 2, 2), labels = c("1: Constant", "4: RS all traits - DTR 9", "4: RS all traits - DTR 12")) +
	ylab("Suitability for transmission - S(T)") + xlab("Temperature (°C)") + xlim(5, 43) +
	theme_classic() +
	theme(legend.position = c(0.2, 0.75), legend.title=element_blank(), axis.title.x = element_blank(),
		  axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
	annotate("text", x = 5, y = 67, label = "D", size = 5)

params_rseffect_suit_plot <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = rseffect_suit_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") + ylim(5, 43) +
	scale_color_manual(values = c(c_constant, c_rstraits09, c_rstraits12)) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 5, label = "H", size = 5)

rseffect_suit_plot_combined <- rseffect_suit_plot / params_rseffect_suit_plot + plot_layout(heights = c(3, 0.75))

#### Combining plots

Fig3_plots <- wrap_plots(emp_suit_plot_combined, rsaccuracy_suit_plot_combined, rslevel_suit_plot_combined, rseffect_suit_plot_combined, nrow = 2)
ggsave('figures/Fig3_plots-Sept2024.pdf', Fig3_plots, width = 12, height = 12)

# Fig3_plots_no_params <- wrap_plots(emp_suit_plot, rsaccuracy_suit_plot, rslevel_suit_plot, nrow = 3)
# ggsave('figures/Fig3_plots-Jan4.pdf', Fig3_plots_no_params, width = 6, height = 12)



all_suit_plot <- all_suit_predictions %>% 
	ggplot() +
	geom_line(aes(x = temperature, y = mean, color = treatment, linetype = treatment), size = 0.8) +
	scale_color_manual(values = c(c_constant, c_emp09, c_emp12, c_emp09, c_emp12, c_emp09, c_emp12, c_emp09, c_emp12), labels = c("Constant", "DTR 9", "DTR 12", "RS 3 - DTR 9", "RS 3 - DTR 12", "RS all - DTR 9", "RS all - DTR 12", "RS S(T) - DTR 9", "RS S(T) - DTR 12")) +
	scale_linetype_manual(values = c(1, 1, 1, 2, 2, 2, 2, 3, 3), labels = c("Constant", "DTR 9", "DTR 12", "RS 3 - DTR 9", "RS 3 - DTR 12", "RS all - DTR 9", "RS all - DTR 12", "RS S(T) - DTR 9", "RS S(T) - DTR 12")) +
	ylab("Suitability for transmission - S(T)") + xlab("Temperature (°C)") + xlim(5, 43) +
	theme_classic() +
	theme(legend.position = c(0.15, 0.75), legend.title=element_blank(), axis.title.x = element_blank(),
		  axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
	annotate("text", x = 5, y = 67, label = "A", size = 5)




##########
###### 8. Packaging output for mapping
##########

head(suit_const_summary)

# suitabilityOutputForMapping <- data.frame(temp = suit_const_summary_rel$temperature,
# 										  constant = suit_const_summary_rel$median,
# 										  empDTR09 = suit_empfluc_dtr09_summary_rel$median,
# 										  empDTR12 = suit_empfluc_dtr12_summary_rel$median,
# 										  rs3TraitsDTR09 = suit_rs_3traits_09_summary_rel$median,
# 										  rs3TraitsDTR12 = suit_rs_3traits_12_summary_rel$median,
# 										  rsAllTraitsDTR09 = suit_rs_traitsall_09_summary_rel$median,
# 										  rsAllTraitsDTR12 = suit_rs_traitsall_12_summary_rel$median,
# 										  rsSuitDTR09 = suit_rs_suit_09_summary_rel$median,
# 										  rsSuitDTR12 = suit_rs_suit_12_summary_rel$median)

# Scale each model's lower CI contour to itself get relative R0/S(T) for mapping
suitabilityOutputForMapping_lowerCI <- data.frame(temp = suit_const_summary$temperature,
												   constant = suit_const_summary$lowerCI/max(suit_const_summary$lowerCI),
												   empDTR12 = suit_empfluc_dtr12_summary$lowerCI/max(suit_empfluc_dtr12_summary$lowerCI),
												   rsAllTraitsDTR12 = suit_rs_traitsall_12_summary$lowerCI/max(suit_rs_traitsall_12_summary$lowerCI),
												   rsSuitDTR12 = suit_rs_suit_12_summary$lowerCI/max(suit_rs_suit_12_summary$lowerCI))


#Check output
plot(constant ~ temp, data = suitabilityOutputForMapping_lowerCI, type = "l")
lines(empDTR09 ~ temp, data = suitabilityOutputForMapping_lowerCI, col = "red")
lines(empDTR12 ~ temp, data = suitabilityOutputForMapping_lowerCI, col = "darkorange")
lines(rs3TraitsDTR09 ~ temp, data = suitabilityOutputForMapping_lowerCI, col = "green")
lines(rs3TraitsDTR12 ~ temp, data = suitabilityOutputForMapping_lowerCI, col = "darkgreen")
lines(rsAllTraitsDTR09 ~ temp, data = suitabilityOutputForMapping_lowerCI, col = "dodgerblue")
lines(rsAllTraitsDTR12 ~ temp, data = suitabilityOutputForMapping_lowerCI, col = "blue")
lines(rsSuitDTR09 ~ temp, data = suitabilityOutputForMapping_lowerCI, col = "violet")
lines(rsSuitDTR12 ~ temp, data = suitabilityOutputForMapping_lowerCI, col = "purple")

#Check output
plot(constant ~ temp, data = suitabilityOutputForMapping, type = "l")
lines(empDTR09 ~ temp, data = suitabilityOutputForMapping, col = "red")
lines(empDTR12 ~ temp, data = suitabilityOutputForMapping, col = "darkorange")
lines(rs3TraitsDTR09 ~ temp, data = suitabilityOutputForMapping, col = "green")
lines(rs3TraitsDTR12 ~ temp, data = suitabilityOutputForMapping, col = "darkgreen")
lines(rsAllTraitsDTR09 ~ temp, data = suitabilityOutputForMapping, col = "dodgerblue")
lines(rsAllTraitsDTR12 ~ temp, data = suitabilityOutputForMapping, col = "blue")
lines(rsSuitDTR09 ~ temp, data = suitabilityOutputForMapping, col = "violet")
lines(rsSuitDTR12 ~ temp, data = suitabilityOutputForMapping, col = "purple")

write.csv(suitabilityOutputForMapping, "data-processed/RateSummationProjectR0forMapping.csv")
write.csv(suitabilityOutputForMapping_lowerCI, "data-processed/RateSummationProjectR0forMapping_lowerCI.csv")

suitabilityOutputForMapping.test <- read.csv("data-processed/RateSummationProjectR0forMapping.csv")
