## Rate Summation Project
## Script for plotting TPC fits
## Written by Marta Shocket & Joey Bernhardt in 2022

##	Table of contents:
##
##	1. Set-up workspace
##	2. Load and process data for Bite Rate (a), Lifespan (lf), and Lifetime eggs (B)
##	3. Plot panels for Bite Rate (a), Lifespan (lf), and Lifetime eggs (B)
##	4. Manuscript Figure 1
##	5. Load and process data for other traits (bc, EIP50, pEA, MDR, gamma)
##	6. Plot panels for other traits (bc, EIP50, pEA, MDR, gamma)
##	7. Manuscript Figure S1



##########
###### 1. Set-up workspace
##########

### Load libraries
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(patchwork)

##### Load functions
source("R-scripts/working-versions-code/00_RSProjectFunctions.R")


##########
###### 2. Load and process data for Bite Rate (a), Lifespan (lf), and Lifetime eggs (B)
##########

####################### Bite rate (a)

### load data - all fluctuation treatments
predictions_bite_rate_constant_summary <- read.csv("data-processed/predictions_bite_rate_constant_summary.csv")
params_bite_rate_constant_summary <- read.csv("data-processed/params_bite_rate_constant_summary.csv")
data_bite_rate_constant_summary <- read.csv("data-processed/data_bite_rate_constant_summary.csv")

predictions_bite_rate_dtr9_summary <- read.csv("data-processed/predictions_bite_rate_dtr9_summary.csv")
params_bite_rate_dtr9_summary <- read.csv("data-processed/params_bite_rate_dtr9_summary.csv")
data_bite_rate_dtr9_summary <- read.csv("data-processed/data_bite_rate_dtr9_summary.csv")

predictions_bite_rate_dtr12_summary <- read.csv("data-processed/predictions_bite_rate_dtr12_summary.csv")
params_bite_rate_dtr12_summary <- read.csv("data-processed/params_bite_rate_dtr12_summary.csv")
data_bite_rate_dtr12_summary <- read.csv("data-processed/data_bite_rate_dtr12_summary.csv")

### edit DTR9 treatment labels so they're alphabetically before DTR12 - ggplot2 assigns categorical variables colors in alphabetical order
predictions_bite_rate_dtr9_summary$treatment <- "bite_rate_dtr09"
params_bite_rate_dtr9_summary$treatment <- "bite_rate_dtr09"
data_bite_rate_dtr9_summary$treatment <- "bite_rate_dtr09"

### combine all fluctuation treatments
all_bite_rate_predictions <- bind_rows(predictions_bite_rate_constant_summary, predictions_bite_rate_dtr9_summary, predictions_bite_rate_dtr12_summary)
all_bite_rate_sums <- bind_rows(data_bite_rate_constant_summary, data_bite_rate_dtr9_summary, data_bite_rate_dtr12_summary)

### combine and merge parameters with different names 
all_bite_rate_params <- bind_rows(params_bite_rate_constant_summary, params_bite_rate_dtr9_summary, params_bite_rate_dtr12_summary) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

####################### Lifespan (lf)

### load data - all fluctuation treatments
predictions_lifespan_constant_summary <- read.csv("data-processed/predictions_lifespan_constant_summary.csv")
params_lifespan_constant_summary <- read.csv("data-processed/params_lifespan_constant_summary.csv")
data_lifespan_constant_summary <- read.csv("data-processed/data_lifespan_constant_summary.csv")

predictions_lifespan_dtr9_summary <- read.csv("data-processed/predictions_lifespan_dtr9_summary.csv")
params_lifespan_dtr9_summary <- read.csv("data-processed/params_lifespan_dtr9_summary.csv")
data_lifespan_dtr9_summary <- read.csv("data-processed/data_lifespan_dtr9_summary.csv")

predictions_lifespan_dtr12_summary <- read.csv("data-processed/predictions_lifespan_dtr12_summary.csv")
params_lifespan_dtr12_summary <- read.csv("data-processed/params_lifespan_dtr12_summary.csv")
data_lifespan_dtr12_summary <- read.csv("data-processed/data_lifespan_dtr12_summary.csv")

### edit DTR9 treatment labels so they're alphabetically before DTR12 - ggplot2 assigns categorical variables colors in alphabetical order
predictions_lifespan_dtr9_summary$treatment <- "lifespan_dtr09"
params_lifespan_dtr9_summary$treatment <- "lifespan_dtr09"
data_lifespan_dtr9_summary$treatment <- "lifespan_dtr09"

### combine all fluctuation treatments
all_lifespan_predictions <- bind_rows(predictions_lifespan_constant_summary, predictions_lifespan_dtr9_summary, predictions_lifespan_dtr12_summary)
all_lifespan_sums <- bind_rows(data_lifespan_constant_summary, data_lifespan_dtr9_summary, data_lifespan_dtr12_summary)

### combine and merge different names 
all_lifespan_params <- bind_rows(params_lifespan_constant_summary, params_lifespan_dtr9_summary, params_lifespan_dtr12_summary) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

####################### Lifetime eggs (B)

### load data - all fluctuation treatments
predictions_eggs_constant_summary <- read.csv("data-processed/predictions_eggs_constant_summary.csv")
params_eggs_constant_summary <- read.csv("data-processed/params_eggs_constant_summary.csv")
data_eggs_constant_summary <- read.csv("data-processed/data_eggs_constant_summary.csv")

predictions_eggs_dtr9_summary <- read.csv("data-processed/predictions_eggs_dtr9_summary.csv")
params_eggs_dtr9_summary <- read.csv("data-processed/params_eggs_dtr9_summary.csv")
data_eggs_dtr9_summary <- read.csv("data-processed/data_eggs_dtr9_summary.csv")

predictions_eggs_dtr12_summary <- read.csv("data-processed/predictions_eggs_dtr12_summary.csv")
params_eggs_dtr12_summary <- read.csv("data-processed/params_eggs_dtr12_summary.csv")
data_eggs_dtr12_summary <- read.csv("data-processed/data_eggs_dtr12_summary.csv")

### edit DTR9 treatment labels so they're alphabetically before DTR12 - ggplot2 assigns categorical variables colors in alphabetical order
predictions_eggs_dtr9_summary$treatment <- "eggs_dtr09"
params_eggs_dtr9_summary$treatment <- "eggs_dtr09"
data_eggs_dtr9_summary$treatment <- "eggs_dtr09"

#### combine all fluctuation treatments
all_eggs_predictions <- bind_rows(predictions_eggs_constant_summary, predictions_eggs_dtr9_summary, predictions_eggs_dtr12_summary)
all_eggs_sums <- bind_rows(data_eggs_constant_summary, data_eggs_dtr9_summary, data_eggs_dtr12_summary)

### combine and merge different names 
all_eggs_params <- bind_rows(params_eggs_constant_summary, params_eggs_dtr9_summary, params_eggs_dtr12_summary) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

##########
###### 3. Plot panels for Bite Rate (a), Lifespan (lf), and Lifetime eggs (B)
##########

### plot for bite rate - xlim() gives warning because not all rows are used in the plot
bite_rate_plot <- all_bite_rate_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment), size = 0.6) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean), shape = 1, data = all_bite_rate_sums, size = 1) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean, color = treatment), data = all_bite_rate_sums, size = 0.8) +
	scale_color_manual(values = c(c_constant, c_emp09, c_emp12), labels = c("Constant", "Empirical DTR 9", "Empirical DTR 12")) +
	scale_fill_manual(values = c(ct_constant, ct_emp09, ct_emp12), labels = c("Constant", "Empirical DTR 9", "Empirical DTR 12")) +
	ylab(parse(text = "Bite~rate~(day^-1)~-~a")) + xlab("Temperature (°C)") +
	theme_classic() +
	theme(axis.title.x = element_blank(), legend.position = c(0.2, 0.8), legend.title=element_blank(), legend.text = element_text(size=11)) + 
	annotate("text", x = 0, y = 0.53, label = "A", size = 5)

params_plot_bite_rate <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = all_bite_rate_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_manual(values = c(c_constant, c_emp09, c_emp12), labels = c("Constant", "DTR 9", "DTR 12")) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 0, label = "D", size = 5)

bite_rate_plot_all <- bite_rate_plot / params_plot_bite_rate + plot_layout(heights = c(3, 0.5))

### plot for lifespan
lifespan_plot <- all_lifespan_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment), size = 0.6) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean), shape = 1, data = all_lifespan_sums, size = 1) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean, color = treatment), data = all_lifespan_sums, size = 0.8) +
	scale_color_manual(values = c(c_constant, c_emp09, c_emp12), labels = c("Constant", "DTR 9", "DTR 12")) +
	scale_fill_manual(values = c(ct_constant, ct_emp09, ct_emp12), labels = c("Constant", "DTR 9", "DTR 12")) +
	ylab("Lifespan (days) - lf") + xlab("Temperature (°C)") + xlim(0,40) +
	theme_classic() +
	theme(axis.title.x = element_blank(), legend.position = "none") +
	annotate("text", x = 0, y = 51, label = "B", size = 5)

params_plot_lifespan <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = all_lifespan_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_manual(values = c(c_constant, c_emp09, c_emp12), labels = c("Constant", "DTR 9", "DTR 12")) + ylim(0,40) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 0, label = "E", size = 5)

lifespan_plot_all <- lifespan_plot / params_plot_lifespan + plot_layout(heights = c(3, 0.5))

### plot for eggs
eggs_plot <- all_eggs_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment)) +
	geom_line(aes(x = temperature, y = mean, color = treatment), size = 0.6) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean), shape = 1, data = all_eggs_sums, size = 1) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean, color = treatment), data = all_eggs_sums, size = 0.8) +
	scale_color_manual(values = c(c_constant, c_emp09, c_emp12), labels = c("Constant", "DTR 9", "DTR 12")) +
	scale_fill_manual(values = c(ct_constant, ct_emp09, ct_emp12), labels = c("Constant", "DTR 9", "DTR 12")) +
	ylab("Lifetime eggs - B") + xlab("Temperature (°C)") +xlim(0, 40) +
	theme_classic() +
	theme(axis.title.x = element_blank(), legend.position = "none") +
	annotate("text", x = 0, y = 450, label = "C", size = 5)

params_plot_eggs <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = all_eggs_params, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_manual(values = c(c_constant, c_emp09, c_emp12), labels = c("Constant", "DTR 9", "DTR 12")) +
	theme_classic() +
	theme(legend.position = "none", axis.title.y = element_blank()) +
	annotate("text", x = 3.2, y = 0, label = "F", size = 5)

eggs_plot_all <- eggs_plot / params_plot_eggs + plot_layout(heights = c(3, 0.5))

##########
###### 4. Manuscript Figure 2
##########

Fig1_plots <- wrap_plots(bite_rate_plot_all, lifespan_plot_all, eggs_plot_all)
ggsave('figures/Fig1_plots-Sept2024.pdf', Fig1_plots, width = 15, height = 6)

#### Alternative Figure Panel Arrangements
#
# Fig1_plots_vert <- wrap_plots(ncol = 1, bite_rate_plot_all, lifespan_plot_all, eggs_plot_all)
# ggsave('figures/Fig1_plots_vert-Jan4.pdf', Fig1_plots_vert, width = 6, height = 17)
# 
# Fig1_plots_sq <- wrap_plots(ncol = 2, bite_rate_plot_all, lifespan_plot_all, eggs_plot_all)
# ggsave('figures/Fig1_plots_sq-Jan4.pdf', Fig1_plots_sq, width = 10, height = 11)
# 
# design <- "AAABBBB
# 		   #CCCC##"
# Fig1_plots_sq_c <- wrap_plots(ncol = 2, design = design, A = bite_rate_plot_all, B = lifespan_plot_all, C = eggs_plot_all)
# ggsave('figures/Fig1_plots_sq_c-Jan4.pdf', Fig1_plots_sq_c, width = 10, height = 11)



##########
###### 5. Load and process data for other traits (bc, EIP50, pEA, MDR, gamma)
##########

### load data
all_bc_predictions <- read.csv("data-processed/predictions_bc_summary.csv")
bc.constant.params <- read.csv("data-processed/params_bc_all.csv")
all_bc_sums <- read.csv("data-processed/data_bc_summary.csv")

### plot for bc
bc_plot <- all_bc_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment), alpha = 0.6) +
	geom_line(aes(x = temperature, y = mean, color = treatment)) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean), shape = 1, data = all_bc_sums, size = 1) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean, color = treatment), data = all_bc_sums, size = 0.8) +
	scale_color_viridis_d(begin = 0.1, end = 0.8) +
	scale_fill_viridis_d(begin = 0.1, end = 0.8) + ylab("bc") + xlab("Temperature (°C)") +xlim(0, 45)

# ggsave("figures/bc_plot.pdf", width = 8, height = 6)

### now get the parameters on there
all_params_bc <- bind_rows(bc.constant.params) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

params_plot_bc <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = all_params_bc, position=position_dodge(width=1)) +
	ylim(0, 45) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_viridis_d(begin = 0.1, end = 0.8) +
	scale_fill_viridis_d(begin = 0.1, end = 0.8)


bc_plot_all <- bc_plot / params_plot_bc + plot_layout(heights = c(3, 0.5))

# ggsave("figures/bc_plot_all.pdf", bc_plot_all, width = 7, height = 5)


##########
###### 7. Extrinsic Incubuation Period 50% (EIP50)  ----------------------------------------------------------
##########

### load data
all_eip50_predictions <- read.csv("data-processed/predictions_eip50_summary.csv")
eip50.constant.params <- read.csv("data-processed/params_eip50_all.csv")
all_eip50_sums <- read.csv("data-processed/data_eip50_summary.csv")

### plot for eip50
eip50_plot <- all_eip50_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment), alpha = 0.6) +
	geom_line(aes(x = temperature, y = mean, color = treatment)) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean), shape = 1, data = all_eip50_sums, size = 1) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean, color = treatment), data = all_eip50_sums, size = 0.8) +
	scale_color_viridis_d(begin = 0.1, end = 0.8) +
	scale_fill_viridis_d(begin = 0.1, end = 0.8) + ylab("eip50") + xlab("Temperature (°C)") +xlim(0, 45)

# ggsave("figures/eip50_plot.pdf", width = 8, height = 6)

### now get the parameters on there
all_params_eip50 <- bind_rows(eip50.constant.params) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

params_plot_eip50 <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = all_params_eip50, position=position_dodge(width=1)) +
	ylim(0, 45) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_viridis_d(begin = 0.1, end = 0.8) +
	scale_fill_viridis_d(begin = 0.1, end = 0.8)


eip50_plot_all <- eip50_plot / params_plot_eip50 + plot_layout(heights = c(3, 0.5))

# ggsave("figures/eip50_plot_all.pdf", eip50_plot_all, width = 7, height = 5)



##########
###### 8. Juvenile survival (pea) ----------------------------------------------------------
##########

### Load data
all_pea_predictions <- read.csv("data-processed/predictions_pea_summary.csv")
pea.constant.params <- read.csv("data-processed/params_pea_all.csv")
all_pea_sums <- read.csv("data-processed/data_pea_summary.csv")

### plot for pea
pea_plot <- all_pea_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment), alpha = 0.6) +
	geom_line(aes(x = temperature, y = mean, color = treatment)) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean), shape = 1, data = all_pea_sums, size = 1) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean, color = treatment), data = all_pea_sums, size = 0.8) +
	scale_color_viridis_d(begin = 0.1, end = 0.8) +
	scale_fill_viridis_d(begin = 0.1, end = 0.8) + ylab("pea") + xlab("Temperature (°C)") +xlim(0, 40)

# ggsave("figures/pea_plot.pdf", width = 8, height = 6)

### now get the parameters on there
all_params_pea <- bind_rows(pea.constant.params) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

params_plot_pea <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = all_params_pea, position=position_dodge(width=1)) +
	ylim(0, 40) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_viridis_d(begin = 0.1, end = 0.8) +
	scale_fill_viridis_d(begin = 0.1, end = 0.8)


pea_plot_all <- pea_plot / params_plot_pea + plot_layout(heights = c(3, 0.5))

# ggsave("figures/pea_plot_all.pdf", pea_plot_all, width = 7, height = 5)



##########
###### 9. Mosquito development rate (MDR) ----------------------------------------------------------
##########

### load data
all_mdr_predictions <- read.csv("data-processed/predictions_mdr_summary.csv")
mdr.constant.params <- read.csv("data-processed/params_mdr_all.csv")
all_mdr_sums <- read.csv("data-processed/data_mdr_summary.csv")

### plot for mdr
mdr_plot <- all_mdr_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment), alpha = 0.6) +
	geom_line(aes(x = temperature, y = mean, color = treatment)) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean), shape = 1, data = all_mdr_sums, size = 1) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean, color = treatment), data = all_mdr_sums, size = 0.8) +
	scale_color_viridis_d(begin = 0.1, end = 0.8) +
	scale_fill_viridis_d(begin = 0.1, end = 0.8) + ylab("mdr") + xlab("Temperature (°C)") +xlim(0, 40)

# ggsave("figures/mdr_plot.pdf", width = 8, height = 6)

### now get the parameters on there
all_params_mdr <- bind_rows(mdr.constant.params) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

params_plot_mdr <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = all_params_mdr, position=position_dodge(width=1)) +
	ylim(0, 40) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_viridis_d(begin = 0.1, end = 0.8) +
	scale_fill_viridis_d(begin = 0.1, end = 0.8)


mdr_plot_all <- mdr_plot / params_plot_mdr + plot_layout(heights = c(3, 0.5))

# ggsave("figures/mdr_plot_all.pdf", mdr_plot_all, width = 7, height = 5)



##########
###### 10. Gamma ----------------------------------------------------------
##########

### load data - separate fluctuation treatments
gamma.constant.predictions <- read.csv("data-processed/predictions_gamma_constant_summary.csv")
gamma.constant.params <- read.csv("data-processed/params_gamma_constant_all.csv")
gamma.constant.sum <- read.csv("data-processed/data_gamma_constant_sum.csv")

gamma.dtr9.predictions <- read.csv("data-processed/predictions_gamma_dtr9_summary.csv")
gamma.dtr9.params <- read.csv("data-processed/params_gamma_dtr9_all.csv")
gamma.dtr9.sum <- read.csv("data-processed/data_gamma_dtr9_sum.csv")

gamma.dtr12.predictions <- read.csv("data-processed/predictions_gamma_dtr12_summary.csv")
gamma.dtr12.params <- read.csv("data-processed/params_gamma_dtr12_all.csv")
gamma.dtr12.sum <- read.csv("data-processed/data_gamma_dtr12_sum.csv")

#### combine all fluctuation treatments
all_gamma_predictions <- bind_rows(gamma.constant.predictions, gamma.dtr9.predictions, gamma.dtr12.predictions)
all_gamma_sums <- bind_rows(gamma.constant.sum, gamma.dtr9.sum, gamma.dtr12.sum)


### plot for gamma
gamma_plot <- all_gamma_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = lowerCI, ymax = upperCI, fill = treatment), alpha = 0.6) +
	geom_line(aes(x = temperature, y = mean, color = treatment)) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean), shape = 1, data = all_gamma_sums, size = 1) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean, color = treatment), data = all_gamma_sums, size = 0.8) +
	scale_color_viridis_d(begin = 0.1, end = 0.8) +
	scale_fill_viridis_d(begin = 0.1, end = 0.8) + ylab("gamma") + xlab("Temperature (°C)") + xlim(0, 45)

# ggsave("figures/gamma_plot.pdf", width = 8, height = 6)

### now get the parameters on there
all_params_gamma <- bind_rows(gamma.constant.params, gamma.dtr9.params,gamma.dtr12.params) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

params_plot_gamma <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = lowerCI, ymax = upperCI, color = treatment), data = all_params_gamma, position=position_dodge(width=1)) +
	ylim(0, 45) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_viridis_d(begin = 0.1, end = 0.8) +
	scale_fill_viridis_d(begin = 0.1, end = 0.8)

gamma_plot_all <- gamma_plot / params_plot_gamma + plot_layout(heights = c(3, 0.5))

# ggsave("figures/gamma_plot_all.pdf", gamma_plot_all, width = 7, height = 5)



##########
###### 11. Manuscript Figure S1 ----------------------------------------------------------
##########

other_plots <- wrap_plots(bc_plot, eip50_plot, mdr_plot, pea_plot,
							nrow = 2, ncol = 2)

all_plots <- wrap_plots(eggs_plot_all, lifespan_plot_all, bite_rate_plot_all)
all_plots <- wrap_plots(eggs_plot, lifespan_plot, bite_rate_plot)
all_plots <- wrap_plots(eggs_plot_all, lifespan_plot_all, bite_rate_plot_all, gamma_plot_all,
						nrow = 2, ncol = 2)

ggsave('figures/all_plots_combined-feb4-3.pdf', all_plots, width = 18, height = 6)
ggsave('figures/all_plots_combined-march22.pdf', all_plots, width = 18, height = 6)
ggsave('figures/all_plots_combined-aug12.pdf', all_plots, width = 18, height = 6)
ggsave('figures/all_plots_combined-aug12.pdf', all_plots, width = 18, height = 6)
