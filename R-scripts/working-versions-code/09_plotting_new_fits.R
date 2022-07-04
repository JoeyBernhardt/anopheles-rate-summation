
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(patchwork)




# New plotting centralized ------------------------------------------------


# ### bc data
# "data-processed/predictions_bc_constant_summary.csv"
# "data-processed/params_bc_constant_all.csv"
# "data-processed/data_bc_constant_sum.csv"

### bite rate constant
bite.rate.constant.predictions <- read_csv("data-processed/predictions_bite_rate_constant_summary.csv")
bite.rate.constant.params <- read_csv("data-processed/params_bite_rate_constant_all.csv")
bite.rate.constant.sum <- read_csv("data-processed/data_bite_rate_constant_sum.csv")

bite.rate.dtr9.predictions <- read_csv("data-processed/predictions_bite_rate_dtr9_summary.csv")
bite.rate.dtr9.params <- read_csv("data-processed/params_bite_rate_dtr9_all.csv")
bite.rate.dtr9.sum <- read_csv("data-processed/data_bite_rate_dtr9_sum.csv")

bite.rate.dtr12.predictions <- read_csv("data-processed/predictions_bite_rate_dtr12_summary.csv")
bite.rate.dtr12.params <- read_csv("data-processed/params_bite_rate_dtr12_all.csv")
bite.rate.dtr12.sum <- read_csv("data-processed/data_bite_rate_dtr12_sum.csv")

#### all bite rates


all_bite_rate_predictions <- bind_rows(bite.rate.constant.predictions, bite.rate.dtr9.predictions, bite.rate.dtr12.predictions)
all_bite_rate_sums <- bind_rows(bite.rate.constant.sum, bite.rate.dtr9.sum, bite.rate.dtr12.sum)


### plot for bite rate
bite_rate_plot <- all_bite_rate_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = `2.5%`, ymax = `97.5%`, fill = treatment), alpha = 0.6) +
	geom_line(aes(x = temperature, y = mean, color = treatment)) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean), shape = 1, data = all_bite_rate_sums, size = 1) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean, color = treatment), data = all_bite_rate_sums, size = 0.8) +
	scale_color_viridis_d(begin = 0.1, end = 0.8) +
	scale_fill_viridis_d(begin = 0.1, end = 0.8) + ylab("Bite rate") + xlab("Temperature (C)")

ggsave("figures/bite_rate_plot.pdf", width = 8, height = 6)

### now get the parameters on there.

all_params_bite_rate <- bind_rows(bite.rate.constant.params, bite.rate.dtr9.params,bite.rate.dtr12.params) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))

params_plot_bite_rate <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = `2.5%`, ymax = `97.5%`, color = treatment), data = all_params_bite_rate, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_viridis_d(begin = 0.1, end = 0.8) +
	scale_fill_viridis_d(begin = 0.1, end = 0.8)

bite_rate_plot_all <- bite_rate_plot / params_plot_bite_rate + plot_layout(heights = c(3, 0.5))

ggsave("figures/bite_rate_plot_all.pdf", bite_rate_plot_all, width = 7, height = 5)


# lifespan ---------------------------------------------------------


lifespan.constant.predictions <- read_csv("data-processed/predictions_lifespan_constant_summary_quad.csv")
lifespan.constant.params <- read_csv("data-processed/params_lifespan_constant_all_quad.csv")
lifespan.constant.sum <- read_csv("data-processed/data_lifespan_constant_sum_quad.csv")

lifespan.dtr9.predictions <- read_csv("data-processed/predictions_lifespan_dtr9_summary_quad.csv")
lifespan.dtr9.params <- read_csv("data-processed/params_lifespan_dtr9_all_quad.csv")
lifespan.dtr9.sum <- read_csv("data-processed/data_lifespan_dtr9_sum_quad.csv")

lifespan.dtr12.predictions <- read_csv("data-processed/predictions_lifespan_dtr12_summary_quad.csv")
lifespan.dtr12.params <- read_csv("data-processed/params_lifespan_dtr12_all_quad.csv")
lifespan.dtr12.sum <- read_csv("data-processed/data_lifespan_dtr12_sum_quad.csv")

#### all lifespan

all_lifespan_predictions <- bind_rows(lifespan.constant.predictions, lifespan.dtr9.predictions, lifespan.dtr12.predictions)
all_lifespan_sums <- bind_rows(lifespan.constant.sum, lifespan.dtr9.sum, lifespan.dtr12.sum)


# View(lifespan.constant.sum )

### plot for lifespan
lifespan_plot <- all_lifespan_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = `2.5%`, ymax = `97.5%`, fill = treatment), alpha = 0.6) +
	geom_line(aes(x = temperature, y = mean, color = treatment)) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean), shape = 1, data = all_lifespan_sums, size = 1) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean, color = treatment), data = all_lifespan_sums, size = 0.8) +
	scale_color_viridis_d(begin = 0.1, end = 0.8) +
	scale_fill_viridis_d(begin = 0.1, end = 0.8) + ylab("lifespan") + xlab("Temperature (C)") + xlim(0,40)

ggsave("figures/lifespan_plot.pdf", width = 8, height = 6)

### now get the parameters on there.

all_params_lifespan <- bind_rows(lifespan.constant.params, lifespan.dtr9.params,lifespan.dtr12.params) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))
params_plot_lifespan <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = `2.5%`, ymax = `97.5%`, color = treatment), data = all_params_lifespan, position=position_dodge(width=1)) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_viridis_d(begin = 0.1, end = 0.8) +
	scale_fill_viridis_d(begin = 0.1, end = 0.8) + ylim(0,40)

lifespan_plot_all <- lifespan_plot / params_plot_lifespan + plot_layout(heights = c(3, 0.5))

ggsave("figures/lifespan_plot_all.pdf", lifespan_plot_all, width = 7, height = 5)

# eggs ----------------------------------------------------------------

eggs.constant.predictions <- read_csv("data-processed/predictions_eggs_constant_summary.csv")
eggs.constant.params <- read_csv("data-processed/params_constant_eggs_all.csv")
eggs.constant.sum <- read_csv("data-processed/data_constant_eggs_sum.csv")

eggs.dtr9.predictions <- read_csv("data-processed/predictions_eggs_dtr9_summary.csv")
eggs.dtr9.params <- read_csv("data-processed/params_dtr9_eggs_all.csv")
eggs.dtr9.sum <- read_csv("data-processed/data_dtr9_eggs_sum.csv")

eggs.dtr12.predictions <- read_csv("data-processed/predictions_eggs_dtr12_summary.csv")
eggs.dtr12.params <- read_csv("data-processed/params_dtr12_eggs_all.csv")
eggs.dtr12.sum <- read_csv("data-processed/data_dtr12_eggs_sum.csv")

#### all eggs

all_eggs_predictions <- bind_rows(eggs.constant.predictions, eggs.dtr9.predictions, eggs.dtr12.predictions)
all_eggs_sums <- bind_rows(eggs.constant.sum, eggs.dtr9.sum, eggs.dtr12.sum)


### plot for eggs
eggs_plot <- all_eggs_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = `2.5%`, ymax = `97.5%`, fill = treatment), alpha = 0.6) +
	geom_line(aes(x = temperature, y = mean, color = treatment)) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean), shape = 1, data = all_eggs_sums, size = 1) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean, color = treatment), data = all_eggs_sums, size = 0.8) +
	scale_color_viridis_d(begin = 0.1, end = 0.8) +
	scale_fill_viridis_d(begin = 0.1, end = 0.8) + ylab("eggs") + xlab("Temperature (C)") +xlim(0, 40)

ggsave("figures/eggs_plot.pdf", width = 8, height = 6)

### now get the parameters on there.

all_params_eggs <- bind_rows(eggs.constant.params, eggs.dtr9.params,eggs.dtr12.params) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))
params_plot_eggs <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = `2.5%`, ymax = `97.5%`, color = treatment), data = all_params_eggs, position=position_dodge(width=1)) +
	ylim(0, 40) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_viridis_d(begin = 0.1, end = 0.8) +
	scale_fill_viridis_d(begin = 0.1, end = 0.8)
	

eggs_plot_all <- eggs_plot / params_plot_eggs + plot_layout(heights = c(3, 0.5))

ggsave("figures/eggs_plot_all.pdf", eggs_plot_all, width = 7, height = 5)




# pea ---------------------------------------------------------------------
pea.constant.predictions <- read_csv("data-processed/predictions_pea_summary.csv")
pea.constant.params <- read_csv("data-processed/params_pea_all.csv")
pea.constant.sum <- read_csv("data-processed/data_pea_sum.csv")

#### all pea

all_pea_predictions <- bind_rows(pea.constant.predictions)
all_pea_sums <- bind_rows(pea.constant.sum)


### plot for pea
pea_plot <- all_pea_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = `2.5%`, ymax = `97.5%`, fill = treatment), alpha = 0.6) +
	geom_line(aes(x = temperature, y = mean, color = treatment)) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean), shape = 1, data = all_pea_sums, size = 1) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean, color = treatment), data = all_pea_sums, size = 0.8) +
	scale_color_viridis_d(begin = 0.1, end = 0.8) +
	scale_fill_viridis_d(begin = 0.1, end = 0.8) + ylab("pea") + xlab("Temperature (C)") +xlim(0, 40)

ggsave("figures/pea_plot.pdf", width = 8, height = 6)

### now get the parameters on there.

all_params_pea <- bind_rows(pea.constant.params) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))
params_plot_pea <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = `2.5%`, ymax = `97.5%`, color = treatment), data = all_params_pea, position=position_dodge(width=1)) +
	ylim(0, 40) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_viridis_d(begin = 0.1, end = 0.8) +
	scale_fill_viridis_d(begin = 0.1, end = 0.8)


pea_plot_all <- pea_plot / params_plot_pea + plot_layout(heights = c(3, 0.5))

ggsave("figures/pea_plot_all.pdf", pea_plot_all, width = 7, height = 5)


# mdr ---------------------------------------------------------------------
mdr.constant.predictions <- read_csv("data-processed/predictions_mdr_summary.csv")
mdr.constant.params <- read_csv("data-processed/params_mdr_all.csv")
mdr.constant.sum <- read_csv("data-processed/data_mdr_sum.csv")

#### all mdr

all_mdr_predictions <- bind_rows(mdr.constant.predictions)
all_mdr_sums <- bind_rows(mdr.constant.sum)


### plot for mdr
mdr_plot <- all_mdr_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = `2.5%`, ymax = `97.5%`, fill = treatment), alpha = 0.6) +
	geom_line(aes(x = temperature, y = mean, color = treatment)) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean), shape = 1, data = all_mdr_sums, size = 1) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean, color = treatment), data = all_mdr_sums, size = 0.8) +
	scale_color_viridis_d(begin = 0.1, end = 0.8) +
	scale_fill_viridis_d(begin = 0.1, end = 0.8) + ylab("mdr") + xlab("Temperature (C)") +xlim(0, 40)

ggsave("figures/mdr_plot.pdf", width = 8, height = 6)

### now get the parameters on there.

all_params_mdr <- bind_rows(mdr.constant.params) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))
params_plot_mdr <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = `2.5%`, ymax = `97.5%`, color = treatment), data = all_params_mdr, position=position_dodge(width=1)) +
	ylim(0, 40) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_viridis_d(begin = 0.1, end = 0.8) +
	scale_fill_viridis_d(begin = 0.1, end = 0.8)


mdr_plot_all <- mdr_plot / params_plot_mdr + plot_layout(heights = c(3, 0.5))

ggsave("figures/mdr_plot_all.pdf", mdr_plot_all, width = 7, height = 5)



# eip50 ---------------------------------------------------------------------
eip50.constant.predictions <- read_csv("data-processed/predictions_eip50_summary.csv")
eip50.constant.params <- read_csv("data-processed/params_eip50_all.csv")
eip50.constant.sum <- read_csv("data-processed/data_eip50_sum.csv")

#### all eip50

all_eip50_predictions <- bind_rows(eip50.constant.predictions)
all_eip50_sums <- bind_rows(eip50.constant.sum)


### plot for eip50
eip50_plot <- all_eip50_predictions %>% 
	ggplot() +
	geom_ribbon(aes(x = temperature, ymin = `2.5%`, ymax = `97.5%`, fill = treatment), alpha = 0.6) +
	geom_line(aes(x = temperature, y = mean, color = treatment)) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean), shape = 1, data = all_eip50_sums, size = 1) +
	geom_pointrange(aes(x = temperature, ymin = mean - std_error, ymax = mean + std_error, y = mean, color = treatment), data = all_eip50_sums, size = 0.8) +
	scale_color_viridis_d(begin = 0.1, end = 0.8) +
	scale_fill_viridis_d(begin = 0.1, end = 0.8) + ylab("eip50") + xlab("Temperature (C)") +xlim(0, 45)

ggsave("figures/eip50_plot.pdf", width = 8, height = 6)

### now get the parameters on there.

all_params_eip50 <- bind_rows(eip50.constant.params) %>% 
	dplyr::filter(term %in% c("Topt", "cf.Tm", "cf.T0")) %>% 
	mutate(term = case_when(term == "cf.Tm" ~ "Tmax",
							term == "cf.T0"~ "Tmin",
							term == "Topt" ~ "Topt"))
params_plot_eip50 <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = `2.5%`, ymax = `97.5%`, color = treatment), data = all_params_eip50, position=position_dodge(width=1)) +
	ylim(0, 40) +
	coord_flip() + ylab("Temperature (°C)") +
	scale_color_viridis_d(begin = 0.1, end = 0.8) +
	scale_fill_viridis_d(begin = 0.1, end = 0.8)


eip50_plot_all <- eip50_plot / params_plot_eip50 + plot_layout(heights = c(3, 0.5))

ggsave("figures/eip50_plot_all.pdf", eip50_plot_all, width = 7, height = 5)

### ok now merge all three plots together! woohoo!


all_plots <- wrap_plots(eggs_plot_all, lifespan_plot_all, bite_rate_plot_all, eip50_plot_all, mdr_plot_all, pea_plot_all)
ggsave('figures/all_plots_combined-feb4-3.pdf', all_plots, width = 18, height = 6)
ggsave('figures/all_plots_combined-march22.pdf', all_plots, width = 18, height = 6)
### OK so the best fits are briere for bite rate and egg production but quadratic for lifespan.


