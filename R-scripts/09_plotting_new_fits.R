
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(patchwork)


#### Bite rate constant

load("saved-posteriors/constant_bite.rate_briere_uniform.Rdata")
bite.rate.constant.preds <-model.out$BUGSoutput$sims.list$z.trait.mu.pred


bite.rate.constant.preds <- model_out_constant_briere$BUGSoutput$sims.list$z.trait.mu.pred



b_params_bite_rate_briere_constant <- as.data.frame(model_out_constant_briere$BUGSoutput$summary[1:5,]) %>% 
	rownames_to_column(var = "term")


predictions_bite_rate_briere_constant <- as.data.frame(model_out_constant_briere$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)  ### columns are temperatures, rows are iterations
colnames(predictions_bite_rate_briere_constant) <- Temp.xs


topt_constant_bite_rate <- predictions_bite_rate_briere_constant %>% 
	mutate(iteration = rownames(.)) %>%
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>% 
	group_by(iteration) %>% 
	top_n(n = 1, wt = growth_rate) %>% 
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(q2.5=quantile(temperature, probs=0.025),
			  q97.5=quantile(temperature, probs=0.975),
			  mean = mean(temperature)) %>% 
	mutate(term = "Topt")


predictions_sub <- predictions_bite_rate_briere_constant %>% 
	mutate(iteration = rownames(.)) %>% 
	dplyr::select(iteration, everything()) %>% 
	sample_n(size = 1000)


predictions_long <- predictions_sub %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1))

predictions_summary <- predictions_bite_rate_briere_constant %>%
	mutate(iteration = rownames(.)) %>% 
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>%
	dplyr::group_by(temperature) %>%  
	summarise(q2.5=quantile(growth_rate, probs=0.025),
			  q97.5=quantile(growth_rate, probs=0.975),
			  mean = mean(growth_rate)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "bite_rate_constant")

View(predictions_summary)

library(beyonce)
(beyonce_palette(18, 3))


#66c2a5
#fc8d62
#8da0cb


bite_rate_constant_plot <- ggplot() + 
	geom_ribbon(aes(x = temperature, ymax = q97.5, ymin = q2.5), data = predictions_summary, fill = "#66c2a5") +
	# geom_line(aes(x = temperature, y = growth_rate, group = iteration), alpha = 0.05, size = 1, data = predictions_long2) +
	geom_point(aes(x = T, y = trait), data = data_constant_bite_rate_2, color = "darkgrey", alpha = 0.7, size = 2) + geom_hline(yintercept = 0) +
	geom_line(aes(x = temperature, y = mean), data = predictions_summary, color = "black") +
	geom_line(aes(x = temperature, y = q2.5), data = predictions_summary, color = "black", linetype = "dashed") +
	geom_line(aes(x = temperature, y = q97.5), data = predictions_summary, color = "black", linetype = "dashed") +
	
	# geom_line(aes(x = temperature - 273.15, y = mean), data = predictions_summary_briere, color = "blue") +
	# geom_line(aes(x = temperature - 273.15, y = q2.5), data = predictions_summary_briere, color = "blue", linetype = "dashed") +
	# geom_line(aes(x = temperature - 273.15, y = q97.5), data = predictions_summary_briere, color = "blue", linetype = "dashed") +
	ylab("Bite rate") + 
	# xlab("Temperature (°C)") +
	xlab("") +
	ylim(0, 1) + xlim(0, 45) 


params_plot_constant_bite_rate <- b_params_bite_rate_briere_constant %>% 
	filter(term %in% c("cf.T0", "cf.Tm")) %>% 
	dplyr::select(term, mean, `2.5%`, `97.5%`) %>% 
	# gather(key = term, value = value, col = c(2:4)) %>% View
	ggplot() + geom_pointrange(aes(x = term, y = mean, ymin = `2.5%`, ymax = `97.5%`)) +
	coord_flip() + ylab("Temperature (°C)") 





library(patchwork)


bite_rate_constant_plot / params_plot_constant_bite_rate +  plot_layout(heights = c(4, 1))


### December 16, come back to how to find Topt




# DTR9 bite rate ----------------------------------------------------------

model.out_bite_rate_dtr9_briere


bite.rate.dtr9.preds <- model.out_bite_rate_dtr9_briere$BUGSoutput$sims.list$z.trait.mu.pred


library(tidyverse)
b_params_dtr9_bite_rate <- as.data.frame(model.out_bite_rate_dtr9_briere$BUGSoutput$summary[1:5,]) %>% 
	rownames_to_column(var = "term")


predictions_dtr9_bite_rate <- as.data.frame(model.out_bite_rate_dtr9_briere$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)  ### columns are temperatures, rows are iterations
colnames(predictions_dtr9_bite_rate) <- Temp.xs





topt_dtr9_bite_rate <- predictions_dtr9_bite_rate %>% 
	mutate(iteration = rownames(.)) %>%
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>% 
	group_by(iteration) %>% 
	top_n(n = 1, wt = growth_rate) %>% 
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(q2.5=quantile(temperature, probs=0.025),
			  q97.5=quantile(temperature, probs=0.975),
			  mean = mean(temperature)) %>% 
	mutate(term = "Topt")



predictions_sub_dtr9_bite_rate <- predictions_dtr9_bite_rate %>% 
	mutate(iteration = rownames(.)) %>% 
	dplyr::select(iteration, everything()) %>% 
	sample_n(size = 1000)





predictions_long_dtr9_bite_rate <- predictions_sub_dtr9_bite_rate %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1))

predictions_summary_dtr9_bite_rate <- predictions_dtr9_bite_rate %>% 
	mutate(iteration = rownames(.)) %>% 
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>% 
	group_by(temperature) %>% 
	summarise(q2.5=quantile(growth_rate, probs=0.025),
			  q97.5=quantile(growth_rate, probs=0.975),
			  mean = mean(growth_rate)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "bite_rate_dtr9")


predictions_long2_dtr9_bite_rate <- predictions_long_dtr9_bite_rate %>% 
	mutate(temperature = as.numeric(temperature))

data_bite_rate_dtr9_2 <- data.specific_bite_rate_dtr9 %>% 
	mutate(treatment = "bite_rate_dtr9")

bite_rate_dtr9_plot <- ggplot() + 
	geom_ribbon(aes(x = temperature, ymax = q97.5, ymin = q2.5), data = predictions_summary_dtr9_bite_rate, fill = "#66c2a5") +
	# geom_line(aes(x = temperature, y = growth_rate, group = iteration), alpha = 0.05, size = 1, data = predictions_long2) +
	geom_point(aes(x = T, y = trait), data = data_bite_rate_dtr9_2, color = "darkgrey", alpha = 0.7, size = 2) + geom_hline(yintercept = 0) +
	geom_line(aes(x = temperature, y = mean), data = predictions_summary_dtr9_bite_rate, color = "black") +
	geom_line(aes(x = temperature, y = q2.5), data = predictions_summary_dtr9_bite_rate, color = "black", linetype = "dashed") +
	geom_line(aes(x = temperature, y = q97.5), data = predictions_summary_dtr9_bite_rate, color = "black", linetype = "dashed") +
	ylab("Bite rate") + 
	# xlab("Temperature (°C)") +
	xlab("") +
	ylim(0, 1) + xlim(0, 45) 


params_plot_dtr9_bite_rate <- b_params_dtr9_bite_rate %>% 
	filter(term %in% c("cf.T0", "cf.Tm")) %>% 
	dplyr::select(term, mean, `2.5%`, `97.5%`) %>% 
	# gather(key = term, value = value, col = c(2:4)) %>% View
	ggplot() + geom_pointrange(aes(x = term, y = mean, ymin = `2.5%`, ymax = `97.5%`)) +
	coord_flip() + ylab("Temperature (°C)") 





library(patchwork)


bite_rate_dtr9_plot / params_plot_dtr9_bite_rate +  plot_layout(heights = c(4, 1))



# Bite rate dir 12 --------------------------------------------------------


model.out_bite_rate_dtr12_briere


bite.rate.dtr12.preds <- model.out_bite_rate_dtr12_briere$BUGSoutput$sims.list$z.trait.mu.pred


library(tidyverse)
b_params_dtr12_bite_rate <- as.data.frame(model.out_bite_rate_dtr12_briere$BUGSoutput$summary[1:5,]) %>% 
	rownames_to_column(var = "term")


predictions_dtr12_bite_rate <- as.data.frame(model.out_bite_rate_dtr12_briere$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)  ### columns are temperatures, rows are iterations
colnames(predictions_dtr12_bite_rate) <- Temp.xs

predictions_sub_dtr12_bite_rate <- predictions_dtr12_bite_rate %>% 
	mutate(iteration = rownames(.)) %>% 
	dplyr::select(iteration, everything()) %>% 
	sample_n(size = 1000)


predictions_long_dtr12_bite_rate <- predictions_sub_dtr12_bite_rate %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1))

predictions_summary_dtr12_bite_rate <- predictions_dtr12_bite_rate %>% 
	mutate(iteration = rownames(.)) %>% 
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>% 
	group_by(temperature) %>% 
	summarise(q2.5=quantile(growth_rate, probs=0.025),
			  q97.5=quantile(growth_rate, probs=0.975),
			  mean = mean(growth_rate)) %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "bite_rate_dtr12")

topt_bite_rate_dtr12 <- predictions_dtr12_bite_rate %>% 
	mutate(iteration = rownames(.)) %>%
	dplyr::select(iteration, everything()) %>% 
	gather(key = temperature, value = growth_rate, 2:(N.Temp.xs + 1)) %>% 
	group_by(iteration) %>% 
	top_n(n = 1, wt = growth_rate) %>% 
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	summarise(q2.5=quantile(temperature, probs=0.025),
			  q97.5=quantile(temperature, probs=0.975),
			  mean = mean(temperature)) %>% 
	mutate(term = "Topt")


predictions_long2_dtr12_bite_rate <- predictions_long_dtr12_bite_rate %>% 
	mutate(temperature = as.numeric(temperature))

data_bite_rate_dtr12_2 <- data_bite_rate_dtr12 %>% 
	mutate(treatment = "bite_rate_dtr12")



bite_rate_dtr12_plot <- ggplot() + 
	geom_ribbon(aes(x = temperature, ymax = q97.5, ymin = q2.5), data = predictions_summary_dtr12_bite_rate, fill = "#8da0cb") +
	# geom_line(aes(x = temperature, y = growth_rate, group = iteration), alpha = 0.05, size = 1, data = predictions_long2) +
	geom_point(aes(x = T, y = trait), data = data_bite_rate_dtr12, color = "#8da0cb", alpha = 0.7, size = 2) + geom_hline(yintercept = 0) +
	geom_line(aes(x = temperature, y = mean), data = predictions_summary_dtr12_bite_rate, color = "black") +
	geom_line(aes(x = temperature, y = q2.5), data = predictions_summary_dtr12_bite_rate, color = "black", linetype = "dashed") +
	geom_line(aes(x = temperature, y = q97.5), data = predictions_summary_dtr12_bite_rate, color = "black", linetype = "dashed") +
	ylab("Bite rate") + 
	# xlab("Temperature (°C)") +
	xlab("") +
	ylim(0, 1) + xlim(0, 45) 


params_plot_dtr12_bite_rate <- b_params_dtr12_bite_rate %>% 
	filter(term %in% c("cf.T0", "cf.Tm")) %>% 
	dplyr::select(term, mean, `2.5%`, `97.5%`) %>% 
	# gather(key = term, value = value, col = c(2:4)) %>% View
	ggplot() + geom_pointrange(aes(x = term, y = mean, ymin = `2.5%`, ymax = `97.5%`)) +
	coord_flip() + ylab("Temperature (°C)") 





library(patchwork)


bite_rate_dtr12_plot / params_plot_dtr12_bite_rate +  plot_layout(heights = c(4, 1))


### now combine constant and fluctuating

bite_rate_dtr12_plot <- ggplot() + 
	geom_ribbon(aes(x = temperature, ymax = q97.5, ymin = q2.5), data = predictions_summary_dtr12_bite_rate, fill = "#66c2a5") +
	geom_point(aes(x = T, y = trait), data = data_bite_rate_dtr12, color = "darkgrey", alpha = 0.7, size = 2) + geom_hline(yintercept = 0) +
	geom_line(aes(x = temperature, y = mean), data = predictions_summary_dtr12_bite_rate, color = "black") +
	geom_line(aes(x = temperature, y = q2.5), data = predictions_summary_dtr12_bite_rate, color = "black", linetype = "dashed") +
	geom_line(aes(x = temperature, y = q97.5), data = predictions_summary_dtr12_bite_rate, color = "black", linetype = "dashed") +
	ylab("Bite rate") + 
	# xlab("Temperature (°C)") +
	xlab("") +
	ylim(0, 1) + xlim(0, 45) 


bite_rate_constant_plot <- ggplot() + 
	geom_ribbon(aes(x = temperature, ymax = q97.5, ymin = q2.5), data = predictions_summary, fill = "#66c2a5") +
	geom_point(aes(x = T, y = trait), data = data, color = "darkgrey", alpha = 0.7, size = 2) + geom_hline(yintercept = 0) +
	geom_line(aes(x = temperature, y = mean), data = predictions_summary, color = "black") +
	geom_line(aes(x = temperature, y = q2.5), data = predictions_summary, color = "black", linetype = "dashed") +
	geom_line(aes(x = temperature, y = q97.5), data = predictions_summary, color = "black", linetype = "dashed") +
	ylab("Bite rate") + 
	# xlab("Temperature (°C)") +
	xlab("") +
	ylim(0, 1) + xlim(0, 45) 


predictions_summary_dtr12_bite_rate2 <- predictions_summary_dtr12_bite_rate 

combined_cons_dtr9_dtr12_bite_rate_plot <- ggplot() + 
	geom_ribbon(aes(x = temperature, ymax = q97.5, ymin = q2.5), data = predictions_summary_dtr12_bite_rate, fill = "#8da0cb") +
	# geom_line(aes(x = temperature, y = growth_rate, group = iteration), alpha = 0.05, size = 1, data = predictions_long2) +
	geom_point(aes(x = T, y = trait), data = data_bite_rate_dtr12, color = "#8da0cb", alpha = 0.7, size = 2) + geom_hline(yintercept = 0) +
	geom_line(aes(x = temperature, y = mean), data = predictions_summary_dtr12_bite_rate, color = "black") +
	geom_line(aes(x = temperature, y = q2.5), data = predictions_summary_dtr12_bite_rate, color = "black", linetype = "dashed") +
	geom_line(aes(x = temperature, y = q97.5), data = predictions_summary_dtr12_bite_rate, color = "black", linetype = "dashed") +
	geom_ribbon(aes(x = temperature, ymax = q97.5, ymin = q2.5), data = predictions_summary, fill = "#66c2a5") +
	geom_point(aes(x = T, y = trait), data = data_constant_bite_rate, color = "#66c2a5", alpha = 0.7, size = 2) + geom_hline(yintercept = 0) +
	geom_line(aes(x = temperature, y = mean), data = predictions_summary, color = "black") +
	geom_line(aes(x = temperature, y = q2.5), data = predictions_summary, color = "black", linetype = "dashed") +
	geom_line(aes(x = temperature, y = q97.5), data = predictions_summary, color = "black", linetype = "dashed") +
	geom_ribbon(aes(x = temperature, ymax = q97.5, ymin = q2.5), data = predictions_summary_dtr9_bite_rate, fill = "#fc8d62") +
	geom_point(aes(x = T, y = trait), data = data.specific_bite_rate_dtr9, color = "#fc8d62", alpha = 0.7, size = 2) + geom_hline(yintercept = 0) +
	geom_line(aes(x = temperature, y = mean), data = predictions_summary_dtr9_bite_rate, color = "black") +
	geom_line(aes(x = temperature, y = q2.5), data = predictions_summary_dtr9_bite_rate, color = "black", linetype = "dashed") +
	geom_line(aes(x = temperature, y = q97.5), data = predictions_summary_dtr9_bite_rate, color = "black", linetype = "dashed") +
	ylab("Bite rate") + 
	# xlab("Temperature (°C)") +
	xlab("") +
	ylim(0, 1) + xlim(0, 45) 


params_dtr9_bite_rate <- b_params_dtr9_bite_rate %>% 
	filter(term %in% c("cf.T0", "cf.Tm")) %>% 
	mutate(term = ifelse(term == "cf.T0", "Tmin", "Tmax")) %>% 
	dplyr::select(term, mean, `2.5%`, `97.5%`) %>% 
	mutate(treatment = "bite_rate_dtr9")
params_dtr12_bite_rate <- b_params_dtr12_bite_rate %>% 
	filter(term %in% c("cf.T0", "cf.Tm")) %>% 
	mutate(term = ifelse(term == "cf.T0", "Tmin", "Tmax")) %>% 
	dplyr::select(term, mean, `2.5%`, `97.5%`) %>% 
	mutate(treatment = "bite_rate_dtr12")

params_constant_bite_rate <- b_params_bite_rate_briere_constant %>% 
	filter(term %in% c("cf.T0", "cf.Tm")) %>% 
	mutate(term = ifelse(term == "cf.T0", "Tmin", "Tmax")) %>% 
	dplyr::select(term, mean, `2.5%`, `97.5%`) %>% 
	mutate(treatment = "bite_rate_constant")

topt_bite_rate_dtr12_2 <- topt_bite_rate_dtr12 %>% 
	mutate(treatment = "bite_rate_dtr12")

topt_bite_rate_dtr9_2 <- topt_dtr9_bite_rate %>% 
	mutate(treatment = "bite_rate_dtr9")

topt_bite_rate_constant2 <- topt_constant_bite_rate %>% 
	mutate(treatment = "bite_rate_constant")

topts <- bind_rows(topt_bite_rate_dtr12_2, topt_bite_rate_dtr9_2,topt_bite_rate_constant2)

params <- bind_rows(params_constant_bite_rate, params_dtr12_bite_rate, params_dtr9_bite_rate) %>% 
	rename(q2.5 = `2.5%`) %>% 
	rename(q97.5 = `97.5%`)

all_params <- bind_rows(topts, params)


# bring all the things to graph together ----------------------------------
library(plotrix)
all_data_bite_rate <- bind_rows(data_bite_rate_dtr12_2, data_bite_rate_dtr9_2, data_constant_bite_rate_2) %>% 
	# mutate(treatment = factor(treatment, levels = c("bite_rate_constant", "bite_rate_dtr9", "bite_rate_dtr12"))) %>% 
	group_by(treatment, `T`) %>% 
	summarise(q2.5=quantile(trait, probs=0.025),
			  q97.5=quantile(trait, probs=0.975),
			  std_error = std.error(trait),
			  mean = mean(trait)) 

View(all_data_bite_rate)
	
all_predictions <- bind_rows(predictions_summary_dtr12_bite_rate2, predictions_summary_dtr9_bite_rate, predictions_summary)

#66c2a5
#fc8d62
#8da0cb

str(all_predictions)

all_predictions2 <- all_predictions 
# %>% 
# 	mutate(treatment = factor(treatment, levels = c("bite_rate_constant", "bite_rate_dtr9", "bite_rate_dtr12")))
	

str(all_predictions2)
levels(all_predictions2$treatment)
levels(all_data_bite_rate$treatment)


### something is wrong here, come back tomorrow (Jan 5)
combined_cons_dtr9_dtr12_bite_rate_plot_mean <- ggplot() + 
	# geom_line(aes(x = temperature, y = mean, color = treatment), data = subset(all_predictions2, treatment == "bite_rate_constant")) +
	# geom_line(aes(x = temperature, y = q2.5), data = all_predictions2, linetype = "dashed") +
	# geom_line(aes(x = temperature, y = q97.5), data = all_predictions2, linetype = "dashed") +
	# geom_line(aes(x = temperature, y = mean), data = all_predictions2) +
	# geom_line(aes(x = temperature, y = q97.5), data = subset(all_predictions2, treatment == "bite_rate_constant"), linetype = "dashed") +
	# geom_line(aes(x = temperature, y = q2.5), data = subset(all_predictions2, treatment == "bite_rate_constant"), linetype = "dashed") +
	geom_line(aes(x = temperature, y = mean, color = treatment), data = all_predictions2) +
	# geom_line(aes(x = temperature, y = q2.5), data = subset(all_predictions2, treatment == "bite_rate_dtr9"), linetype = "dashed") +
	# geom_line(aes(x = temperature, y = q97.5), data = subset(all_predictions2, treatment == "bite_rate_dtr9"), linetype = "dashed") +
	# geom_line(aes(x = temperature, y = mean), data = subset(all_predictions2, treatment == "bite_rate_dtr9")) +
	# geom_line(aes(x = temperature, y = q2.5), data = subset(all_predictions2, treatment == "bite_rate_dtr12"),  linetype = "dashed") +
	# geom_line(aes(x = temperature, y = q97.5), data = subset(all_predictions2, treatment == "bite_rate_dtr12"),  linetype = "dashed") +
	# geom_line(aes(x = temperature, y = mean), data = subset(all_predictions2, treatment == "bite_rate_dtr12")) +
	geom_ribbon(aes(x = temperature, ymax = q97.5, ymin = q2.5, fill = treatment), data = all_predictions2, alpha = 0.5) +
	# geom_ribbon(aes(x = temperature, ymax = q97.5, ymin = q2.5, fill = treatment), data = subset(all_predictions2, treatment == "bite_rate_dtr9"), alpha = 0.5) +
	geom_pointrange(aes(x = `T`, y = mean, ymin = mean - std_error, ymax = mean + std_error, color = treatment), data = all_data_bite_rate, alpha = 0.7) + geom_hline(yintercept = 0) + 
	ylab("Bite rate") + 
	xlab("") +
	ylim(0, 0.55) + xlim(0, 45) 
ggsave("figures/all_bite_rate_facet.pdf", width = 7, height = 4)

levels(all_data_bite_rate$treatment)

	params_plot_bite_rate <- ggplot() +
	geom_pointrange(aes(x = term, y = mean, ymin = `2.5%`, ymax = `97.5%`), data = params_constant_bite_rate, color = "#66c2a5") +
	geom_pointrange(aes(x = term, y = mean, ymin = `2.5%`, ymax = `97.5%`), data = params_dtr9_bite_rate, color = "#fc8d62") +
	geom_pointrange(aes(x = term, y = mean, ymin = `2.5%`, ymax = `97.5%`), data = params_dtr12_bite_rate, color = "#8da0cb") +
	geom_pointrange(aes(x = term, y = mean, ymin = q2.5, ymax = q97.5), data = topt_bite_rate_dtr12, color = "#8da0cb") +
	geom_pointrange(aes(x = term, y = mean, ymin = q2.5, ymax = q97.5), data = topt_bite_rate_dtr9_2, color = "#fc8d62") +
	geom_pointrange(aes(x = term, y = mean, ymin = q2.5, ymax = q97.5), data = topt_bite_rate_constant2, color = "#66c2a5") +
	coord_flip() + ylab("Temperature (°C)") 
	
	params_plot_bite_rate <- ggplot() +
		geom_pointrange(aes(x = term, y = mean, ymin = q2.5, ymax = q97.5, color = treatment), data = all_params, position=position_dodge(width=1)) +
		coord_flip() + ylab("Temperature (°C)") 
	
	View(all_params)
	
	
	bite_rate_plot <- combined_cons_dtr9_dtr12_bite_rate_plot_mean / params_plot_bite_rate + plot_layout(heights = c(3, 0.5))
	
ggsave("figures/bite_rate_plot.pdf", bite_rate_plot, width = 7, height = 5)


### now onto the other metrics of performance!



# New plotting centralized ------------------------------------------------


### bc data
"data-processed/predictions_bc_constant_summary.csv"
"data-processed/params_bc_constant_all.csv"
"data-processed/data_bc_constant_sum.csv"

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
library(cowplot)
theme_set(theme_cowplot())

View(bite.rate.dtr12.predictions)

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

### Ok now move onto the other parameters



# actual lifespan ---------------------------------------------------------


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


View(lifespan.constant.sum )

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


### ok now merge all three plots together! woohoo!
### ok something is wrong with the lifespan plot -- we are missing the constant data and the fits look wrong

all_plots <- wrap_plots(eggs_plot_all, lifespan_plot_all, bite_rate_plot_all)
ggsave('figures/all_plots_combined-feb4-3.pdf', all_plots, width = 18, height = 6)
### OK so the best fits are briere for bite rate and egg production but quadratic for lifespan.
