###### Load packages
library(tidyverse)
library(R2jags)
library(coda)
library(cowplot)
library(mcmcplots)

######  Load raw trait data
cdata <- read_csv("data-raw/constant.individual.trait.csv")
f9data <- read_csv("data-raw/fluc9.individual.trait.csv")
f12data <- read_csv("data-raw/fluc12.individual.trait.csv")

###### Summarizing + plotting data - method #1 (lifespan only)

# Sumarize dataframes
cdata.sum <- cdata %>%
	group_by(Treatment) %>%
	summarise(DTR = "constant", lifespan.mean = mean(lifespan), lifespan.SE = sd(lifespan)/sqrt(n()),
			  bite.rate.mean = mean(lifespan), bite.rate.SE = sd(bite.rate)/sqrt(n()),
			  lifetime.eggs.mean = mean(lifetime.eggs), lifetime.eggs.SE = sd(lifetime.eggs)/sqrt(n()))

f9data.sum <- f9data %>%
	group_by(Treatment) %>%
	summarise(DTR = "+/- 9", lifespan.mean = mean(lifespan), lifespan.SE = sd(lifespan)/sqrt(n()),
			  bite.rate.mean = mean(lifespan), bite.rate.SE = sd(bite.rate)/sqrt(n()),
			  lifetime.eggs.mean = mean(lifetime.eggs), lifetime.eggs.SE = sd(lifetime.eggs)/sqrt(n()))

f12data.sum <- f12data %>%
	group_by(Treatment) %>%
	summarise(DTR = "+/- 12", lifespan.mean = mean(lifespan), lifespan.SE = sd(lifespan)/sqrt(n()),
			  bite.rate.mean = mean(lifespan), bite.rate.SE = sd(bite.rate)/sqrt(n()),
			  lifetime.eggs.mean = mean(lifetime.eggs), lifetime.eggs.SE = sd(lifetime.eggs)/sqrt(n()))

# Plot the separate dataframes with manual jittering
ggplot() +
	geom_point(data = cdata.sum, aes(x = Treatment, y = lifespan.mean)) + ylim(0,55) +
	geom_errorbar(data = cdata.sum, aes(x = Treatment, ymin = lifespan.mean - lifespan.SE, ymax = lifespan.mean + lifespan.SE), width = 0) +
	geom_point(data = f9data.sum, aes(x = Treatment + 0.2, y = lifespan.mean), color = "blue") +
	geom_errorbar(data = f9data.sum, aes(x = Treatment + 0.2, ymin = lifespan.mean - lifespan.SE, ymax = lifespan.mean + lifespan.SE), width = 0,
				  color = "blue") +
	geom_point(data = f12data.sum, aes(x = Treatment + 0.4, y = lifespan.mean), color = "dodgerblue") +
	geom_errorbar(data = f12data.sum, aes(x = Treatment + 0.4, ymin = lifespan.mean - lifespan.SE, ymax = lifespan.mean + lifespan.SE), width = 0,
				  color = "dodgerblue")

###### Plotting data - method #2 (aka, Marta slowly learning dplyr/ggplot stuff) (lifespan only)

# Combine into a single dataframe, add jittering, plot
alldata.sum <- rbind(cdata.sum, f9data.sum, f12data.sum)
alldata.sum$jitter <- c(rep(0,6), rep(0.2,5), rep(0.4,5))
ggplot() +
	geom_point(data = alldata.sum, aes(x = Treatment + jitter, y = lifespan.mean, color = DTR)) +
	geom_errorbar(data = alldata.sum, aes(x = Treatment + jitter, ymin = lifespan.mean - lifespan.SE, ymax = lifespan.mean + lifespan.SE, color = DTR), width = 0) +
	scale_color_manual(values = c("dodgerblue", "blue", "black"))

###### Summarizing + plotting data - method #3 (aka, Marta slowly learning dplyr/ggplot stuff) (all 3 traits: lifespan, biting rate, fecundity)

# Add in DTR Treatment columns, combine in a single dataframe
cdata$DTR <- "constant"
f9data$DTR <- "+/- 9"
f12data$DTR <- "+/- 12"
alldata <- rbind(cdata, f9data, f12data)

alldata <- alldata %>% mutate(TempK = Treatment + 273.15)

# Summarize dataframe, add jitter
alldata.sum.2 <- alldata %>%
	group_by(DTR, Treatment) %>%
	summarise(lifespan.mean = mean(lifespan), lifespan.SE = sd(lifespan)/sqrt(n()),
			  bite.rate.mean = mean(bite.rate), bite.rate.SE = sd(bite.rate)/sqrt(n()),
			  lifetime.eggs.mean = mean(lifetime.eggs), lifetime.eggs.SE = sd(lifetime.eggs)/sqrt(n()))
alldata.sum.2$jitter <- c(rep(0.4,5), rep(0.2,5), rep(0,6))

# Plot lifespan
p.lf <- ggplot() + 
	geom_point(data = alldata.sum.2, aes(x = Treatment + jitter, y = lifespan.mean, color = DTR)) +
	geom_errorbar(data = alldata.sum.2, aes(x = Treatment + jitter, ymin = lifespan.mean - lifespan.SE, ymax = lifespan.mean + lifespan.SE, color = DTR), width = 0) +
	scale_color_manual(values = c("dodgerblue", "blue", "black")) + theme(legend.position = c(0.6,0.8)) +
	ylab("Lifespan (days)") + xlab("Temperature (°C)") 

# Plot biting rate
p.br <- ggplot() + 
	geom_point(data = alldata.sum.2, aes(x = Treatment + jitter, y = bite.rate.mean, color = DTR)) +
	geom_errorbar(data = alldata.sum.2, aes(x = Treatment + jitter, ymin = bite.rate.mean - bite.rate.SE, ymax = bite.rate.mean + bite.rate.SE, color = DTR), width = 0) +
	scale_color_manual(values = c("dodgerblue", "blue", "black")) + theme(legend.position = "none") +
	ylab(expression(paste("Biting rate (days"^-1,")"))) + xlab("Temperature (°C)") 

# Plot fecundity
p.f <- ggplot() + 
	geom_point(data = alldata.sum.2, aes(x = Treatment + jitter, y = lifetime.eggs.mean, color = DTR)) +
	geom_errorbar(data = alldata.sum.2, aes(x = Treatment + jitter, ymin = lifetime.eggs.mean - lifetime.eggs.SE, ymax = lifetime.eggs.mean + lifetime.eggs.SE, color = DTR), width = 0) +
	scale_color_manual(values = c("dodgerblue", "blue", "black")) + theme(legend.position = "none") +
	ylab("Lifetime fecundity (eggs)") + xlab("Temperature (°C)") 

plot_grid(p.lf, p.br, p.f, nrow = 1)
ggsave("figures/SummarizedTraitData.pdf")
