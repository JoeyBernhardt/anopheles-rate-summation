################################# 1. Set-up

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

###### Add in DTR Treatment columns, combine in a single dataframe
cdata$DTR <- "constant"
f9data$DTR <- "+/- 9"
f12data$DTR <- "+/- 12"
alldata <- rbind(cdata, f9data, f12data)
flucdata <- rbind(f9data, f12data)

######  Load treatment means
meandata <- data.frame(read.csv("data-raw/fluc.program.trait.means.csv", colClasses = c("DTR" ="character")))


################################# 2. Plot trait distributions using individual-level data

################# By treatment type (constant, +/- 9, +/- 12)

##### constant temperature
cdata %>%
	filter(Treatment != 36) %>%
	ggplot(aes(x = lifespan, fill = as.factor(Treatment))) +
	geom_histogram(binwidth = 3, position = "identity") +
	facet_grid(Treatment ~ .)

cdata %>%
	filter(Treatment != 36) %>%
	ggplot(aes(x = bite.rate, fill = as.factor(Treatment))) +
	geom_histogram(binwidth = 0.03, position = "identity") +
	facet_grid(Treatment ~ .)

cdata %>%
	filter(Treatment != 36) %>%
	ggplot(aes(x = lifetime.eggs, fill = as.factor(Treatment))) +
	geom_histogram(binwidth = 50, position = "identity") +
	facet_grid(Treatment ~ .)

##### fluctuating temperature - 9
ggplot(f9data, aes(x = lifespan, fill = as.factor(Treatment))) +
	geom_histogram(binwidth = 3, position = "identity") +
	facet_grid(Treatment ~ .)

ggplot(f9data, aes(x = bite.rate, fill = as.factor(Treatment))) +
	geom_histogram(binwidth = 0.03, position = "identity") +
	facet_grid(Treatment ~ .)

ggplot(f9data, aes(x = lifetime.eggs, fill = as.factor(Treatment))) +
	geom_histogram(binwidth = 50, position = "identity") +
	facet_grid(Treatment ~ .)

##### fluctuating temperature - 12
ggplot(f12data, aes(x = lifespan, fill = as.factor(Treatment))) +
	geom_histogram(binwidth = 3, position = "identity") +
	facet_grid(Treatment ~ .)

ggplot(f12data, aes(x = bite.rate, fill = as.factor(Treatment))) +
	geom_histogram(binwidth = 0.03, position = "identity") +
	facet_grid(Treatment ~ .)

ggplot(f12data, aes(x = lifetime.eggs, fill = as.factor(Treatment))) +
	geom_histogram(binwidth = 50, position = "identity") +
	facet_grid(Treatment ~ .)


################################# 3. ANOVA and t-tests

# Subset data by temperature
data.32 <- subset(alldata, Treatment == 32)
data.28 <- subset(alldata, Treatment == 28)
data.24 <- subset(alldata, Treatment == 24)
data.20 <- subset(alldata, Treatment == 20)
data.16 <- subset(alldata, Treatment == 16)

fdata.32 <- subset(flucdata, Treatment == 32)
fdata.28 <- subset(flucdata, Treatment == 28)
fdata.24 <- subset(flucdata, Treatment == 24)
fdata.20 <- subset(flucdata, Treatment == 20)
fdata.16 <- subset(flucdata, Treatment == 16)


##### Biting rate
aov.br.16 <- aov(bite.rate ~ DTR, data = data.16)
summary(aov.br.16)
t.br.16 <- t.test(bite.rate ~ DTR, data = fdata.16)
t.br.16

aov.br.20 <- aov(bite.rate ~ DTR, data = data.20)
summary(aov.br.20)
t.br.20 <- t.test(bite.rate ~ DTR, data = fdata.20)
t.br.20

aov.br.24 <- aov(bite.rate ~ DTR, data = data.24)
summary(aov.br.24)
t.br.24 <- t.test(bite.rate ~ DTR, data = fdata.24)
t.br.24

aov.br.28 <- aov(bite.rate ~ DTR, data = data.28)
summary(aov.br.28)
t.br.28 <- t.test(bite.rate ~ DTR, data = fdata.28)
t.br.28

aov.br.32 <- aov(bite.rate ~ DTR, data = data.32)
summary(aov.br.32)
t.br.32 <- t.test(bite.rate ~ DTR, data = fdata.32)
t.br.32


##### Lifespan
aov.lf.16 <- aov(lifespan ~ DTR, data = data.16)
summary(aov.lf.16)
t.lf.16 <- t.test(lifespan ~ DTR, data = fdata.16)
t.lf.16

aov.lf.20 <- aov(lifespan ~ DTR, data = data.20)
summary(aov.lf.20)
t.lf.20 <- t.test(lifespan ~ DTR, data = fdata.20)
t.lf.20

aov.lf.24 <- aov(lifespan ~ DTR, data = data.24)
summary(aov.lf.24)
t.lf.24 <- t.test(lifespan ~ DTR, data = fdata.24)
t.lf.24

aov.lf.28 <- aov(lifespan ~ DTR, data = data.28)
summary(aov.lf.28)
t.lf.28 <- t.test(lifespan ~ DTR, data = fdata.28)
t.lf.28

aov.lf.32 <- aov(lifespan ~ DTR, data = data.32)
summary(aov.lf.32)
t.lf.32 <- t.test(lifespan ~ DTR, data = fdata.32)
t.lf.32


##### Lifetime eggs
aov.e.16 <- aov(lifetime.eggs ~ DTR, data = data.16)
summary(aov.e.16)
t.e.16 <- t.test(lifetime.eggs ~ DTR, data = fdata.16)
t.e.16

aov.e.20 <- aov(lifetime.eggs ~ DTR, data = data.20)
summary(aov.e.20)
t.e.20 <- t.test(lifetime.eggs ~ DTR, data = fdata.20)
t.e.20

aov.e.24 <- aov(lifetime.eggs ~ DTR, data = data.24)
summary(aov.e.24)
t.e.24 <- t.test(lifetime.eggs ~ DTR, data = fdata.24)
t.e.24

aov.e.28 <- aov(lifetime.eggs ~ DTR, data = data.28)
summary(aov.e.28)
t.e.28 <- t.test(lifetime.eggs ~ DTR, data = fdata.28)
t.e.28

aov.e.32 <- aov(lifetime.eggs ~ DTR, data = data.32)
summary(aov.e.32)
t.e.32 <- t.test(lifetime.eggs ~ DTR, data = fdata.32)
t.e.32


################################# 4. Plotting treatment means

ggplot(meandata, aes(Treatment, bite.rate, color = DTR)) + 
	geom_point() + scale_color_manual(values=c("#0066FF", "skyblue", "black"))

ggplot(meandata, aes(Treatment, lifespan, color = DTR)) + 
	geom_point() + scale_color_manual(values=c("#0066FF", "skyblue", "black"))

ggplot(meandata, aes(Treatment, lifetime.eggs, color = DTR)) + 
	geom_point() + scale_color_manual(values=c("#0066FF", "skyblue", "black"))

