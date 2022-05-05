###simplistic code to make the evaluation of RS run faster
#Need faster run time to (1)replace values for complete Bayesian model output instead of the mean summary outputs
#                            (1-rationale)- Need to generated error around the rate summation estimates, I feel more comfortable summarizing over the full model evaluation
#Need faster run time to (2)use with raster objects when generating a Parton-Logan model from a raster cell (Tave, Tmax, Tmin [Tdtr]) and applying rate summation using that temperature profile

#Dataframes
#1. Df of temperature profiles (aka, a LP model with a dtr of 9 or 12 for each mean temp of seq(0,50,0.1) for each hr of the day)
#     1-cols = temp program for that mean temp / dtr combo
#     1-rows = temperature at that hour of the day (24hr period)
#Final Dimensions = [24, 1003]

br.c.out <- read.csv("br.c.out.ME.csv")

#1. Simplified Dims[24,3]
temp.programs.example <- data.frame(dtr.9.mean.20 = c(17.0,16.7,16.5,16.3,16.2,16.1,16.0,17.9,19.7,21.3,22.7,23.8,24.6,25.0,25.0,24.6,23.8,22.7,21.3,20.1,19.1,18.4,17.8,17.3),
                                    dtr.12.mean.32.4 = c(28.4,28.0,27.7,27.5,27.3,27.2,27.1,29.6,32.0,34.1,36.0,37.5,38.5,39.0,39.0,38.5,37.5,36.0,34.1,32.5,31.2,30.2,29.5,28.8),
                                    dtr.9.mean.28.5 = c(25.5,25.2,25.0,24.8,24.7,24.6,24.5,26.4,28.2,29.8,31.2,32.3,33.1,33.5,33.5,33.1,32.3,31.2,29.8,28.6,27.6,26.9,26.3,25.8))

#Visual of temp programs
plot(temp.programs.example$dtr.9.mean.20, col ="black", pch = 16, ylim = c(16,40), ylab = "Temprature (C)", xlab = "Time of day (hour)")
points(temp.programs.example$dtr.12.mean.32.4, col = "darkcyan", pch = 16)
points(temp.programs.example$dtr.9.mean.28.5, col = "orchid", pch =16)

#2. Df of Bayesian model mean outputs over a temperature gradient
#       2-cols = [1]; temp, [2] mean.out
#       2-rows = temperatures

#1. Simplified Dims[501,1]
temp.gradient.example <- seq(0,50,0.1)
bite.rate.out.example <- br.c.out$median
#load csv
#br.c.out <- read.csv("br.c.out.ME.csv")
model.out.example <- data.frame( temp = temp.gradient.example, br = bite.rate.out.example)

#Visual of trait model
plot(model.out.example$temp, model.out.example$br, col = "seagreen3", ylab = "Trait value", xlab = "Temperature (C)", main= "Bayesian output")


##################################################
# A. Current workflow that works over the Bayesian model summary outputs, but is time consuming
#################################################
hourly.estimates = function(x, temp.performance){
  x <- as.numeric(round(x,4)) #to prevent any structural discrepancies
  value.df <- temp.performance[which(apply(as.matrix(temp.performance$temp), MARGIN = c(1), all.equal, x)=="TRUE"),] #use all.equal to circumvent matching error with ===, extracts the row in temp.performance which temp matches temp.program value
  value <- value.df$median[1] #Median is the second column in the dataframe
  return(value)
}

hourly.estimates2 = function(x, temp.performance){
  x <- as.numeric(x)
  median.rate <- temp.performance$median[match(x, temp.performance$temp)]
  return(median.rate)
}

br.estimate2 <- data.frame(apply(temp.programs.example, MARGIN = 2, hourly.estimates2, temp.performance = br.c.out))

#Uses the function on each trait of interest
br.estimate <- data.frame(apply(temp.programs.example, MARGIN = c(1,2), hourly.estimates, temp.performance = br.c.out))

#compare
library(microbenchmark)
time.test <- microbenchmark::microbenchmark(apply(temp.programs.example, MARGIN = 2, hourly.estimates2, temp.performance = br.c.out), apply(temp.programs.example, MARGIN = c(1,2), hourly.estimates, temp.performance = br.c.out), 
                                            times = 10, check = 'equivalent')

time.test
#Plot results
plot(br.estimate$dtr.9.mean.20, col = "black", ylab = "Trait value", xlab = "Time of day (hour)")

###############################################################################
#####Note I have all the temperature program columns because my subsequent code follows as such
#1) Integrate trait values over a program to generate trait performance with RS
br.estimate.integrated <- colMeans(br.estimate)
#2) Plot these ouputs for a dtr of 9 and 12 across the full temperature gradient (x-axis = Temp (0-50C), y-axis = Estimated trait performance)


