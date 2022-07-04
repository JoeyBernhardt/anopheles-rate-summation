
## Kerri Miazgowicz, University of Georgia
## October 16, 2018; modified May 22,2019 to code for a continous Parton-Logan model
##
## Purpose: Generate rate summation estimates based on thermal performance at constant temperature (Bayesian fits) 
##
## Contents:
# (1) Load data and packages
# (2) Evaluate trait median curve values at 0.1C (what the temp dataframe is (Note: I can do at 0.01C intervals but I'll have to refit the Bayesian models))
# (4)a  Create a function which rounds the the temperature program df to the nearest 0.1C,
# (4)b Creates a new dataframe which lists trait performance at each time interval
# (4)c Computes the average of the hourly estimates to generate a single value
# (4)d Place these values into a dataframe which consists of 'temp' & 'traitvalue'
# (5) A provisional way to geneate a TPC for the fluctuaiton estimate (Bayesian fit); once I can code the Logan-Parton model, I can incorporate that to generate estimates with fluctuations around all temperatures at 0.1C intervals

##########
###### 1. Load data and packages
##########

# Set working directory
mainDir = "C:/Users/Kerri/Desktop/Fluctuation_BayesianFits"
setwd(mainDir)

# Check whether there's a folder in the directory for saving plots
# If not, create one
subDir = "plots"
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)

# Load libraties for plotting traits
library(plyr) # Slices and dices data
library(plotrix) # For standard error function
library(coda)
library(dplyr)

############ Load traits fits
# Fits from Miazgowicz constant dataset
load("saved posteriors/constant_lifespan_quadratic_uniform.Rdata")
lifespan.constant.preds <- model.out$BUGSoutput$sims.list$z.trait.mu.pred
load("saved posteriors/constant_lifetime.eggs_quadratic_uniform.Rdata")
lifetime.eggs.constant.preds <-model.out$BUGSoutput$sims.list$z.trait.mu.pred
load("saved posteriors/constant_bite.rate_briere_uniform.Rdata")
bite.rate.constant.preds <-model.out$BUGSoutput$sims.list$z.trait.mu.pred

# Temperature sequences that was used to make all calculations
load("saved posteriors/temps.Rdata")
N.Temp.xs <-length(Temp.xs)


##The Parton-Logan model -> formula taken from the Murdock Lab Fluctuation Program Calculator Excel Sheet (Creator: Krijn Paaijmans)
##A sinusoial  relationship during the day followed by an exponential decay throughout the course of the day

##Part 1: Generate a function to calculate the needed parameters based on sunrise, sunset, L:D hours
#Note: to simplify the coding I will use military time 0/24 = midnight; 1 = 1am; 14 = 2pm ect.
night.length = 12 #[$0$3]
day.length = 24 - night.length
time.sunrise = 6
time.sunset = time.sunrise + day.length  #[$O$5]
nocturnal.constant.tau = 4 #[$O$7]
time.Tmax.afterNoon = 1.5

#Make a loop that creates an array with columns containing the program values for row hour values (1hr increments)***
###########
######################
#initializing values of interest to create Parton-Logan models for;
dtr = c(9,12)
mean.temp = seq(0,45,1) #will be a vector of integers from 0 to 45C when implemented (n = 46)
temp.programs.2 <- data.frame("hr" = seq(0,23,1))

####Loop through dtrs of interest
for ( a in 1:length(dtr)){
####Loop through mean temps of interest
  for(b in 1:length(mean.temp)){

############################

correction.factor = -0.0575824  #only applies for the day.length of 12

median.temp = mean.temp[b] - correction.factor*dtr[a]

#The below variables will need to be calculated for each 'dtr' and 'mean.temp' used [contained either in an arrary or a dataframe]
amplitude = dtr[a]/2
Tmin = median.temp - amplitude
Tmax = median.temp + amplitude
Tsunset = Tmin+(Tmax-Tmin)*sin(pi*day.length/(day.length+2*time.Tmax.afterNoon)) 

#Calculate hourly temperatures  --> Note: to follow suit with the code provided in the excel sheet you start at 6am and the hr goes to 29.5/30; the true hour values are reassigned after the loop implementation
#Step 1: Vector of hours
arbitary.hr <- seq(6,29,1) #starts at 6am for if statment to work properly.... 
#create an empty vector for hr.temp to be populated into of the same length of arbitary.hr

temp <- rep(NA, length(arbitary.hr))

#Step 2: Make a vector of temperatures corresponding to each [] of hr vector [A18] uusing the Parton-Logan model
for(i in 1:length(arbitary.hr)){
  if(arbitary.hr[i]<time.sunset) {hr.temp = Tmin +(Tmax-Tmin)*sin(pi*(arbitary.hr[i]-12+ day.length/2)/(day.length+2*time.Tmax.afterNoon))}  #day sin
  else {hr.temp = (Tmin-Tsunset*exp(-night.length/nocturnal.constant.tau)+(Tsunset-Tmin)*exp(-(arbitary.hr[i]-time.sunset)/nocturnal.constant.tau))/(1-exp(-night.length/nocturnal.constant.tau))}  # night exponetial decay
 temp[i] = hr.temp
   }

#Step 3: Change arbitary.hr to actual.hr 'hr' 
#Combine vectors together into a df
temp.program <- as.data.frame(cbind(arbitary.hr, temp))
#add in actual hour by looping through each element in the vector (loop?)
#create an empty vector for hr to be populated into of the same length of arbitary.hr
hr <- rep(NA, length(arbitary.hr))

for(j in 1: nrow(temp.program)){
 ifelse(temp.program$arbitary.hr[j] > 23, actual.hr <- (temp.program$arbitary.hr[j]-24), actual.hr <- temp.program$arbitary.hr[j])
 hr[j] = actual.hr
}

#create a new column with the correct vector values
temp.program$hr <- hr
#Reorder to start at time = 0/24 (midnight)
temp.program <- temp.program[order(temp.program$hr),] #note row name label are not changed using this method..but the index is fixed
hr <- temp.program$hr
#rename 'temp' column to have a meaningful name for downstream applications prior to appending it to the dataframe
new.name = paste0("X", mean.temp[b], ".dtr", dtr[a])
colnames(temp.program)[colnames(temp.program)=="temp"] <- new.name

#add 'temp' values to master temp.programs.2 df (cbind)
temp.programs.2 <- cbind(temp.programs.2, temp.program[2]) #To circumvent having to directly reference the new col name I use the colm index reference 

} #end of b loop (associated with mean.temp)
   } #end of a loop (associated with dtr)

# Round to nearest 0.1C
temp.programs.2 = round(temp.programs.2,1)

write.csv(temp.programs.2, "data/temp.programs.continuous.csv", row.names = FALSE) #list of all temp programs
###################
########################
###############


##########
###### 2. Evaluate trait median curve values at 0.1C
##########

###### calcPostQuants - Function to calculate quantiles of derived quantity posteriors (for plotting)
# Arguments:
#   input       data frame with posterior distributions of trait predicted over a temperature gradient (i.e., 'trait.preds')
#   temp.list   the temperature gradient itself (i.e., 'Temp.xs')

calcPostQuants = function(input, temp.list) {
  
  # Get length of gradient
  N.grad.xs <- length(temp.list)
  
  # Create output dataframe
  output.df <- data.frame("mean" = numeric(N.grad.xs), "median" = numeric(N.grad.xs), "lowerCI" = numeric(N.grad.xs), "upperCI" = numeric(N.grad.xs))
  
  # Calculate mean & quantiles
  for(i in 1:N.grad.xs){
    output.df$mean[i] <- mean(input[ ,i])
    output.df$lowerCI[i] <- quantile(input[ ,i], 0.025)
    output.df$upperCI[i] <- quantile(input[ ,i], 0.975)
    output.df$median[i] <- quantile(input[ ,i], 0.5)
  }
  
  output.df # return output
  
}

###### Calculate trait posteriors and quantiles
bite.constant.out = calcPostQuants(bite.rate.constant.preds, Temp.xs)
lifetime.eggs.constant.out = calcPostQuants(lifetime.eggs.constant.preds, Temp.xs)
lifespan.constant.out = calcPostQuants(lifespan.constant.preds, Temp.xs)

###### Modify outputs into a format useable for performing rate summation
bite.constant.out$temp = Temp.xs
lifetime.eggs.constant.out$temp = Temp.xs
lifespan.constant.out$temp = Temp.xs


##########
###### 4 (a) Round the the temperature program df to the nearest 0.1C, 
######   (b) Creates a new dataframe which lists trait performance at each time interval ('median' from trait.out)
######   (c) Computes the average of the hourly estimates to generate a single value
######   (d) Place these values into a dataframe which consists of 'temp' & 'traitvalue'
##########
# Round to nearest 0.1C
#temp.programs = round(temp.programs,1)

#Place into a list for apply statements later (Note; use hr for corresponding hr of temp.program)
#temp.programs.list <-  list( temp_program = temp.programs[2:length(temp.programs)])  #removes hr, creates a 1 unit list of the dataframe containing all of the hourly temperatures for both deters.

#A function which returns the value of the hourly estimate of trait performance based on the values from temperature.program 
# Arguments:
#   temp.performance       a data frame with posterior distributions of trait predicted over a temperature gradient (i.e., 'trait.preds')
hourly.estimates = function(x, temp.performance){
  x <- as.numeric(x) #to prevent any structural discrepancies
  value.df <- temp.performance[which(apply(as.matrix(temp.performance$temp), MARGIN = c(1), all.equal, x)=="TRUE"),] #use all.equal to circumvent matching error with ===, extracts the row in temp.performance which temp matches temp.program value
  value <- value.df$median[1] #Median is the second column in the dataframe
  return(value)
  }

#Uses the function on each trait of interest
bite.estimate <- data.frame(apply(temp.programs.2, MARGIN = c(1,2), hourly.estimates, temp.performance = bite.constant.out))
lifespan.estimate <-data.frame(apply(temp.programs.2, MARGIN = c(1,2), hourly.estimates, temp.performance =lifespan.constant.out))
lifetime.eggs.estimate <- data.frame(apply(temp.programs.2, MARGIN = c(1,2), hourly.estimates, temp.performance = lifetime.eggs.constant.out))

##Filling in NAs with 0 (aka negative temperatues have no performance values)
bite.estimate[is.na(bite.estimate)] <- 0
lifespan.estimate[is.na(lifespan.estimate)] <- 0
lifetime.eggs.estimate[is.na(lifetime.eggs.estimate)] <- 0

##Fix hour column
bite.estimate$hr <- hr # hr here should still correlate to the reordered hr value from the temp.programs.2
lifespan.estimate$hr <- hr
lifetime.eggs.estimate$hr <- hr

####################
#############
########Visualizing the outputs generated to double-check everything worked as expected

#put all estimates into a list for use of the function below
component.estimation <- list(bite.estimate, lifespan.estimate, lifetime.eggs.estimate)
#test <- component.estimation[[1]] #test$X = hr

plots_hourly <- list() #initializing an object to store all the plots

#Use the non-integrated dataframes
#Have point symbol and linetype correspond to progam -> extracted from the col name
#Note: I code to have a different linetype and shape for dtr 9 or 12; if more dtrs are used this will need to be changed
for(c in 1:length(component.estimation)) {  #outter loop start
   temporary.df <- component.estimation[[c]] #extracts out each df in the list
   for(e in 2:(ncol(temporary.df))) { #inner loop start
    ydata <- colnames(temporary.df)[e]
    dtr.point.symbol = sub("^[^.]*", "",colnames(temporary.df)[e])
    dtr.line.type = sub("^[^.]*", "",colnames(temporary.df)[e])
    if (e == 2)
      { z <- ggplot(temporary.df, aes_string(hr, as.name(paste(ydata)))) +
          geom_point(colour = e, shape = ifelse(dtr.point.symbol == ".dtr9", 1, 4))+
          geom_line(colour = e, linetype = ifelse(dtr.line.type == ".dtr9", 1,2)) +
          theme_classic()
        
      } 
    if(e > 2)
      { z <- z + 
        geom_point(data = temporary.df, aes_string(hr, as.name(paste(ydata))),  shape = ifelse(dtr.point.symbol == ".dtr9", 1, 4), colour = e)+
        geom_line(data = temporary.df, aes_string(hr, as.name(paste(ydata))), linetype = ifelse(dtr.line.type == ".dtr9", 1,2), colour = e)
    }
   plots_hourly[[c]] <- z
   } #inner loop close
} #outter loop close
#Look at hourly estimate plots for
plots_hourly[[1]] #bite
plots_hourly[[2]] #lifespan
plots_hourly[[3]] #lifetime.eggs

###################
######################
########################

#To avoid further issues with data conversion; I will export these as csv files and then reload them in R.
write.csv(bite.estimate, "data/bite.estimate.csv", row.names = FALSE)
write.csv(lifespan.estimate, "data/lifespan.estimate.csv", row.names = FALSE)
write.csv(lifetime.eggs.estimate, "data/lifetime.eggs.estimate.csv", row.names = FALSE)

#Reimport
bite.estimate <- read.csv("data/bite.estimate.csv")
lifespan.estimate <- read.csv("data/lifespan.estimate.csv")
lifetime.eggs.estimate <- read.csv("data/lifetime.eggs.estimate.csv")


#Remove hour column
bite.estimate = bite.estimate %>% dplyr::select(-c( hr))
lifespan.estimate =  lifespan.estimate %>% dplyr::select(-c( hr))
lifetime.eggs.estimate=  lifetime.eggs.estimate %>% dplyr::select(-c( hr))

#Modify this to be correct now that I have continuous data
traits.integrated2 <- data.frame (mean.temp = as.numeric(substr(colnames(bite.estimate),2,3)),
                                 dtr = sub("^[^.]*", "",colnames(bite.estimate)),
                                 program = colnames(bite.estimate),
                                 bite.rate = as.vector(colMeans(bite.estimate, na.rm =TRUE)),
                                 lifespan =  as.vector(colMeans(lifespan.estimate, na.rm = TRUE)),
                                 lifetime.eggs = as.vector(colMeans(lifetime.eggs.estimate, na.rm = TRUE)))
#fix dtr column to remove '.dtr' and change to a numeric data type
traits.integrated2$dtr <- substr(traits.integrated2$dtr, 5,6)
traits.integrated2$dtr <- as.numeric(traits.integrated2$dtr)

#Export final dataframe for use in second bayesian fits (Note: in alternative approach where I estimate around each mean temperature (in 0.1C intervals, I will avoid the need to re-fit the estimations with a Bayesian appraoch, and can instead use the rate summation on the origianl upper and lower CI))
write.csv(traits.integrated2, "data/traits.integrated2.csv", row.names = FALSE)

###Note: If I wanted I could use the above code to employ rate summation on the other components of R0 (MDR, pea, PDR, bc, ect.)

###########Create a graph for traits.integrated2 to show rate summation estimates for each trait--> will use later when I overlay the bayesian fits of observed constant, dtr9 , and dtr12 data
#bite rate estimates
ggplot(traits.integrated2, aes(x= mean.temp, y = bite.rate, colour = factor(dtr), linetype = factor(dtr)))+
  scale_linetype_manual(values = c("longdash", "dotted"))+
  scale_color_manual(values = c("dodgerblue", "purple"))+
  #geom_point()+
  geom_line(size = 1.1)+
  theme_classic()

#lifespan estimates
ggplot(traits.integrated2, aes(x= mean.temp, y = lifespan, colour = factor(dtr), linetype = factor(dtr)))+
  scale_linetype_manual(values = c("longdash", "dotted"))+
  scale_color_manual(values = c("dodgerblue", "purple"))+
  #geom_point()+
  geom_line(size = 1.1)+
  theme_classic()

#lifetime egg production estimates
ggplot(traits.integrated2, aes(x= mean.temp, y = lifetime.eggs, colour = factor(dtr), linetype = factor(dtr)))+
  scale_linetype_manual(values = c("longdash", "dotted"))+
  scale_color_manual(values = c("dodgerblue", "purple"))+
  #geom_point()+
  geom_line(size = 1.1)+
  theme_classic()

