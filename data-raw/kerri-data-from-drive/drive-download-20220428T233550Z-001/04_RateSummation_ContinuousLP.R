
## Kerri Miazgowicz, University of Georgia
## Updated May 07, 2020

## Purpose: Generate rate summation estimates based on thermal performance at constant temperature (Bayesian fits) 
##
## Contents:
# (1) Load data and packages
# (2) Evaluate trait median curve values at 0.1C (what the temp dataframe is)
# (4)a Create a function which rounds the the temperature program df to the nearest 0.1C,
# (4)b Creates a new dataframe which lists trait performance at each time interval
# (4)c Computes the average of the hourly estimates to generate a single value
# (4)d Place these values into a dataframe which consists of 'temp' & 'traitvalue'
# (5) A provisional way to geneate a TPC for the fluctuaiton estimate (Bayesian fit); once I can code the Logan-Parton model, I can incorporate that to generate estimates with fluctuations around all temperatures at 0.1C intervals

##########
###### 1. Load data and packages
##########

#Set working directory
mainDir = "C:/Users/Kerri/Desktop/Chapter2 InProgress"
setwd(mainDir)

# Load libraties for plotting traits
library(plyr) # Slices and dices data
library(plotrix) # For standard error function
library(coda)
library(dplyr)

#Load the fits
load("saved posteriors/temps.Rdata")

#Constant temperature fits
load("saved posteriors/constant_lf_c.quad_T.uniform.Rdata")
# Calculate summary statistics of trait trajectories for plotting
lf.c.quad_T.summary <- trait.trajs.quad(lf.c.quad_T, Temp.xs, summary = FALSE)
load("saved posteriors/constant_br.c.briere_T.uniform.Rdata")
br.c.briere_T.summary <- trait.trajs.briere(br.c.briere_T, Temp.xs, summary = FALSE)
load("saved posteriors/constant_le.c.briere_T.uniform.Rdata")
le.c.briere_T.summary <- trait.trajs.briere(le.c.briere_T, Temp.xs, summary = FALSE)
load("saved posteriors/constant_gamma_c.quad_T.uniform.Rdata")
g.c.quad_T.summary<- trait.trajs.quad(gamma.c.quad_T,Temp.xs, summary = FALSE)


###literature based traits
load("saved posteriors/bc_quad_withT_uniform.Rdata")
bc.c.quad_T.summary <- trait.trajs.quad(bc.c.quad_T,Temp.xs,summary = FALSE)
load("saved posteriors/MDR_briere_withT_uniform.Rdata")
MDR.c.briere_T.summary <- trait.trajs.briere(MDR.c.briere_T, Temp.xs, summary = FALSE)
load("saved posteriors/PDR_briere_withT_uniform.Rdata")
PDR.c.briere_T.summary <- trait.trajs.briere(PDR.c.briere_T, Temp.xs,summary= FALSE)
load("saved posteriors/pea_quad_withT_uniform.Rdata")
pea.c.quad_T.summary <- trait.trajs.quad(pea.c.quad_T, Temp.xs, summary = FALSE)


#Make sure the 03_Calc_R0_constant.R code file was ran and the files are still accessible. (need the matrix object for R0.constant)

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
mean.temp = seq(0,50,.1) #will be a vector of integers from 0 to 45C when implemented (n = 46) #changed to 0.1 for working with rate summation estimates
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

#write.csv(temp.programs.2, "data/temp.programs2.continuous.csv", row.names = FALSE) #list of all temp programs


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
  
  return(output.df) # return output
}

#Note need to transpose all of the summary outputs while using except the prior R0 file
###### Calculate trait posteriors and quantiles
br.c.out <- calcPostQuants(t(br.c.briere_T.summary), Temp.xs)
le.c.out <- calcPostQuants(t(le.c.briere_T.summary), Temp.xs)
lf.c.out <- calcPostQuants(t(lf.c.quad_T.summary), Temp.xs)
g.c.out <- calcPostQuants(t(g.c.quad_T.summary), Temp.xs)
bc.c.out <- calcPostQuants(t(bc.c.quad_T.summary), Temp.xs)
pea.c.out <- calcPostQuants(t(pea.c.quad_T.summary), Temp.xs)
mdr.c.out <- calcPostQuants(t(MDR.c.briere_T.summary), Temp.xs)
R0.c.out <- calcPostQuants(R0.constant, Temp.xs) #Remember no t() here needed, pre-transformed


###### Modify outputs into a format useable for performing rate summation
br.c.out$temp  <-  Temp.xs
le.c.out$temp <- Temp.xs
lf.c.out$temp <-Temp.xs
g.c.out$temp <- Temp.xs
bc.c.out$temp <- Temp.xs
pea.c.out$temp <- Temp.xs
mdr.c.out$temp <- Temp.xs
R0.c.out$temp<- Temp.xs


########## 4
######   (c) Computes the average of the hourly estimates to generate a single value
######   (d) Place these values into a dataframe which consists of 'temp' & 'traitvalue'
##########



#A function which returns the value of the hourly estimate of trait performance based on the values from temperature.program 
# Arguments:
#   x                      a data frame containing the temperature profiles calculated by the LP model
#   temp.performance       a data frame with posterior distributions of trait predicted over a temperature gradient (i.e., 'trait.preds')

hourly.estimates = function(x, temp.performance){
  x<- as.numeric(round(x, 4))
  temp.performance$temp <- as.numeric(round(temp.performance$temp,4))
  median.rate <- temp.performance$median[match(x, temp.performance$temp)]
  return(median.rate)
}

##If temperature is greater or lower than the CTMIN or CTMAX for each trait replace with a negative value (aka -1)
##This way hourly temps function will return an NA Vavlue..
#Thus I'll have a seperate temp.df for each trait....(this will be easy to integrate with the full posteriors....)
br.temp.programs <- le.temp.programs<- lf.temp.programs<- temp.programs.2
g.temp.programs<-bc.temp.programs<-pea.temp.programs<-mdr.temp.programs<-R0.temp.programs<-temp.programs.2
#br.temp.programs[br.temp.programs < 1.2] <- -1  #CTmin 1.2 - don't apply a cutoff to ctmin
br.temp.programs[br.temp.programs >42.1] <- -1 #CTmax 42.1
lf.temp.programs[lf.temp.programs>36.4] <- -1 #CTmax 36.4
le.temp.programs[le.temp.programs> 33.2] <- -1 #CTmax 33.2
g.temp.programs[g.temp.programs>43.2] <- -1
bc.temp.programs[bc.temp.programs>38.1]<- -1
pea.temp.programs[pea.temp.programs>37.1]<- -1
mdr.temp.programs[mdr.temp.programs>35.9]<- -1
R0.temp.programs[R0.temp.programs>38]<- -1

#Uses the function on each trait of interest
br.estimate <- data.frame(apply(br.temp.programs, MARGIN = 2, hourly.estimates, temp.performance = br.c.out))
lf.estimate <-data.frame(apply(lf.temp.programs, MARGIN = 2, hourly.estimates, temp.performance =lf.c.out))
le.estimate <- data.frame(apply(le.temp.programs, MARGIN = 2, hourly.estimates, temp.performance = le.c.out))
g.estimate <- data.frame(apply(temp.programs.2, MARGIN = 2, hourly.estimates, temp.performance = g.c.out))
bc.estimate <- data.frame(apply(temp.programs.2, MARGIN = 2, hourly.estimates, temp.performance = bc.c.out))
pea.estimate <- data.frame(apply(temp.programs.2, MARGIN = 2, hourly.estimates, temp.performance = pea.c.out))
mdr.estimate <- data.frame(apply(temp.programs.2, MARGIN = 2, hourly.estimates, temp.performance = mdr.c.out))
R0.estimate <- data.frame(apply(temp.programs.2, MARGIN = 2, hourly.estimates, temp.performance = R0.c.out))

##Filling in NAs with 0 (aka negative temperatues have no performance values)

##NAs at the low end mean the temp program was negative (0 performance)
##NAs at the high end mean the temp excceded the CTMAX (assign an unique value to use in a conditional statement later; aka keep as NA)
#Leave the NAs in the vector, and when you average, na.rm = FALSE
#br.estimate[is.na(br.estimate)] <- 0
#lf.estimate[is.na(lf.estimate)] <- 0
#le.estimate[is.na(le.estimate)] <- 0
#g.estimate[is.na(g.estimate)] <- 0
#bc.estimate[is.na(bc.estimate)] <- 0
#pea.estimate[is.na(pea.estimate)] <- 0
#mdr.estimate[is.na(mdr.estimate)] <- 0
#R0.estimate[is.na(R0.estimate)] <- 0




##Fix hour column
br.estimate$hr <- hr # hr here should still correlate to the reordered hr value from the exp.temp.programs
lf.estimate$hr <-hr
le.estimate$hr <- hr
g.estimate$hr<- hr
bc.estimate$hr <- hr
pea.estimate$hr<- hr
mdr.estimate$hr <- hr
R0.estimate$hr <- hr


####################
#############
########Visualizing the outputs generated to double-check everything worked as expected

#put all estimates into a list for use of the function below
component.estimation <- list(br.estimate,
                             lf.estimate,
                             le.estimate,
                             g.estimate,
                             bc.estimate,
                             pea.estimate,
                             mdr.estimate,
                             R0.estimate)

plot(br.estimate$hr, br.estimate$X20.dtr12, ylim = c(0,1)) #visualize test -> could plot all the columns together
plot(lf.estimate$hr, lf.estimate$X18.8.dtr9,ylim = c(0,50))
plot(le.estimate$hr, le.estimate$X24.5.dtr9, ylim = c(0,600))

#fix dtr column to remove '.dtr' and change to a numeric data type
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#Takes all of the hourly performance estimates and averages across them
##NOTE: THIS IS WHERE I APPLY THE CTMAX THRESHOLD... changed the na.rm to false
traits.integrated <- data.frame (mean.temp = substr(colnames(br.estimate),2,5),
                                 dtr = sub("^[^.]*", "", substrRight(colnames(br.estimate),6)),
                                 program = colnames(br.estimate),
                                 br = as.vector(colMeans(br.estimate,na.rm = FALSE)),
                                 lf = as.vector(colMeans(lf.estimate,na.rm = FALSE)),
                                 le = as.vector(colMeans(le.estimate,na.rm = FALSE)),
                                 g = as.vector(colMeans(g.estimate,na.rm = FALSE)),
                                 bc = as.vector(colMeans(bc.estimate,na.rm = FALSE)),
                                 pea = as.vector(colMeans(pea.estimate,na.rm = FALSE)),
                                 mdr = as.vector(colMeans(mdr.estimate,na.rm = FALSE)),
                                 R0= as.vector(colMeans(R0.estimate,na.rm = FALSE)))


###Note: If I wanted I could use the above code to employ rate summation on the other components of R0 (MDR, pea, PDR, bc, ect.)
#remove the first row to fix hr issues
traits.integrated <- traits.integrated[-c(1),]
traits.integrated$mean.temp <- seq(0,50,0.1)
#Export dataframe if needed later for use in second bayesian fits (Note: in alternative approach where I estimate around each mean temperature (in 0.1C intervals, I will avoid the need to re-fit the estimations with a Bayesian appraoch, and can instead use the rate summation on the origianl upper and lower CI))
#write.csv(traits.integrated, "data/traits.integrated.CTmaxcutoff.csv", row.names = FALSE)


###########Create a graph for traits.integrated2 to show rate summation estimates for each trait--> will use later when I overlay the bayesian fits of observed constant, dtr9 , and dtr12 data
estimate.plot <- function(traitname, integrated.df,cols ="black", add = "false", labs = ""){
  # Specify the  trait to plot
  par(mar=c(5,6,5,1))
  temp = integrated.df$mean.temp
  traitplot = (integrated.df[ , which(colnames(integrated.df)== traitname)])
  
  # plot the data
  if(add== "false"){
    plot(traitplot ~ temp, xlim = c(0, 50), ylim = c(0, max(traitplot, na.rm = T)*1.02), 
         ylab = labs, xlab = expression(paste("Temperature (",degree,"C)")), cex.lab = 1.4, cex.axis = 1.4, pch = 1, col = cols)
  }else{
    points(traitplot ~ temp, data = df, pch = pch, col = cols)
  }
}

#Replace NA in traits.integrated with 0s for plotting
traits.integrated[is.na(traits.integrated)] <- 0

estimate.plot("br",traits.integrated,"red", labs = "br") #right now shows dtr 9 and 12 at the same time...
estimate.plot("lf", traits.integrated,"blue", labs = "lf")
estimate.plot("le", traits.integrated, "purple", labs ="le")
#plotting function needs work....


