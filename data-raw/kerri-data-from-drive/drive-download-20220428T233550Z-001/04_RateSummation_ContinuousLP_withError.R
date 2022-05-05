
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
load("saved posteriors/temps.Rdata") #501 length

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
load("saved posteriors/MDR_briere_neg_uniform.Rdata")
MDR.c.briere_neg.summary<- trait.trajs.neg.briere(MDR.c.briere_neg, Temp.xs, summary =FALSE)


####Checking how RS differs if the constant temperature TPC is allowed to go negative
load("saved posteriors/pea_neg_quad_uniform.Rdata")
pea.c.quad_neg.summary <- trait.trajs.neg.quad(pea.c.quad_neg, Temp.xs, summary = FALSE)
load("saved posteriors/constant_br.c.neg.test.briere.uniform.Rdata")
br.c.briere_neg.test.summary<-trait.trajs.neg.test.briere(br.c.neg.test.briere, Temp.xs, summary = FALSE)
br.c.neg.test.briere.summary.expanded <- trait.trajs.neg.test.briere(br.c.neg.test.briere, Temp.xs.2, summary = FALSE)

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
mean.temp = seq(0,50,.1) #will be a vector of integers from 0 to 50C when implemented. #changed to 0.1 for working with rate summation estimates
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

###################
#Make a dataframe of only the temperature programs needed for this experiment...
exp.temp.programs <- data.frame(hr = temp.programs.2$hr,
                                X16.dtr9 = temp.programs.2$X16.dtr9,
                                X20.dtr9 =temp.programs.2$X20.dtr9,
                                X24.dtr9 =temp.programs.2$X24.dtr9,
                                X28.dtr9 =temp.programs.2$X28.dtr9,
                                X32.dtr9 =temp.programs.2$X32.dtr9,
                                X16.dtr12 = temp.programs.2$X16.dtr12,
                                X20.dtr12 =temp.programs.2$X20.dtr12,
                                X24.dtr12 =temp.programs.2$X24.dtr12,
                                X28.dtr12 =temp.programs.2$X28.dtr12,
                                X32.dtr12 =temp.programs.2$X32.dtr12)

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

#effect of negative functions
#br.c.neg.out <- calcPostQuants(t(br.c.briere_neg.summary), Temp.xs)
#mdr.c.neg.out <- calcPostQuants(t(MDR.c.briere_neg.summary), Temp.xs)

#pea.c.neg.out <-calcPostQuants(t(pea.c.quad_neg.summary), Temp.xs)
#br.c.neg.test.out <- calcPostQuants(t(br.c.briere_neg.test.summary), Temp.xs)#maybe dont need transformation here
#br.c.neg.test.out.expanded <-calcPostQuants(t(br.c.neg.test.briere.summary.expanded), Temp.xs.2)

###### Modify outputs into a format useable for performing rate summation
br.c.out$temp  <-  Temp.xs
le.c.out$temp <- Temp.xs
lf.c.out$temp <-Temp.xs
g.c.out$temp <- Temp.xs
bc.c.out$temp <- Temp.xs
pea.c.out$temp <- Temp.xs
mdr.c.out$temp <- Temp.xs
R0.c.out$temp<- Temp.xs

#br.c.neg.out$temp <- Temp.xs
#mdr.c.neg.out$temp <- Temp.xs
#pea.c.neg.out$temp <-Temp.xs

#br.c.neg.test.out$temp <- Temp.xs
#br.c.neg.test.out.expanded$temp <- Temp.xs.2

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



##########################################
############################
#Modifying the function to work over all the predictions [1:501, 1:7500]
#Function that converts each mcmc simulation to a dataframe

mcmc.out.func <- function(summary.col){
  temps <- as.numeric(round(seq(0,50,0.1),4)) #can create this externally for use internal?
  mcmc.col.df <- data.frame(temps = temps, sim = summary.col)
  return(mcmc.col.df)
}

##Since the datastructure at the end is changing, I may need to use a for loop and call hourly.performance in it
#instead of temp.performance = br.c.out, use temp.performance = test.fun[[k]]$sim

hourly.estimates3 = function(x, temp.performance){
  x<- as.numeric(round(x, 4))
  temp.performance$temps <- as.numeric(round(temp.performance$temps,4))
  median.rate <- temp.performance$sim[match(x, temp.performance$temps)]
  return(median.rate)
}
################################################################



###################################

#All the traits
#bite rate first
trait.summary <- br.c.briere_T.summary #Matrix [1:501,1:7500] elements
test.fun <- apply(trait.summary, MARGIN = 2, mcmc.out.func) #Generates a list of 7500 elements [containing dfs[501 obs of 2 varibles]; col[1] temps, col[2] perfromance values]
#le,lf,g,bc,pea,mdr,R0
le.fun <- apply(le.c.briere_T.summary, MARGIN = 2, mcmc.out.func) #Generates a list of 7500 elements [containing dfs[501 obs of 2 varibles]; col[1] temps, col[2] perfromance values]
lf.fun <- apply(lf.c.quad_T.summary, MARGIN = 2, mcmc.out.func) #Generates a list of 7500 elements [containing dfs[501 obs of 2 varibles]; col[1] temps, col[2] perfromance values]
g.fun <- apply(g.c.quad_T.summary, MARGIN = 2, mcmc.out.func) #Generates a list of 7500 elements [containing dfs[501 obs of 2 varibles]; col[1] temps, col[2] perfromance values]
bc.fun <- apply(bc.c.quad_T.summary, MARGIN = 2, mcmc.out.func) #Generates a list of 7500 elements [containing dfs[501 obs of 2 varibles]; col[1] temps, col[2] perfromance values]
pea.fun <- apply(pea.c.quad_T.summary, MARGIN = 2, mcmc.out.func) #Generates a list of 7500 elements [containing dfs[501 obs of 2 varibles]; col[1] temps, col[2] perfromance values]
mdr.fun <- apply(MDR.c.briere_T.summary, MARGIN = 2, mcmc.out.func) #Generates a list of 7500 elements [containing dfs[501 obs of 2 varibles]; col[1] temps, col[2] perfromance values]
R0.fun <- apply(t(R0.constant), MARGIN = 2, mcmc.out.func) #Generates a list of 7500 elements [containing dfs[501 obs of 2 varibles]; col[1] temps, col[2] perfromance values]

br.neg.fun <- apply(br.c.briere_neg.summary, MARGIN =2, mcmc.out.func)
mdr.neg.fun <- apply(MDR.c.briere_neg.summary, MARGIN = 2, mcmc.out.func)
pea.neg.fun <- apply(pea.c.quad_neg.summary, MARGIN = 2, mcmc.out.func)

br.neg.test.fun <- apply(br.c.briere_neg.test.summary, MARGIN = 2, mcmc.out.func)

#####If not going negative remove this function
#mcmc.out.func.expanded <- function(summary.col){
#  temps <- as.numeric(round(seq(-6,56,0.1),4)) #can create this externally for use internal?
#  mcmc.col.df <- data.frame(temps = temps, sim = summary.col)
#  return(mcmc.col.df)
#}

#br.neg.test.fun.expanded <- apply(br.c.neg.test.briere.summary.expanded, MARGIN = 2, mcmc.out.func.expanded)
#plot(br.neg.test.fun[[60]]$sim ~br.neg.test.fun[[60]]$temps)

#Remove all large rjags objects from global environment to free up some memory
remove(bc.c.quad_T)
remove(br.c.briere_T)
remove(gamma.c.quad_T)
remove(le.c.briere_T)
remove(lf.c.quad_T)
remove(MDR.c.briere_T)
remove(PDR.c.briere_T)
remove(pea.c.quad_T)
remove(pea.c.quad_neg)
remove(br.c.neg.test.briere)

#initilize item as a list
#This takes a while to run.... (~1 hr run time)
br.estimate2<-NULL
le.estimate2<- lf.estimate2<-g.estimate2<-bc.estimate2<-pea.estimate2<-mdr.estimate2<-R0.estimate<-NULL
pea.estimate2<-NULL

#Save a version of the original for use in later code
temp.programs.negative <- temp.programs.2

#To avoid NA's in the below evaluations change all negative valus in temp.programs.2 to 0
temp.programs.2[temp.programs.2 < 0] <- 0 #I dont think I want to do this when I'm letting the functions go negative

#There will be negative temperature values early on that will be evaluated as Nas for the Briere function....
#I'll need to handle the Na's below the Tmin differently than Na's above the Tmax 

#pea.neg.estimate2<-NULL
#br.neg.estimate2<- NULL
#mdr.neg.estimate2<- NULL
#br.neg.test.estimate2<- NULL


#actually run this seperately for each trait, my computer memory is overwhelmed when I do it all at once
#remove the function after everyblock to free up memory

#Only use 'hourly.estimates3' function if I wanted the functions to go negative
#Else apply 'hourly.estimates' and after the fact cutoff values at the estimated Tmax for that simulation.

for(k in 1:length(pea.neg.fun)){
  pea.neg.estimate2[[k]]<- data.frame(apply(temp.programs.2, MARGIN =2, hourly.estimates3, temp.performance= pea.neg.fun[[k]]))
  print(k)
}
remove(pea.neg.fun)

for(k in 1:length(br.neg.test.fun.expanded)){
  br.neg.test.estimate2[[k]]<- data.frame(apply(temp.programs.2, MARGIN =2, hourly.estimates3, temp.performance= br.neg.test.fun.expanded[[k]]))
  print(k)
}
remove(br.neg.test.fun.expanded)


for(k in 1:length(br.neg.test.fun)){ 
  br.neg.test.estimate2[[k]]<- data.frame(apply(temp.programs.2, MARGIN =2, hourly.estimates3, temp.performance= br.neg.test.fun[[k]]))
  print(k)
}
remove(br.neg.test.fun)


for(k in 1:length(mdr.neg.fun)){
  mdr.neg.estimate2[[k]] <- data.frame(apply(temp.programs.2, MARGIN =2, hourly.estimates3, temp.performance= mdr.neg.fun[[k]]))
  print(k)
  }
remove(mdr.neg.fun)
for(k in 1:length(test.fun)){ 
br.estimate2[[k]] <- data.frame(apply(temp.programs.2, MARGIN = 2, hourly.estimates3, temp.performance = test.fun[[k]]))
print(k)
}
remove(test.fun)

for(k in 1:length(le.fun)){ 
  le.estimate2[[k]] <- data.frame(apply(temp.programs.2, MARGIN = 2, hourly.estimates3, temp.performance = le.fun[[k]]))
   print(k)
}
remove(le.fun)

for(k in 1:length(lf.fun)){ 
  lf.estimate2[[k]] <-data.frame(apply(temp.programs.2, MARGIN = 2, hourly.estimates, temp.performance = lf.fun[[k]]))
   print(k)
}
remove(lf.fun)


####remaining traits needed for S calculation
for(k in 1:length(g.fun)){ 
  g.estimate2[[k]] <-data.frame(apply(temp.programs.2, MARGIN = 2, hourly.estimates3, temp.performance = g.fun[[k]]))
   print(k)
}
remove(g.fun)

for(k in 1:length(bc.fun)){ 
  bc.estimate2[[k]] <-data.frame(apply(temp.programs.2, MARGIN = 2, hourly.estimates3, temp.performance = bc.fun[[k]]))
   print(k)
}
remove(bc.fun)

for(k in 1:length(pea.fun)){ 
  pea.estimate2[[k]] <-data.frame(apply(temp.programs.2, MARGIN = 2, hourly.estimates3, temp.performance = pea.fun[[k]]))
   print(k)
}
remove(pea.fun)

for(k in 1:length(mdr.fun)){ 
  mdr.estimate2[[k]] <-data.frame(apply(temp.programs.2, MARGIN = 2, hourly.estimates3, temp.performance = mdr.fun[[k]]))
   print(k)
}
remove(mdr.fun)

###This one can be run independently.... used in Figure 3
for(k in 1:length(R0.fun)){ 
  R0.estimate[[k]] <-data.frame(apply(temp.programs.2, MARGIN = 2, hourly.estimates3, temp.performance = R0.fun[[k]]))
  print(k)
}
remove(R0.fun)

###################################################
###########################
################
#For each trait
calcThreshPeakDists = function(input, temp.list) {
  
  # Create output dataframe
  output.df <- data.frame("peak" = numeric(nrow(input)), "T0" = numeric(nrow(input)), "Tmax" = numeric(nrow(input)))
  
  for (i in 1:nrow(input)) { # loop through each row of the input (MCMC step)
    
    output.df$peak[i] <- temp.list[which.max(input[i, ])] # Calculate index of R0 peak & store corresponding temperature
    
    # Makes lists of T0 and Tmax
    index.list <- which(input[i, ] > 0) # Create vector of list of indices where R0 > 1
    length.index.list <- length(index.list)
    output.df$T0[i] <- temp.list[index.list[1] - 1] # Store T0 (index prior to first value in index.list)
    output.df$Tmax[i] <- temp.list[index.list[length.index.list] + 1] # Store Tmax (index after to last value in index.list)
  }
  return(output.df) # return
}

####
###################
########################
lf.c.quad_T.preds <- trait.trajs.quad(lf.c.quad_T, Temp.xs, summary = FALSE)
lf.c.TPdists <- calcThreshPeakDists(t(lf.c.quad_T.preds), Temp.xs)
#both 'lf.c.quad_T.preds' & 'lf.c.quad_T.summary (FALSE)' are the same thing

#Cutoff Boolean df function
#A function that makes a df of the LP models on if an indexed simulation exceeds the critial CTmin or CTmax values
#daily.temps is 'temp.programs.negative'
#samps.T.thres is the trait TPdists df calculated from the function 'calcThresPeakDists'
#RETURNS: boolean.df ; 2cols (CTmax, and CTmin)-with 1 indicated that the critical thermal threshold has been exceeded
  # rows correspond to the mcmc sample index

boolean.cutoff.df <- function(daily.temps, samps.T.thres){
  boolean.df<- data.frame(matrix(NA, ncol = ncol(daily.temps), nrow = nrow(samps.T.thres))) #initialize new df
  max.value <- NA
  min.value <- NA
  vector.vals <- NA
  for(i in 1:nrow(samps.T.thres)){ #loops through each posteior sample
    for (k in 2: length(daily.temps)){#loops through each mean temp of the PL model; ignorning the hr col
            ifelse(any(daily.temps > samps.T.thres$Tmax[i]),max.value <- 1, max.value <- 0) #if CTmax is exceeded TRUE, ELSE, FALSE
            ifelse(any(daily.temps < samps.T.thres$T0[i]),min.value <- 1, min.value <- 0) #if CTmin is exceeded TRUE, ELSE, FALSE
            vector.vals  <- boolean.df[i,1]  #, max.value = boolean.df[i,2])
            boolean.df[i,k]<- vector.vals
            print(i)
                   }#end i for loop
        print(k)
            }#end k for loop
    boolean.df[,1]<- NA #ignores the hr colum in the program for all simulations.
  return(boolean.df)
}#end function

lf.boolean.df <- boolean.cutoff.df(temp.programs.negative, lf.c.TPdists)

#A function that uses the Tmax from each posterior distribution sampling; If the PL temperatures exceed that Tmax, than the mean RS estimate is reassigned to be 0 
##Takes in the LP temperature programs
##Takes in TPdists for the posterior distributions
Threshold_cutoff <- function(boolean.cutoff.df){
  
}

####start running code again here....
#so I don't have to rerun the long code if I mess up the list while debugging code
br.estimate3 <- br.estimate2
lf.estimate3<- lf.estimate2

remove(br.estimate2)
#rename br.estimate3 list items names to avoid confusion
names(br.estimate3) <- c(seq(1,7500,1))#c(seq(1,7500,1))
names(lf.estimate2) <- c(seq(1,7500,1))
names(le.estimate2) <-c(seq(1,7500,1))
names(g.estimate2)<- names(bc.estimate2) <- names(pea.estimate2)<- names(mdr.estimate2) <- c(seq(1,7500,1))#c(seq(1,7500,1))


#Try plotting some of this
#plot(hr, br.estimate3$`4`$X1.7.dtr9) #hr by performance for a single mcmc step
plot(hr, br.estimate3[[4]]$X0.9.dtr9) #hr by performance for a single mcmc step
plot(hr, pea.neg.estimate2[[4]]$X7.dtr9)
plot(hr, br.neg.test.estimate2[[4]]$X40.dtr9, ylim= c(-1, .6))

#For each list element [7500], average the values over the temp program columns in the nested df
#Output: List, 1:7500, where each list element contains a df 1obs of 1003 variables; df variable name correponds to the temp.program
pea.neg.mcmc.perf.integration <- lapply(pea.neg.estimate2, FUN = colMeans)
mdr.neg.mcmc.perf.integration <- lapply(mdr.neg.estimate2, FUN = colMeans)
br.mcmc.perf.integration <- lapply(br.estimate3,FUN= colMeans) #but needs to be applied to each list element internal df column

br.neg.mcmc.perf.integration <-lapply(br.neg.test.estimate2, FUN = colMeans)

#save memory by removing the estimate lists after use
remove(br.estimate3)
le.mcmc.perf.rs.integration <- lapply(le.estimate2,FUN= colMeans)
remove(le.estimate2)
lf.mcmc.perf.rs.integration <- lapply(lf.estimate2,FUN= colMeans)
remove(lf.estimate2)
g.mcmc.perf.rs.integration <- lapply(g.estimate2,FUN= colMeans)
remove(g.estimate2)
bc.mcmc.perf.rs.integration <- lapply(bc.estimate2,FUN= colMeans)
remove(bc.estimate2)
pea.mcmc.perf.rs.integration <- lapply(pea.estimate2,FUN= colMeans)
remove(pea.estimate2)
mdr.mcmc.perf.rs.integration <- lapply(mdr.estimate2,FUN= colMeans)
remove(mdr.estimate2)
R0.mcmc.perf.rs.integration <- lapply(R0.estimate,FUN= colMeans)
remove(R0.estimate)

plot(pea.neg.mcmc.perf.integration[[360]])
plot(mdr.neg.mcmc.perf.integration[[360]])
plot(br.neg.mcmc.perf.integration[[360]])


plot(br.mcmc.perf.integration[[200]]) #the outlier at index 1 is hr , which I remove from the dataset later on
plot(br.mcmc.perf.integration[[360]]) #the outlier at index 1 is hr , which I remove from the dataset later on
plot(lf.mcmc.perf.rs.integration[[280]])
plot(le.mcmc.perf.rs.integration[[300]])
plot(R0.mcmc.perf.rs.integration[[300]])
#visualize the output- 1 simulation- 1 dtr set
#set NAs to 0

pea.neg.mcmc.perf.integration <-lapply(pea.neg.mcmc.perf.integration, function(x) replace(x, is.na(x%%1==0),0))
mdr.neg.mcmc.perf.integration <- lapply(mdr.neg.mcmc.perf.integration, function(x) replace (x, is.na(x%%1==0),0))
br.neg.mcmc.perf.integration <- lapply(br.neg.mcmc.perf.integration, function(x) replace (x, is.na(x%%1==0),0))


br.mcmc.perf.integration<-lapply(br.mcmc.perf.integration, function(x) replace(x, is.na(x%%1==0), 0))
le.mcmc.perf.rs.integration<-lapply(le.mcmc.perf.rs.integration, function(x) replace(x, is.na(x%%1==0), 0))
lf.mcmc.perf.rs.integration<-lapply(lf.mcmc.perf.rs.integration, function(x) replace(x, is.na(x%%1==0), 0))
g.mcmc.perf.rs.integration<-lapply(g.mcmc.perf.rs.integration, function(x) replace(x, is.na(x%%1==0), 0))
bc.mcmc.perf.rs.integration<-lapply(bc.mcmc.perf.rs.integration, function(x) replace(x, is.na(x%%1==0), 0))
pea.mcmc.perf.rs.integration<-lapply(pea.mcmc.perf.rs.integration, function(x) replace(x, is.na(x%%1==0), 0))
mdr.mcmc.perf.rs.integration<-lapply(mdr.mcmc.perf.rs.integration, function(x) replace(x, is.na(x%%1==0), 0))

R0.mcmc.perf.rs.integration<-lapply(R0.mcmc.perf.rs.integration, function(x) replace(x, is.na(x%%1==0), 0))

#test plot
plot(Temp.xs,pea.neg.mcmc.perf.integration[[3]][2:502], ylab ="pea RS value dtr 9") #technicall should not be any NAs
plot(Temp.xs,br.neg.mcmc.perf.integration[[3]][2:502])
plot(Temp.xs,br.neg.mcmc.perf.integration[[3]][503:1003])


plot(Temp.xs, mdr.neg.mcmc.perf.integration[[3]][2:502])

plot(Temp.xs,br.mcmc.perf.integration[[3]][2:502], ylab = "br RS value dtr9") #which two colmns correspond with hr?[1,503]
plot(Temp.xs,br.mcmc.perf.integration[[3]][503:1003], ylab = "br RS value dtr12") #which two colmns correspond with hr?[1,503]

###Custom function to return the subsetted dataframe
##Arguments:
#x- a df that the indexing will be performed on [each list element]
#index - a vector containing the start and indexes used to truncate the df #assumes continuious numbering
#returns the subsetted df
df.subset<- function(df, index){
  df.out <- df[index[1]:index[2]]
  return(df.out)
}
#Step2: Take the values for all mcmc sims for each temp program & calculate calcPostQuants to determine the mean,median, and error intervals
#Will have to split up the dtr9 and dtr12 predictions to make it compatible for the calcPostQuants, also may be better for plotting to split up
#2.i - subset the dfs within the list in half 1:501 instead of 1:1003 -> dtr 9 & dtr 12 
pea.neg.dtr9.mcmc.sims <- lapply(pea.neg.mcmc.perf.integration, FUN = df.subset, index = c(2,502))
pea.neg.dtr12.mcmc.sims <- lapply(pea.neg.mcmc.perf.integration, FUN = df.subset, index = c(503,1003))
br.neg.dtr9.mcmc.sims <- lapply(br.neg.mcmc.perf.integration, FUN = df.subset, index = c(2,502)) #change index based on the temp sequence
br.neg.dtr12.mcmc.sims<- lapply(br.neg.mcmc.perf.integration, FUN = df.subset, index = c(503,1003))
mdr.neg.dtr9.mcmc.sims <- lapply(mdr.neg.mcmc.perf.integration, FUN = df.subset, index = c(2,502))
mdr.neg.dtr12.mcmc.sims <- lapply(mdr.neg.mcmc.perf.integration, FUN = df.subset, index = c(503,1003))


br.dtr9.mcmc.sims <- lapply(br.mcmc.perf.integration, FUN = df.subset, index= c(2,502)) #all list elements only internal elements 2:502
#save memory by removing the mcmc.perf.integration files after use
br.dtr12.mcmc.sims <- lapply(br.mcmc.perf.integration, FUN = df.subset, index= c(503,1003))#all list elements only internal elements 503:1003
remove(br.mcmc.perf.integration)

le.dtr9.mcmc.sims <-lapply(le.mcmc.perf.rs.integration, FUN = df.subset, index= c(2,502)) 
le.dtr12.mcmc.sims <-lapply(le.mcmc.perf.rs.integration, FUN = df.subset, index= c(503,1003))
remove(le.mcmc.perf.rs.integration) 
 lf.dtr9.mcmc.sims <-lapply(lf.mcmc.perf.rs.integration, FUN = df.subset, index= c(2,502)) 
  lf.dtr12.mcmc.sims <-lapply(lf.mcmc.perf.rs.integration, FUN = df.subset, index= c(503,1003))
  remove(lf.mcmc.perf.rs.integration)  
  g.dtr9.mcmc.sims <-lapply(g.mcmc.perf.rs.integration, FUN = df.subset, index= c(2,502)) 
  g.dtr12.mcmc.sims <-lapply(g.mcmc.perf.rs.integration, FUN = df.subset, index= c(503,1003))
  remove(g.mcmc.perf.rs.integration)    
  bc.dtr9.mcmc.sims <-lapply(bc.mcmc.perf.rs.integration, FUN = df.subset, index= c(2,502)) 
  bc.dtr12.mcmc.sims <-lapply(bc.mcmc.perf.rs.integration, FUN = df.subset, index= c(503,1003))
  remove(bc.mcmc.perf.rs.integration)       
         pea.dtr9.mcmc.sims <-lapply(pea.mcmc.perf.rs.integration, FUN = df.subset, index= c(2,502)) 
         pea.dtr12.mcmc.sims <-lapply(pea.mcmc.perf.rs.integration, FUN = df.subset, index= c(503,1003))
         remove(pea.mcmc.perf.rs.integration)   
            mdr.dtr9.mcmc.sims <-lapply(mdr.mcmc.perf.rs.integration, FUN = df.subset, index= c(2,502)) 
            mdr.dtr12.mcmc.sims <-lapply(mdr.mcmc.perf.rs.integration, FUN = df.subset, index= c(503,1003))
            remove(mdr.mcmc.perf.rs.integration)    
                
            R0.dtr9.mcmc.sims <-lapply(R0.mcmc.perf.rs.integration, FUN = df.subset, index= c(2,502)) 
                R0.dtr12.mcmc.sims <-lapply(R0.mcmc.perf.rs.integration, FUN = df.subset, index= c(503,1003))

#Step 3: arrange the data to be compatible for use with the calcPostQuants function
#Make 2 new matrixs! (for dtr9 and dtr12) [501] with elements containing a df of all RS values for that temp.program (1 obs of 7500 variables or vice versa) 
###i.e. convert from 7500 (matrix dim 2) ; list elements to 501 (matrix dim 2; internal list objects)
pea.neg.dtr9.matrix.sims <- matrix(unlist(pea.neg.dtr9.mcmc.sims), ncol = 501, byrow =TRUE)
pea.neg.dtr12.matrix.sims <- matrix(unlist(pea.neg.dtr12.mcmc.sims), ncol = 501, byrow =TRUE)
br.neg.dtr9.matrix.sims <- matrix(unlist(br.neg.dtr9.mcmc.sims), ncol =501, byrow =TRUE)
br.neg.dtr12.matrix.sims <- matrix(unlist(br.neg.dtr12.mcmc.sims), ncol =501, byrow =TRUE)
mdr.neg.dtr9.matrix.sims <- matrix(unlist(mdr.neg.dtr9.mcmc.sims), ncol = 501, byrow=TRUE)
mdr.neg.dtr12.matrix.sims <- matrix(unlist(mdr.neg.dtr12.mcmc.sims), ncol = 501, byrow=TRUE)

               
br.dtr9.matrix.sims <- matrix(unlist(br.dtr9.mcmc.sims), ncol = 501,  byrow =TRUE)  #vector 1:2004 [501; so if I populate the matrix by]
remove(br.dtr9.mcmc.sims)
br.dtr12.matrix.sims <- matrix(unlist(br.dtr12.mcmc.sims), ncol  =501, byrow =TRUE)
remove(br.dtr12.mcmc.sims)

le.dtr9.matrix.sims<- matrix(unlist(le.dtr9.mcmc.sims), ncol= 501, byrow=TRUE)
remove(le.dtr9.mcmc.sims)
le.dtr12.matrix.sims<- matrix(unlist(le.dtr12.mcmc.sims), ncol= 501, byrow=TRUE)
remove(le.dtr12.mcmc.sims)
lf.dtr9.matrix.sims<- matrix(unlist(lf.dtr9.mcmc.sims), ncol= 501, byrow=TRUE)
remove(lf.dtr9.mcmc.sims)
lf.dtr12.matrix.sims<- matrix(unlist(lf.dtr12.mcmc.sims), ncol= 501, byrow=TRUE)
remove(lf.dtr12.mcmc.sims)
g.dtr9.matrix.sims<- matrix(unlist(g.dtr9.mcmc.sims), ncol= 501, byrow=TRUE)
remove(g.dtr9.mcmc.sims)
g.dtr12.matrix.sims<- matrix(unlist(g.dtr12.mcmc.sims), ncol= 501, byrow=TRUE)
remove(g.dtr12.mcmc.sims)
bc.dtr9.matrix.sims<- matrix(unlist(bc.dtr9.mcmc.sims), ncol= 501, byrow=TRUE)
remove(bc.dtr9.mcmc.sims)
bc.dtr12.matrix.sims<- matrix(unlist(bc.dtr12.mcmc.sims), ncol= 501, byrow=TRUE)
remove(bc.dtr12.mcmc.sims)

pea.dtr9.matrix.sims<- matrix(unlist(pea.dtr9.mcmc.sims), ncol= 501, byrow=TRUE)
pea.dtr9.matrix.sims[is.na(pea.dtr9.matrix.sims)] = 0

remove(pea.dtr9.mcmc.sims)
pea.dtr12.matrix.sims<- matrix(unlist(pea.dtr12.mcmc.sims), ncol= 501, byrow=TRUE)
remove(pea.dtr12.mcmc.sims)
mdr.dtr9.matrix.sims<- matrix(unlist(mdr.dtr9.mcmc.sims), ncol= 501, byrow=TRUE)
remove(mdr.dtr9.mcmc.sims)
mdr.dtr12.matrix.sims<- matrix(unlist(mdr.dtr12.mcmc.sims), ncol= 501, byrow=TRUE)
remove(mdr.dtr12.mcmc.sims)

plot(pea.dtr9.matrix.sims[2050,])
plot(br.neg.dtr9.matrix.sims[2050,])#The temp evals for the 2050 mcmc posterior
R0.dtr9.matrix.sims<- matrix(unlist(R0.dtr9.mcmc.sims), ncol= 501, byrow=TRUE)
R0.dtr12.matrix.sims<- matrix(unlist(R0.dtr12.mcmc.sims), ncol= 501, byrow=TRUE)


###### calcPostQuants2 - Function to calculate quantiles of derived quantity posteriors (for plotting)
# Arguments:
#   input       data frame with posterior distributions of trait predicted over a temperature gradient (i.e., 'trait.preds')
#   temp.list   the temperature gradient itself (i.e., 'Temp.xs')
calcPostQuants2 = function(input, temp.list) {
  # Get length of gradient
  N.grad.xs <- length(temp.list)
  
  # Create output dataframe
  output.df <- data.frame("mean" = numeric(N.grad.xs), "median" = numeric(N.grad.xs), "lowerCI" = numeric(N.grad.xs), "upperCI" = numeric(N.grad.xs))
  
  # Calculate mean & quantiles
  for(i in 1:N.grad.xs){
    output.df$mean[i] <- mean(input[,i],na.rm = T)
    output.df$lowerCI[i] <- quantile(input[,i], 0.025,na.rm = T)
    output.df$upperCI[i] <- quantile(input[,i], 0.975,na.rm = T)
    output.df$median[i] <- quantile(input[,i], 0.5, na.rm = T)
  }
  return(output.df) # return output
}
pea.neg.dtr9.PostQuants <- calcPostQuants2(pea.neg.dtr9.matrix.sims, Temp.xs)
pea.neg.dtr12.PostQuants <- calcPostQuants2(pea.neg.dtr12.matrix.sims, Temp.xs)
pea.neg.dtr9.PostQuants$temps <- Temp.xs
pea.neg.dtr12.PostQuants$temps <- Temp.xs
br.neg.dtr9.PostQuants <- calcPostQuants2(br.neg.dtr9.matrix.sims, Temp.xs)
br.neg.dtr12.PostQuants <- calcPostQuants2(br.neg.dtr12.matrix.sims, Temp.xs)
br.neg.dtr9.PostQuants$temps <- Temp.xs
br.neg.dtr12.PostQuants$temps <- Temp.xs
mdr.neg.dtr9.PostQuants <- calcPostQuants2(mdr.neg.dtr9.matrix.sims, Temp.xs)
mdr.neg.dtr12.PostQuants <- calcPostQuants2(mdr.neg.dtr12.matrix.sims, Temp.xs)
mdr.neg.dtr9.PostQuants$temps <- Temp.xs
mdr.neg.dtr12.PostQuants$temps <- Temp.xs



br.dtr9.sims.PostQuants <- calcPostQuants2(br.dtr9.matrix.sims, Temp.xs) 
br.dtr12.sims.PostQuants <-calcPostQuants2(br.dtr12.matrix.sims, Temp.xs) 

le.dtr9.sims.PostQuants<- calcPostQuants2(le.dtr9.matrix.sims, Temp.xs)
le.dtr12.sims.PostQuants<- calcPostQuants2(le.dtr12.matrix.sims, Temp.xs)
lf.dtr9.sims.PostQuants<- calcPostQuants2(lf.dtr9.matrix.sims, Temp.xs)
lf.dtr12.sims.PostQuants<- calcPostQuants2(lf.dtr12.matrix.sims, Temp.xs)

###do I need to do this for the remaining traits? YES for plotting over the constant temp TPCs
g.dtr9.sims.PostQuants<- calcPostQuants2(g.dtr9.matrix.sims, Temp.xs)
g.dtr12.sims.PostQuants<- calcPostQuants2(g.dtr12.matrix.sims, Temp.xs)
bc.dtr9.sims.PostQuants<- calcPostQuants2(bc.dtr9.matrix.sims, Temp.xs)
bc.dtr12.sims.PostQuants<- calcPostQuants2(bc.dtr12.matrix.sims, Temp.xs)
pea.dtr9.sims.PostQuants<- calcPostQuants2(pea.dtr9.matrix.sims, Temp.xs)
pea.dtr12.sims.PostQuants<- calcPostQuants2(pea.dtr12.matrix.sims, Temp.xs)
mdr.dtr9.sims.PostQuants<- calcPostQuants2(mdr.dtr9.matrix.sims, Temp.xs)
mdr.dtr12.sims.PostQuants<- calcPostQuants2(mdr.dtr12.matrix.sims, Temp.xs)
R0.dtr9.sims.PostQuants<- calcPostQuants2(R0.dtr9.matrix.sims, Temp.xs)
R0.dtr12.sims.PostQuants<- calcPostQuants2(R0.dtr12.matrix.sims, Temp.xs)

#Step4: Use dtr9.sims.PostQuants to plot mean, and 95%CI for each estimate at 0.1C, for each fluc program dtr9,12
br.dtr9.sims.PostQuants$temps <- Temp.xs
br.dtr12.sims.PostQuants$temps <- Temp.xs
lf.dtr9.sims.PostQuants$temps <- Temp.xs
lf.dtr12.sims.PostQuants$temps <- Temp.xs
le.dtr9.sims.PostQuants$temps <- Temp.xs
le.dtr12.sims.PostQuants$temps <- Temp.xs
#do the same for the rest of the files
g.dtr9.sims.PostQuants$temps<-Temp.xs
g.dtr12.sims.PostQuants$temps<-Temp.xs
bc.dtr9.sims.PostQuants$temps<-Temp.xs
bc.dtr12.sims.PostQuants$temps<-Temp.xs
pea.dtr9.sims.PostQuants$temps<-Temp.xs
pea.dtr12.sims.PostQuants$temps<-Temp.xs
mdr.dtr9.sims.PostQuants$temps<-Temp.xs
mdr.dtr12.sims.PostQuants$temps<-Temp.xs
   
R0.dtr9.sims.PostQuants$temps<- R0.dtr12.sims.PostQuants$temps<- Temp.xs


#br dtr9
plot(Temp.xs,br.mcmc.perf.integration$`3`[2:502],xlim = c(0,50),ylim=c(0,1), ylab = "br RS value dtr9") #which two colmns correspond with hr?[1,503]
points(br.dtr9.sims.PostQuants$temps, br.dtr9.sims.PostQuants$mean, xlim = c(0,50),col='purple')
points(br.dtr9.sims.PostQuants$temps, br.dtr9.sims.PostQuants$median, col='grey')
points(br.dtr9.sims.PostQuants$temps, br.dtr9.sims.PostQuants$upperCI, xlim = c(0,50),col ="red")
points(br.dtr9.sims.PostQuants$temps, br.dtr9.sims.PostQuants$lowerCI, xlim = c(0,50), col ="blue")
#br dtr12
plot(Temp.xs,br.mcmc.perf.integration$`3`[503:1003],xlim = c(0,50),ylim=c(0,1), ylab = "br RS value dtr12") #which two colmns correspond with hr?[1,503]
points(br.dtr12.sims.PostQuants$temps, br.dtr12.sims.PostQuants$mean, col = "purple")
points(br.dtr12.sims.PostQuants$temps, br.dtr12.sims.PostQuants$median, col='grey')
points(br.dtr12.sims.PostQuants$temps, br.dtr12.sims.PostQuants$upperCI, xlim = c(0,50),col ="red")
points(br.dtr12.sims.PostQuants$temps, br.dtr12.sims.PostQuants$lowerCI, xlim = c(0,50), col ="blue")

#great now I'm ready to plot
#le dtr12
plot(Temp.xs,le.mcmc.perf.rs.integration$`3`[503:1003],xlim = c(0,50),ylim=c(0,400), ylab = "le RS value dtr12") #which two colmns correspond with hr?[1,503]
points(le.dtr12.sims.PostQuants$temps, le.dtr12.sims.PostQuants$mean, col = "purple")
points(le.dtr12.sims.PostQuants$temps, le.dtr12.sims.PostQuants$median, col='grey')
points(le.dtr12.sims.PostQuants$temps, le.dtr12.sims.PostQuants$upperCI, xlim = c(0,50),col ="red")
points(le.dtr12.sims.PostQuants$temps, le.dtr12.sims.PostQuants$lowerCI, xlim = c(0,50), col ="blue")



####################################################################################################
####################################################################################################
### FIGURE 2- RS ESTIMATES OVER THE EXPERIMENTAL DATA
####################################################################################################
####################################################################################################
#adjust my previous plotting function
trait.plots.poly2 = function(traitname, summary,df, add = "false", cols = "black", pch = 1, labs = ""){
  temp = df$Treatment
  traitplot = (df[ , which(colnames(df)==traitname)])
  if (add == "false") {
    par(mar=c(5,5,5,1))
    plot(traitplot ~ temp, lty = 1,xlim = c(0,45), xlab = expression(paste("Temperature (",degree,"C)")), ylab = labs, cex.axis = 1.4, cex.lab = 1.4, ylim = c(0, max(traitplot, na.rm = T)*1.2), cex.main = 1.5, col = cols, pch = pch, type="n")
    polygon(c(summary$temp, rev(summary$temp)), c(summary$upper, rev(summary$lower)), col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  } else {
    polygon(c(summary$temp, rev(summary$temp)), c(summary$upper, rev(summary$lower)), col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  }
  #points(traitplot~temp, col = cols, pch = pch)
  lines(summary$temp, summary$mean, lty = 1, lwd = 2, col =  adjustcolor(cols, alpha.f = 1))
}


###Use the code of the observed data for DTR 9 and 12 and then add the RS estimates on top of.
####Main Figure 2a- lifespan, 
# Plot thermal performance curves for observed data over complete RS estimates
pdf("plots/Main_Fig2_lf_overlay_RSwithError.pdf", width = 8, height = 8)
par(mfrow = c(3,2))
trait.plots.poly2("lifespan",lf.f9.quad_T.summary , f9.data.block.temp, expression(paste("Lifespan (days)")), cols = "darkcyan", pch = 16, add = "false")
trait.plots.poly2("lifespan",lf.f12.quad_T.summary , f12.data.block.temp, expression(paste("Lifespan (days)")), cols = "#4a3fa8", pch = 16, add = "true")
lines(lf.dtr9.sims.PostQuants$temps, lf.dtr9.sims.PostQuants$median, col='darkcyan', lty = 5, lwd = 2)
lines(lf.dtr9.sims.PostQuants$temps, lf.dtr9.sims.PostQuants$upperCI, xlim = c(0,50), lty = 3, lwd = 1, col = adjustcolor("darkcyan", alpha.f = 0.5))
lines(lf.dtr9.sims.PostQuants$temps, lf.dtr9.sims.PostQuants$lowerCI, xlim = c(0,50), lty = 3, lwd = 1, col = adjustcolor("darkcyan", alpha.f = 0.5))
polygon(c(lf.dtr9.sims.PostQuants$temps, rev(lf.dtr9.sims.PostQuants$temps)), c(lf.dtr9.sims.PostQuants$upperCI, rev(lf.dtr9.sims.PostQuants$lowerCI)), col = adjustcolor("darkcyan", alpha.f = 0.2), border = NA)
lines(lf.dtr12.sims.PostQuants$temps, lf.dtr12.sims.PostQuants$median, col='#4a3fa8', lty = 5, lwd = 2)
lines(lf.dtr12.sims.PostQuants$temps, lf.dtr12.sims.PostQuants$upperCI, xlim = c(0,50), lty = 3, lwd = 1, col = adjustcolor("#4a3fa8", alpha.f = 0.5))
lines(lf.dtr12.sims.PostQuants$temps, lf.dtr12.sims.PostQuants$lowerCI, xlim = c(0,50), lty = 3, lwd = 1, col = adjustcolor("#4a3fa8", alpha.f = 0.5))
polygon(c(lf.dtr12.sims.PostQuants$temps, rev(lf.dtr12.sims.PostQuants$temps)), c(lf.dtr12.sims.PostQuants$upperCI, rev(lf.dtr12.sims.PostQuants$lowerCI)), col = adjustcolor("#4a3fa8", alpha.f = 0.2), border = NA)
mtext("a", side = 3, at = -5, cex = 1.8, line = 2)
par(mfrow = c(1,1))
dev.off()

####Main Figure 2b- bite rate, 
# Plot thermal performance curves for observed data over complete RS estimates
pdf("plots/Main_Fig2_br_overlay_RSwithError.pdf", width = 8, height = 8)
par(mfrow = c(3,2))
trait.plots.poly2("bite.rate",br.f9.briere_T.summary , f9.data.block.temp, expression(paste("Bite rate (days"^"-1", ")")), cols = "darkcyan", pch = 16, add = "false")
trait.plots.poly2("bite.rate",br.f12.briere_T.summary , f12.data.block.temp, expression(paste("Bite rate (days"^"-1", ")")), cols = "#4a3fa8", pch = 16, add = "true")
lines(br.dtr9.sims.PostQuants$temps, br.dtr9.sims.PostQuants$median, col='darkcyan', lty = 5, lwd = 2)
lines(br.dtr9.sims.PostQuants$temps, br.dtr9.sims.PostQuants$upperCI, xlim = c(0,50), lty = 3, lwd = 1, col = adjustcolor("darkcyan", alpha.f = 0.5))
lines(br.dtr9.sims.PostQuants$temps, br.dtr9.sims.PostQuants$lowerCI, xlim = c(0,50), lty = 3, lwd = 1, col = adjustcolor("darkcyan", alpha.f = 0.5))
polygon(c(br.dtr9.sims.PostQuants$temps, rev(br.dtr9.sims.PostQuants$temps)), c(br.dtr9.sims.PostQuants$upperCI, rev(br.dtr9.sims.PostQuants$lowerCI)), col = adjustcolor("darkcyan", alpha.f = 0.2), border = NA)
lines(br.dtr12.sims.PostQuants$temps, br.dtr12.sims.PostQuants$median, col='#4a3fa8', lty = 5, lwd = 2)
lines(br.dtr12.sims.PostQuants$temps, br.dtr12.sims.PostQuants$upperCI, xlim = c(0,50), lty = 3, lwd = 1, col = adjustcolor("#4a3fa8", alpha.f = 0.5))
lines(br.dtr12.sims.PostQuants$temps, br.dtr12.sims.PostQuants$lowerCI, xlim = c(0,50), lty = 3, lwd = 1, col = adjustcolor("#4a3fa8", alpha.f = 0.5))
polygon(c(br.dtr12.sims.PostQuants$temps, rev(br.dtr12.sims.PostQuants$temps)), c(br.dtr12.sims.PostQuants$upperCI, rev(br.dtr12.sims.PostQuants$lowerCI)), col = adjustcolor("#4a3fa8", alpha.f = 0.2), border = NA)
mtext("b", side = 3, at = -5, cex = 1.8, line = 2)
par(mfrow = c(1,1))
dev.off()

####Main Figure 2c- lifetime eggs, 
# Plot thermal performance curves for observed data over complete RS estimates
pdf("plots/Main_Fig2_le_overlay_RSwithError.pdf", width = 8, height = 8)
par(mfrow = c(3,2))
trait.plots.poly2("lifetime.eggs",le.f9.briere_T.summary , f9.data.block.temp, expression(paste("Lifetime eggs")), cols = "darkcyan", pch = 16, add = "false")
trait.plots.poly2("lifetime.eggs",le.f12.briere_T.summary , f12.data.block.temp, expression(paste("Lifetime eggs")), cols = "#4a3fa8", pch = 16, add = "true")
lines(le.dtr9.sims.PostQuants$temps, le.dtr9.sims.PostQuants$median, col='darkcyan', lty = 5, lwd = 2)
lines(le.dtr9.sims.PostQuants$temps, le.dtr9.sims.PostQuants$upperCI, xlim = c(0,50), lty = 3, lwd = 1, col = adjustcolor("darkcyan", alpha.f = 0.5))
lines(le.dtr9.sims.PostQuants$temps, le.dtr9.sims.PostQuants$lowerCI, xlim = c(0,50), lty = 3, lwd = 1, col = adjustcolor("darkcyan", alpha.f = 0.5))
polygon(c(le.dtr9.sims.PostQuants$temps, rev(le.dtr9.sims.PostQuants$temps)), c(le.dtr9.sims.PostQuants$upperCI, rev(le.dtr9.sims.PostQuants$lowerCI)), col = adjustcolor("darkcyan", alpha.f = 0.2), border = NA)
lines(le.dtr12.sims.PostQuants$temps, le.dtr12.sims.PostQuants$median, col='#4a3fa8', lty = 5, lwd = 2)
lines(le.dtr12.sims.PostQuants$temps, le.dtr12.sims.PostQuants$upperCI, xlim = c(0,50), lty = 3, lwd = 1, col = adjustcolor("#4a3fa8", alpha.f = 0.5))
lines(le.dtr12.sims.PostQuants$temps, le.dtr12.sims.PostQuants$lowerCI, xlim = c(0,50), lty = 3, lwd = 1, col = adjustcolor("#4a3fa8", alpha.f = 0.5))
polygon(c(le.dtr12.sims.PostQuants$temps, rev(le.dtr12.sims.PostQuants$temps)), c(le.dtr12.sims.PostQuants$upperCI, rev(le.dtr12.sims.PostQuants$lowerCI)), col = adjustcolor("#4a3fa8", alpha.f = 0.2), border = NA)
mtext("c", side = 3, at = -5, cex = 1.8, line = 2)
par(mfrow = c(1,1))
dev.off()

trait.plots.poly = function(traitname, summary,df, add = "false", cols = "black", pch = 1, labs = ""){
  temp = df$Treatment
  traitplot = (df[ , which(colnames(df)==traitname)])
  if (add == "false") {
    par(mar=c(5,5,5,1))
    plot(traitplot ~ temp, lty = 1,xlim = c(0,45), xlab = expression(paste("Temperature (",degree,"C)")), ylab = labs, cex.axis = 1.4, cex.lab = 1.4, ylim = c(0, max(traitplot, na.rm = T)*1.5), cex.main = 1.5, col = cols, pch = pch, type="n")
    polygon(c(summary$temp, rev(summary$temp)), c(summary$upper, rev(summary$lower)), col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  } else {
    polygon(c(summary$temp, rev(summary$temp)), c(summary$upper, rev(summary$lower)), col = adjustcolor(cols, alpha.f = 0.2), border = NA)
  }
  points(traitplot~temp, col = cols, pch = pch)
  lines(summary$temp, summary$mean, lty = 1, lwd = 2, col =  adjustcolor(cols, alpha.f = 1))
}

##########Plot RS estimates over the constant temperature TPC
# Plot thermal performance curves with data they were fit over
pdf("plots/allTRAITS_RS_over_constantTPC.pdf", width = 8, height = 10)
#pdf("plots/lf_RS_over_constantTPC.pdf", width = 8, height = 10)
par(mfrow = c(4,3))
trait.plots.poly("lifespan",lf.c.quad_T.summary, c.data.block.temp, expression(paste("Lifespan (days)")), cols = "black", pch = 16, add = "false")
lines(lf.dtr9.sims.PostQuants$temps, lf.dtr9.sims.PostQuants$median, col='darkcyan', lty = 5, lwd = 2)
polygon(c(lf.dtr9.sims.PostQuants$temps, rev(lf.dtr9.sims.PostQuants$temps)), c(lf.dtr9.sims.PostQuants$upperCI, rev(lf.dtr9.sims.PostQuants$lowerCI)), col = adjustcolor("darkcyan", alpha.f = 0.2), border = NA)
lines(lf.dtr12.sims.PostQuants$temps, lf.dtr12.sims.PostQuants$median, col='#4a3fa8', lty = 5, lwd = 2)
polygon(c(lf.dtr12.sims.PostQuants$temps, rev(lf.dtr12.sims.PostQuants$temps)), c(lf.dtr12.sims.PostQuants$upperCI, rev(lf.dtr12.sims.PostQuants$lowerCI)), col = adjustcolor("#4a3fa8", alpha.f = 0.2), border = NA)
mtext("a", side = 3, at = -5, cex = 1.8, line = 2)
#par(mfrow = c(1,1))
#dev.off()

#pdf("plots/br_RS_over_constantTPC.pdf", width = 8, height = 10)
#par(mfrow = c(4,3))
trait.plots.poly("bite.rate",br.c.briere_T.summary, c.data.block.temp, expression(paste("Bite rate (days"^"-1", ")")), cols = "black", pch = 16, add = "false")
lines(br.dtr9.sims.PostQuants$temps, br.dtr9.sims.PostQuants$median, col='darkcyan', lty = 5, lwd = 2)
polygon(c(br.dtr9.sims.PostQuants$temps, rev(br.dtr9.sims.PostQuants$temps)), c(br.dtr9.sims.PostQuants$upperCI, rev(br.dtr9.sims.PostQuants$lowerCI)), col = adjustcolor("darkcyan", alpha.f = 0.2), border = NA)
lines(br.dtr12.sims.PostQuants$temps, br.dtr12.sims.PostQuants$median, col='#4a3fa8', lty = 5, lwd = 2)
polygon(c(br.dtr12.sims.PostQuants$temps, rev(br.dtr12.sims.PostQuants$temps)), c(br.dtr12.sims.PostQuants$upperCI, rev(br.dtr12.sims.PostQuants$lowerCI)), col = adjustcolor("#4a3fa8", alpha.f = 0.2), border = NA)
mtext("b", side = 3, at = -5, cex = 1.8, line = 2)
#par(mfrow = c(1,1))
#dev.off()

#summary needs to be true ## plot for revised br RS with negative c.TPC
br.c.briere_neg.summary2 <- trait.trajs.neg.test.briere(br.c.neg.test.briere, Temp.xs, summary = TRUE)
plot( br.c.briere_neg.summary2$median~br.c.briere_neg.summary2$temp, ylab =expression(paste("Bite rate (a)")), col = "black", pch = 16, ylim = c(0,0.5), yaxs="i",type = "n")
lines(br.c.briere_neg.summary2$median~br.c.briere_neg.summary2$temp, col = "black", lty = 1, lwd =2)
lines(br.neg.dtr9.PostQuants$temps, br.neg.dtr9.PostQuants$median, col='darkcyan', lty = 5, lwd = 2)
polygon(c(br.neg.dtr9.PostQuants$temps, rev(br.neg.dtr9.PostQuants$temps)), c(br.neg.dtr9.PostQuants$upperCI, rev(br.neg.dtr9.PostQuants$lowerCI)), col = adjustcolor("darkcyan", alpha.f = 0.2), border = NA)
lines(br.neg.dtr12.PostQuants$temps, br.neg.dtr12.PostQuants$median, col='#4a3fa8', lty = 5, lwd = 2)
polygon(c(br.neg.dtr12.PostQuants$temps, rev(br.neg.dtr12.PostQuants$temps)), c(br.neg.dtr12.PostQuants$upperCI, rev(br.neg.dtr12.PostQuants$lowerCI)), col = adjustcolor("#4a3fa8", alpha.f = 0.2), border = NA)
mtext("b", side = 3, at = -5, cex = 1.8, line = 2)
#Ok but why are these values converging to the same Tmax with the estimates?
#Add in lines for the observed data fits
lines(br.f9.briere_T.summary$median ~br.f9.briere_T.summary$temp, col ="darkcyan", lwd =2)
lines(br.f12.briere_T.summary$median ~br.f12.briere_T.summary$temp, col ="#4a3fa8", lwd =2)

br.neg.dtr9.PostQuants$median
br.neg.dtr9.PostQuants$temps[10] #min values 


#pdf("plots/le_RS_over_constantTPC.pdf", width = 8, height = 10)
#par(mfrow = c(4,3))
trait.plots.poly("lifetime.eggs",le.c.briere_T.summary, c.data.block.temp, expression(paste("Lifetime eggs")), cols = "black", pch = 16, add = "false")
lines(le.dtr9.sims.PostQuants$temps, le.dtr9.sims.PostQuants$median, col='darkcyan', lty = 5, lwd = 2)
polygon(c(le.dtr9.sims.PostQuants$temps, rev(le.dtr9.sims.PostQuants$temps)), c(le.dtr9.sims.PostQuants$upperCI, rev(le.dtr9.sims.PostQuants$lowerCI)), col = adjustcolor("darkcyan", alpha.f = 0.2), border = NA)
lines(le.dtr12.sims.PostQuants$temps, le.dtr12.sims.PostQuants$median, col='#4a3fa8', lty = 5, lwd = 2)
polygon(c(le.dtr12.sims.PostQuants$temps, rev(le.dtr12.sims.PostQuants$temps)), c(le.dtr12.sims.PostQuants$upperCI, rev(le.dtr12.sims.PostQuants$lowerCI)), col = adjustcolor("#4a3fa8", alpha.f = 0.2), border = NA)
mtext("c", side = 3, at = -5, cex = 1.8, line = 2)
#par(mfrow = c(1,1))
#dev.off()

#pdf("plots/g_RS_over_constantTPC.pdf", width = 8, height = 10)
#par(mfrow = c(4,3))
trait.plots.poly("gamma",g.c.quad_T.summary, c.gamma.data, expression(paste("Gamma")), cols = "black", pch = 16, add = "false")
lines(g.dtr9.sims.PostQuants$temps, g.dtr9.sims.PostQuants$median, col='darkcyan', lty = 5, lwd = 2)
polygon(c(g.dtr9.sims.PostQuants$temps, rev(g.dtr9.sims.PostQuants$temps)), c(g.dtr9.sims.PostQuants$upperCI, rev(g.dtr9.sims.PostQuants$lowerCI)), col = adjustcolor("darkcyan", alpha.f = 0.2), border = NA)
lines(g.dtr12.sims.PostQuants$temps, g.dtr12.sims.PostQuants$median, col='#4a3fa8', lty = 5, lwd = 2)
polygon(c(g.dtr12.sims.PostQuants$temps, rev(g.dtr12.sims.PostQuants$temps)), c(g.dtr12.sims.PostQuants$upperCI, rev(g.dtr12.sims.PostQuants$lowerCI)), col = adjustcolor("#4a3fa8", alpha.f = 0.2), border = NA)
mtext("d", side = 3, at = -5, cex = 1.8, line = 2)
#par(mfrow = c(1,1))
#dev.off()

#pdf("plots/bc_RS_over_constantTPC.pdf", width = 8, height = 10)
#par(mfrow = c(4,3))
trait.plots.poly("bc",bc.c.quad_T.summary, bc.PDR.data, expression(paste("Vector competence")), cols = "black", pch = 16, add = "false")
lines(bc.dtr9.sims.PostQuants$temps, bc.dtr9.sims.PostQuants$median, col='darkcyan', lty = 5, lwd = 2)
polygon(c(bc.dtr9.sims.PostQuants$temps, rev(bc.dtr9.sims.PostQuants$temps)), c(bc.dtr9.sims.PostQuants$upperCI, rev(bc.dtr9.sims.PostQuants$lowerCI)), col = adjustcolor("darkcyan", alpha.f = 0.2), border = NA)
lines(bc.dtr9.sims.PostQuants$temps, bc.dtr9.sims.PostQuants$median, col='#4a3fa8', lty = 5, lwd = 2)
polygon(c(bc.dtr9.sims.PostQuants$temps, rev(bc.dtr9.sims.PostQuants$temps)), c(bc.dtr12.sims.PostQuants$upperCI, rev(bc.dtr12.sims.PostQuants$lowerCI)), col = adjustcolor("#4a3fa8", alpha.f = 0.2), border = NA)
mtext("e", side = 3, at = -5, cex = 1.8, line = 2)
#par(mfrow = c(1,1))
#dev.off()
##########Calculate the TPdists for the RS estimates

#add plot for pea
#pdf("plots/pea_RS_over_constantTPC.pdf", width = 8, height = 10)
#par(mfrow = c(4,3))
trait.plots.poly("Pea",pea.c.quad_T.summary, pea.MDR.data, expression(paste("Prob. egg to adult survival")), cols = "black", pch = 16, add = "false")
lines(pea.dtr9.sims.PostQuants$temps, pea.dtr9.sims.PostQuants$median, col='darkcyan', lty = 5, lwd = 2)
polygon(c(pea.dtr9.sims.PostQuants$temps, rev(pea.dtr9.sims.PostQuants$temps)), c(pea.dtr9.sims.PostQuants$upperCI, rev(pea.dtr9.sims.PostQuants$lowerCI)), col = adjustcolor("darkcyan", alpha.f = 0.2), border = NA)
lines(pea.dtr12.sims.PostQuants$temps, pea.dtr12.sims.PostQuants$median, col='#4a3fa8', lty = 5, lwd = 2)
polygon(c(pea.dtr12.sims.PostQuants$temps, rev(pea.dtr12.sims.PostQuants$temps)), c(pea.dtr12.sims.PostQuants$upperCI, rev(pea.dtr12.sims.PostQuants$lowerCI)), col = adjustcolor("#4a3fa8", alpha.f = 0.2), border = NA)
mtext("f", side = 3, at = -5, cex = 1.8, line = 2)
#par(mfrow = c(1,1))
#dev.off()
#add plot for mdr

#summary needs to be true ## plot for revised pea RS with negative c.TPC
pea.c.quad_neg.summary2 <- trait.trajs.neg.quad(pea.c.quad_T_neg, Temp.xs, summary = TRUE)
trait.plots.poly2("Pea",pea.c.quad_neg.summary2, pea.MDR.data,expression(paste("Prob. egg to adult survival")), cols = "black", pch = 16, add = "false")
lines(pea.neg.dtr9.PostQuants$temps, pea.neg.dtr9.PostQuants$median, col='darkcyan', lty = 5, lwd = 2)
polygon(c(pea.neg.dtr9.PostQuants$temps, rev(pea.neg.dtr9.PostQuants$temps)), c(pea.neg.dtr9.PostQuants$upperCI, rev(pea.neg.dtr9.PostQuants$lowerCI)), col = adjustcolor("darkcyan", alpha.f = 0.2), border = NA)
lines(pea.neg.dtr12.PostQuants$temps, pea.neg.dtr12.PostQuants$median, col='#4a3fa8', lty = 5, lwd = 2)
polygon(c(pea.neg.dtr12.PostQuants$temps, rev(pea.neg.dtr12.PostQuants$temps)), c(pea.neg.dtr12.PostQuants$upperCI, rev(pea.neg.dtr12.PostQuants$lowerCI)), col = adjustcolor("#4a3fa8", alpha.f = 0.2), border = NA)
mtext("f", side = 3, at = -5, cex = 1.8, line = 2)


#pdf("plots/MDR_RS_over_constantTPC.pdf", width = 8, height = 10)
#par(mfrow = c(4,3))
trait.plots.poly("MDR",MDR.c.briere_T.summary, pea.MDR.data, expression(paste("Mosquito development rate")), cols = "black", pch = 16, add = "false")
lines(mdr.dtr9.sims.PostQuants$temps, mdr.dtr9.sims.PostQuants$median, col='darkcyan', lty = 5, lwd = 2)
polygon(c(mdr.dtr9.sims.PostQuants$temps, rev(mdr.dtr9.sims.PostQuants$temps)), c(mdr.dtr9.sims.PostQuants$upperCI, rev(mdr.dtr9.sims.PostQuants$lowerCI)), col = adjustcolor("darkcyan", alpha.f = 0.2), border = NA)
lines(mdr.dtr12.sims.PostQuants$temps, mdr.dtr12.sims.PostQuants$median, col='#4a3fa8', lty = 5, lwd = 2)
polygon(c(mdr.dtr12.sims.PostQuants$temps, rev(mdr.dtr12.sims.PostQuants$temps)), c(mdr.dtr12.sims.PostQuants$upperCI, rev(mdr.dtr12.sims.PostQuants$lowerCI)), col = adjustcolor("#4a3fa8", alpha.f = 0.2), border = NA)
mtext("g", side = 3, at = -5, cex = 1.8, line = 2)
par(mfrow = c(1,1))
dev.off()
#add plot for R0 directly

#summary needs to be true ## plot for revised pea RS with negative c.TPC
mdr.c.briere_neg.summary2 <- trait.trajs.neg.briere(MDR.c.briere_neg, Temp.xs, summary = TRUE)
plot( mdr.c.briere_neg.summary2$mean~mdr.c.briere_neg.summary2$temp, ylab ="MDR",ylim=c(0,.2), col = "black", pch = 16)
lines(mdr.neg.dtr9.PostQuants$temps, mdr.neg.dtr9.PostQuants$median, col='darkcyan', lty = 5, lwd = 2)
polygon(c(mdr.neg.dtr9.PostQuants$temps, rev(mdr.neg.dtr9.PostQuants$temps)), c(mdr.neg.dtr9.PostQuants$upperCI, rev(mdr.neg.dtr9.PostQuants$lowerCI)), col = adjustcolor("darkcyan", alpha.f = 0.2), border = NA)
lines(mdr.neg.dtr12.PostQuants$temps, mdr.neg.dtr12.PostQuants$median, col='#4a3fa8', lty = 5, lwd = 2)
polygon(c(mdr.neg.dtr12.PostQuants$temps, rev(mdr.neg.dtr12.PostQuants$temps)), c(mdr.neg.dtr12.PostQuants$upperCI, rev(mdr.neg.dtr12.PostQuants$lowerCI)), col = adjustcolor("#4a3fa8", alpha.f = 0.2), border = NA)
mtext("b", side = 3, at = -5, cex = 1.8, line = 2)



##########
###### 4. Function to calculate distributions of T0, Tmax, and peak R0
##########
# Arguments:
#   input       data frame with posterior distributions of R0 (or other quantity) predicted over a temperature gradient
#           input structure [:,:]
#   temp.list   the temperature gradient itself (i.e., 'Temp.xs')

calcThreshPeakDists = function(input, temp.list) {
  
  # Create output dataframe
  output.df <- data.frame("peak" = numeric(nrow(input)), "T0" = numeric(nrow(input)), "Tmax" = numeric(nrow(input)))
  
  for (i in 1:nrow(input)) { # loop through each row of the input (MCMC step)
    
    output.df$peak[i] <- temp.list[which.max(input[i, ])] # Calculate index of R0 peak & store corresponding temperature
    
    # Makes lists of T0 and Tmax
    index.list <- which(input[i, ] > 0) # Create vector of list of indices where R0 > 0
    length.index.list <- length(index.list)
    #create a condition for if index.list includes the first temperatue element
    ifelse(1 %in% index.list,
           output.df$T0[i] <- temp.list[index.list[1]], # Store T0 (index prior to first value in index.list)
           output.df$T0[i] <- temp.list[index.list[1] - 1] # Store T0 (index prior to first value in index.list)
          )
    output.df$Tmax[i] <- temp.list[index.list[length.index.list] + 1] # Store Tmax (index after to last value in index.list)
    #print(i)
    }
  return(output.df) # return
}

#bite rate
br.RS.dtr9.TPdists <- calcThreshPeakDists(br.dtr9.matrix.sims, Temp.xs)
br.RS.dtr12.TPdists <- calcThreshPeakDists(br.dtr12.matrix.sims, Temp.xs)

lf.RS.dtr9.TPdists <- calcThreshPeakDists(lf.dtr9.matrix.sims, Temp.xs)
lf.RS.dtr12.TPdists <- calcThreshPeakDists(lf.dtr12.matrix.sims, Temp.xs)

le.RS.dtr9.TPdists <- calcThreshPeakDists(le.dtr9.matrix.sims, Temp.xs)
le.RS.dtr12.TPdists <- calcThreshPeakDists(le.dtr12.matrix.sims, Temp.xs)

R0.RS.dtr9.direct.TPdists <- calcThreshPeakDists(R0.dtr9.matrix.sims, Temp.xs)
R0.RS.dtr12.direct.TPdists <- calcThreshPeakDists(R0.dtr12.matrix.sims, Temp.xs)


break.n <- seq(0,50,0.4)
# br- Tmin
pRS1 <- hist(br.RS.dtr9.TPdists$T0, breaks = break.n,freq= F)
pRS2 <- hist(br.RS.dtr12.TPdists$T0, breaks = break.n,freq= F)
#br-peak
pRS3 <- hist(br.RS.dtr9.TPdists$peak, breaks = break.n,freq= F)
pRS4 <- hist(br.RS.dtr12.TPdists$peak, breaks = break.n,freq= F)
#br-Tmax
pRS5 <- hist(br.RS.dtr9.TPdists$Tmax, breaks = break.n,freq= F)
pRS6 <- hist(br.RS.dtr12.TPdists$Tmax, breaks = break.n,freq= F)

# lf- Tmin
pRS7 <- hist(lf.RS.dtr9.TPdists$T0, breaks = break.n,freq= F)
pRS8 <- hist(lf.RS.dtr12.TPdists$T0, breaks = break.n,freq= F)
#lf-peak
pRS9 <- hist(lf.RS.dtr9.TPdists$peak, breaks = break.n,freq= F)
pRS10 <- hist(lf.RS.dtr12.TPdists$peak, breaks = break.n,freq= F)
#lf-Tmax
pRS11 <- hist(lf.RS.dtr9.TPdists$Tmax, breaks = break.n,freq= F)
pRS12 <- hist(lf.RS.dtr12.TPdists$Tmax, breaks = break.n,freq= F)

# le- Tmin
pRS13 <- hist(le.RS.dtr9.TPdists$T0, breaks = break.n,freq= F)
pRS14 <- hist(le.RS.dtr12.TPdists$T0, breaks = break.n,freq= F)
#le-peak
pRS15 <- hist(le.RS.dtr9.TPdists$peak, breaks = break.n,freq= F)
pRS16 <- hist(le.RS.dtr12.TPdists$peak, breaks = break.n,freq= F)
#le-Tmax
pRS17 <- hist(le.RS.dtr9.TPdists$Tmax, breaks = break.n,freq= F)
pRS18 <- hist(le.RS.dtr12.TPdists$Tmax, breaks = break.n,freq= F)


# S- Tmin
pRS19 <- hist(R0.RS.dtr9.allpost.TPdists$T0, breaks = break.n,freq= F)
pRS20 <- hist(R0.RS.dtr12.allpost.TPdists$T0, breaks = break.n,freq= F)
#S-peak
pRS21 <- hist(R0.RS.dtr9.allpost.TPdists$peak, breaks = break.n,freq= F)
pRS22 <- hist(R0.RS.dtr12.allpost.TPdists$peak, breaks = break.n,freq= F)
#S-Tmax
pRS23 <- hist(R0.RS.dtr9.allpost.TPdists$Tmax, breaks = break.n,freq= F)
pRS24 <- hist(R0.RS.dtr12.allpost.TPdists$Tmax, breaks = break.n,freq= F)





#########Now plot with the constant temperature density plots

library(plotrix)
###############################Ok make a seperate plot for each trait
pdf("plots/Thresholds_density_RS_lf_a_B_S.pdf", width = 12, height = 12) #width = 8, height = 16) #alter dims to fit onto a single page// 4 by 8/// 2 by 4// 1 x 2//
#pdf("plots/Thresholds_density_RS_lf.pdf", width = 12, height = 12) #width = 8, height = 16) #alter dims to fit onto a single page// 4 by 8/// 2 by 4// 1 x 2//
par(mar=c(5,6,5,1))
par(mfrow = c(4,3))
#lifespan -Tmin
plot(p2, col= scales::alpha('darkcyan',0.4),freq= F,border=scales::alpha('darkcyan',0.6),xlim= c(0,15),ylim = c(0,0.5), main = expression(italic(paste("T"[min])-lf)),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")),cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5 )
plot(p3, col= scales::alpha('#4a3fa8',0.4),freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
plot(pRS7, col= scales::alpha('darkcyan',0.4),angle = 45,density=25,freq= F,border=scales::alpha('darkcyan',0.6), ylab = "Density",add = T)
plot(pRS8, col= scales::alpha('#4a3fa8',0.4),angle = -45,density=25,freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
legend("topright", legend = c("data: DTR 9C", "data: DTR 12C","RS: DTR 9C", "RS: DTR 12C"),
       fill= c(scales::alpha('darkcyan',0.4),alpha('#4a3fa8',0.4),"darkcyan","#4a3fa8"),
       #c(scales::alpha('darkcyan',0.4),scales::alpha('#4a3fa8',0.4), "white","white"),
       angle= c(NA,NA,45,-45),
       density = c(NA,NA,25,25),
       border = c("darkcyan","#4a3fa8", "darkcyan","#4a3fa8"),
       pch = NULL,
       bty = 'n')
mtext("a", side = 3, at = -2.5, cex = 1.8, line = 2)
#peak
plot(p5, col= scales::alpha('darkcyan',0.4),freq= F,border=scales::alpha('darkcyan',0.6),xlim= c(15,24),ylim = c(0,0.5),main = expression(italic(paste("T"[opt])-lf)),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")),cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5 )
plot(p6, col= scales::alpha('#4a3fa8',0.4),freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
plot(pRS9, col= scales::alpha('darkcyan',0.4),angle = 45,density=25,freq= F,border=scales::alpha('darkcyan',0.6), ylab = "Density",add = T)
plot(pRS10, col= scales::alpha('#4a3fa8',0.4),angle = -45,density=25,freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
#Tmax
plot(p8, col= scales::alpha('darkcyan',0.4),freq= F,border=scales::alpha('darkcyan',0.6),xlim= c(24,45),ylim = c(0,0.4),main = expression(italic(paste("T"[max])-lf)),xlab = expression(paste("Temperature (",degree,"C)")),cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5  ,ylab = "Density",cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5 )
plot(p9, col= scales::alpha('#4a3fa8',0.4),freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
plot(pRS11, col= scales::alpha('darkcyan',0.4),angle = 45,density=25,freq= F,border=scales::alpha('darkcyan',0.6), ylab = "Density",add = T)
plot(pRS12, col= scales::alpha('#4a3fa8',0.4),angle = -45,density=25,freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
#par(mfrow = c(1,1))
#dev.off()

################################3
#pdf("plots/Thresholds_density_RS_a.pdf", width = 12, height = 12) #width = 8, height = 16) #alter dims to fit onto a single page// 4 by 8/// 2 by 4// 1 x 2//
#par(mar=c(5,6,5,1))
#par(mfrow = c(4,3))
#bite rate
plot(p11, col= scales::alpha('darkcyan',0.4),freq= F,border=scales::alpha('darkcyan',0.6),xlim= c(0,8),ylim = c(0,1), main = expression(italic(paste("T"[min])-a)),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")),cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5 )
plot(p12, col= scales::alpha('#4a3fa8',0.4),freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
plot(pRS1, col= scales::alpha('darkcyan',0.4),angle = 45,density=25,freq= F,border=scales::alpha('darkcyan',0.6), ylab = "Density",add = T)
plot(pRS2, col= scales::alpha('#4a3fa8',0.4),angle = -45,density=25,freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
legend("topright", legend = c("data: DTR 9C", "data: DTR 12C","RS: DTR 9C", "RS: DTR 12C"),
                  fill= c(scales::alpha('darkcyan',0.4),alpha('#4a3fa8',0.4),"darkcyan","#4a3fa8"),
                  #c(scales::alpha('darkcyan',0.4),scales::alpha('#4a3fa8',0.4), "white","white"),
                  angle= c(NA,NA,45,-45),
                  density = c(NA,NA,25,25),
                  border = c("darkcyan","#4a3fa8", "darkcyan","#4a3fa8"),
                  pch = NULL,
                  bty = 'n')
mtext("b", side = 3, at = -1.25, cex = 1.8, line = 2)
#peak
plot(p14, col= scales::alpha('darkcyan',0.4),freq= F,border=scales::alpha('darkcyan',0.6),xlim= c(26,36),ylim = c(0,0.6),main = expression(italic(paste("T"[opt])-a)),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")),cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5 )
plot(p15, col= scales::alpha('#4a3fa8',0.4),freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
plot(pRS3, col= scales::alpha('darkcyan',0.4),angle = 45,density=25,freq= F,border=scales::alpha('darkcyan',0.6), ylab = "Density",add = T)
plot(pRS4, col= scales::alpha('#4a3fa8',0.4),angle = -45,density=25,freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
#Tmax
plot(p17, col= scales::alpha('darkcyan',0.4),freq= F,border=scales::alpha('darkcyan',0.6),xlim= c(32,46),ylim = c(0,0.6),main = expression(italic(paste("T"[max])-a)),xlab = expression(paste("Temperature (",degree,"C)")),cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5  ,ylab = "Density",cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5 )
plot(p18, col= scales::alpha('#4a3fa8',0.4),freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
plot(pRS5, col= scales::alpha('darkcyan',0.4),angle = 45,density=25,freq= F,border=scales::alpha('darkcyan',0.6), ylab = "Density",add = T)
plot(pRS6, col= scales::alpha('#4a3fa8',0.4),angle = -45,density=25,freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
#par(mfrow = c(1,1))
#dev.off()
########################################LIFESPAN


########################################LIFetime eggs
#pdf("plots/Thresholds_density_RS_le.pdf", width = 12, height = 12) #width = 8, height = 16) #alter dims to fit onto a single page// 4 by 8/// 2 by 4// 1 x 2//
#par(mar=c(5,6,5,1))
#par(mfrow = c(4,3))
#lifetime eggs -Tmin
plot(p21, col= scales::alpha('darkcyan',0.4),freq= F,border=scales::alpha('darkcyan',0.6),xlim= c(0,20),ylim = c(0,0.2), main = expression(italic(paste("T"[min])-B)),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")),cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5 )
plot(p22, col= scales::alpha('#4a3fa8',0.4),freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
plot(pRS13, col= scales::alpha('darkcyan',0.4),angle = 45,density=25,freq= F,border=scales::alpha('darkcyan',0.6), ylab = "Density",add = T)
plot(pRS14, col= scales::alpha('#4a3fa8',0.4),angle = -45,density=25,freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
legend("topright", legend = c("data: DTR 9C", "data: DTR 12C","RS: DTR 9C", "RS: DTR 12C"),
       fill= c(scales::alpha('darkcyan',0.4),alpha('#4a3fa8',0.4),"darkcyan","#4a3fa8"),
       #c(scales::alpha('darkcyan',0.4),scales::alpha('#4a3fa8',0.4), "white","white"),
       angle= c(NA,NA,45,-45),
       density = c(NA,NA,25,25),
       border = c("darkcyan","#4a3fa8", "darkcyan","#4a3fa8"),
       pch = NULL,
       bty = 'n')
mtext("c", side = 3, at = -4, cex = 1.8, line = 2)
#peak
plot(p23, col= scales::alpha('darkcyan',0.4),freq= F,border=scales::alpha('darkcyan',0.6),xlim= c(20,32),ylim = c(0,0.5),main = expression(italic(paste("T"[opt])-B)),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")),cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5 )
plot(p24, col= scales::alpha('#4a3fa8',0.4),freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
plot(pRS15, col= scales::alpha('darkcyan',0.4),angle = 45,density=25,freq= F,border=scales::alpha('darkcyan',0.6), ylab = "Density",add = T)
plot(pRS16, col= scales::alpha('#4a3fa8',0.4),angle = -45,density=25,freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
#Tmax
plot(p8, col= scales::alpha('darkcyan',0.4),freq= F,border=scales::alpha('darkcyan',0.6),xlim= c(25,45),ylim = c(0,0.4),main = expression(italic(paste("T"[max])-B)),xlab = expression(paste("Temperature (",degree,"C)")),cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5  ,ylab = "Density",cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5 )
plot(p9, col= scales::alpha('#4a3fa8',0.4),freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
plot(pRS17, col= scales::alpha('darkcyan',0.4),angle = 45,density=25,freq= F,border=scales::alpha('darkcyan',0.6), ylab = "Density",add = T)
plot(pRS18, col= scales::alpha('#4a3fa8',0.4),angle = -45,density=25,freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
#par(mfrow = c(1,1))
#dev.off()


########################################LIFetime eggs
#pdf("plots/Thresholds_density_RS_S.pdf", width = 12, height = 12) #width = 8, height = 16) #alter dims to fit onto a single page// 4 by 8/// 2 by 4// 1 x 2//
#par(mar=c(5,6,5,1))
#par(mfrow = c(4,3))
#Suitability -Tmin
plot(p31, col= scales::alpha('darkcyan',0.4),freq= F,border=scales::alpha('darkcyan',0.6),xlim= c(8,20),ylim = c(0,1), main = expression(italic(paste("T"[min])-S)),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")),cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5 )
plot(p32, col= scales::alpha('#4a3fa8',0.4),freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
plot(pRS19, col= scales::alpha('darkcyan',0.4),angle = 45,density=25,freq= F,border=scales::alpha('darkcyan',0.6), ylab = "Density",add = T)
plot(pRS20, col= scales::alpha('#4a3fa8',0.4),angle = -45,density=25,freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
legend("topright", legend = c("data: DTR 9C", "data: DTR 12C","RS: DTR 9C", "RS: DTR 12C"),
       fill= c(scales::alpha('darkcyan',0.4),alpha('#4a3fa8',0.4),"darkcyan","#4a3fa8"),
       #c(scales::alpha('darkcyan',0.4),scales::alpha('#4a3fa8',0.4), "white","white"),
       angle= c(NA,NA,45,-45),
       density = c(NA,NA,25,25),
       border = c("darkcyan","#4a3fa8", "darkcyan","#4a3fa8"),
       pch = NULL,
       bty = 'n')
mtext("d", side = 3, at = 6, cex = 1.8, line = 2)
#suitabilit-peak
plot(p33, col= scales::alpha('darkcyan',0.4),freq= F,border=scales::alpha('darkcyan',0.6),xlim= c(24,28),ylim = c(0,1),main = expression(italic(paste("T"[opt])-S)),ylab = "Density",xlab = expression(paste("Temperature (",degree,"C)")),cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5 )
plot(p34, col= scales::alpha('#4a3fa8',0.4),freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
plot(pRS21, col= scales::alpha('darkcyan',0.4),angle = 45,density=25,freq= F,border=scales::alpha('darkcyan',0.6), ylab = "Density",add = T)
plot(pRS22, col= scales::alpha('#4a3fa8',0.4),angle = -45,density=25,freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
#suitability-Tmax
plot(p36, col= scales::alpha('darkcyan',0.4),freq= F,border=scales::alpha('darkcyan',0.6),xlim= c(28,42),ylim = c(0,.5),main = expression(italic(paste("T"[max])-S)),xlab = expression(paste("Temperature (",degree,"C)")),cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5  ,ylab = "Density",cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5 )
plot(p37, col= scales::alpha('#4a3fa8',0.4),freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
plot(pRS23, col= scales::alpha('darkcyan',0.4),angle = 45,density=25,freq= F,border=scales::alpha('darkcyan',0.6), ylab = "Density",add = T)
plot(pRS24, col= scales::alpha('#4a3fa8',0.4),angle = -45,density=25,freq= F,border=scales::alpha('#4a3fa8',0.6), ylab = "Density",add = T)
par(mfrow = c(1,1))
dev.off()
####################################################333
###################################################33
#################################################3
###############################################


summary.fn = function(df){
  med = median(df, na.rm = T)
  int = c(HPDinterval(mcmc(df)))
  out = c(med, int)
  names(out) = c("median", "lowerCI", "upperCI")
  out
}

#Apply summary.fn to each TPdists df
br.RS.dtr9.interval <- apply(br.RS.dtr9.TPdists, MARGIN =2, summary.fn)
br.RS.dtr12.interval <- apply(br.RS.dtr12.TPdists, MARGIN = 2, summary.fn)
br.RS.TThresholds <- list (br.RS.dtr9.interval=br.RS.dtr9.interval,
                           br.RS.dtr12.interval=br.RS.dtr12.interval)
#sink("RS_trait_thresholds_a.txt") #file name
#print(br.RS.TThresholds) #object to be exported as text
#sink() #exporting to source file location

lf.RS.dtr9.interval <- apply(lf.RS.dtr9.TPdists, MARGIN =2, summary.fn)
lf.RS.dtr12.interval <- apply(lf.RS.dtr12.TPdists, MARGIN =2, summary.fn)
lf.RS.TThresholds <- list (lf.RS.dtr9.interval=lf.RS.dtr9.interval,
                           lf.RS.dtr12.interval=lf.RS.dtr12.interval)
#sink("RS_trait_thresholds_lf.txt") #file name
#print(lf.RS.TThresholds) #object to be exported as text
#sink() #exporting to source file location

le.RS.dtr9.interval <- apply(le.RS.dtr9.TPdists, MARGIN =2, summary.fn)
le.RS.dtr12.interval <- apply(le.RS.dtr12.TPdists, MARGIN =2, summary.fn)
le.RS.TThresholds <- list (le.RS.dtr9.interval=le.RS.dtr9.interval,
                           le.RS.dtr12.interval=le.RS.dtr12.interval)
#sink("RS_trait_thresholds_le.txt") #file name
#print(le.RS.TThresholds) #object to be exported as text
#sink() #exporting to source file location

R0.RS.dtr9.direct.interval <- apply(R0.RS.dtr9.TPdists, MARGIN = 2, summary.fn)
R0.RS.dtr12.direct.interval <- apply(R0.RS.dtr12.TPdists, MARGIN = 2, summary.fn)
R0.RS.direct.TThresholds <- list (R0.RS.dtr9.direct.interval=R0.RS.dtr9.direct.interval,
                                  R0.RS.dtr12.direct.interval=R0.RS.dtr12.direct.interval)
#sink("RS_trait_thresholds_S_direct.txt") #file name
#print(R0.RS.direct.TThresholds) #object to be exported as text
#sink() #exporting to source file location

##############################################################
####################################################33
#############################################################
# SUPPLEMENTAL FIGURE x
# SHOWING THE DIFFERENTIAL BETWEEN OBS TRAIT VALUE AND RS ESTIMATE ACROSS TEMP
# T.summary files need to be summary TRUE option
################################################
###########################################33
##Lifespan
#pdf("plots/SIFIG6_differential.pdf", width = 8, height = 12)#Run to make the full figure at once if all of the code is loaded
pdf("plots/lf_RS_differential.pdf", width = 8, height = 8)
par(mar=c(5,6,5,1))
par(mfrow = c(4,2))
plot((lf.dtr9.sims.PostQuants$mean- lf.c.quad_T.summary$mean) ~ lf.dtr9.sims.PostQuants$temps,xlim = c(0,45),ylim =c(-20,20),main = expression(paste("Lifespan (lf)")),ylab = "Relative to constant temperature", xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, col = "darkcyan",type="n")
abline(v= 20.5, col =scales::alpha("black",1), lty =1, lwd=2)
lines((lf.dtr9.sims.PostQuants$mean-lf.c.quad_T.summary$mean) ~ lf.dtr9.sims.PostQuants$temps, lwd =2,lty = 2, col = "darkcyan")
lines(( lf.dtr12.sims.PostQuants$mean-lf.c.quad_T.summary$mean) ~ lf.dtr12.sims.PostQuants$temps, lwd =2,lty=2, col = "#4a3fa8")
lines(( lf.f9.quad_T.summary$mean-lf.c.quad_T.summary$mean) ~ lf.dtr9.sims.PostQuants$temps, lwd =2,lty = 1, col = "darkcyan")
lines(( lf.f12.quad_T.summary$mean-lf.c.quad_T.summary$mean) ~ lf.dtr12.sims.PostQuants$temps, lwd =2,lty=1, col = "#4a3fa8")
lines(rep(0, length( lf.dtr9.sims.PostQuants$temps))~ lf.dtr9.sims.PostQuants$temps, col ="black", lty=3, lwd =1.5)
text(17,-15, expression(italic(paste("T"[opt]))), cex =1.2)
legend("topright", legend = c("data: DTR 9C", "data: DTR 12C","RS: DTR 9C", "RS: DTR 12C"),lwd =2, lty = c(1,1,2,2), col = c("darkcyan", "#4a3fa8"), bty = 'n')
mtext("a", side = 3, at = -5, cex = 1.8, line = 2)
par(mfrow = c(1,1))
dev.off()

#Lifespan RS estimates relative to each respective fluc temp
pdf("plots/lf_RS_diff_to_observed.pdf", width = 8, height = 8)
par(mar=c(5,6,5,1))
par(mfrow = c(2,2))
plot((lf.dtr9.sims.PostQuants$mean- lf.c.quad_T.summary$mean) ~ lf.dtr9.sims.PostQuants$temps,xlim = c(0,45),ylim =c(-20,20),main = expression(paste("Lifespan (lf)")),ylab = "Relative to fluctuating temperature", xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, col = "darkcyan",type="n")
abline(v= 18.4, col ="darkcyan", lty =1, lwd=2)
abline(v= 17.9, col ="#4a3fa8", lty =1, lwd=2)
abline(v= 20.5, col ="black", lty =1, lwd=2)
lines((lf.dtr9.sims.PostQuants$mean- lf.f9.quad_T.summary$mean) ~ lf.dtr9.sims.PostQuants$temps, lwd =2,lty = 2, col = "darkcyan")
lines(( lf.dtr12.sims.PostQuants$mean-lf.f12.quad_T.summary$mean) ~ lf.dtr12.sims.PostQuants$temps, lwd =2,lty=2, col = "#4a3fa8")
lines(rep(0, length( lf.dtr9.sims.PostQuants$temps))~ lf.dtr9.sims.PostQuants$temps, col ="black", lty=3, lwd =1.5)
text(13,-15, expression(italic(paste("T"[opt]))), cex =1.2, col = "black")
#legend("topright", legend = c("RS: DTR 9C", "RS: DTR 12C"),lwd =2, lty = c(2,2), col = c("darkcyan", "#4a3fa8"), bty = 'n')
mtext("b", side = 3, at = -5, cex = 1.8, line = 2)
par(mfrow = c(1,1))
dev.off()


##bite rate
pdf("plots/br_RS_differential.pdf", width = 8, height = 8)
par(mar=c(5,6,5,1))
par(mfrow = c(2,2))
plot(( br.dtr9.sims.PostQuants$mean-br.c.briere_T.summary$mean) ~ br.dtr9.sims.PostQuants$temps,xlim = c(0,45),ylim =c(-.4,.4),main = expression(paste("Bite rate (a)")),ylab = "Relative to constant temperature", xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, col = "darkcyan",type="n")
abline(v= 33.8, col =scales::alpha("black",1), lty =1, lwd=2)
lines(( br.dtr9.sims.PostQuants$mean-br.c.briere_T.summary$mean) ~ br.dtr9.sims.PostQuants$temps, lwd =2,lty = 2, col = "darkcyan")
lines(( br.dtr12.sims.PostQuants$mean-br.c.briere_T.summary$mean) ~ br.dtr12.sims.PostQuants$temps, lwd =2,lty=2, col = "#4a3fa8")
lines(( br.f9.briere_T.summary$mean-br.c.briere_T.summary$mean) ~ br.dtr9.sims.PostQuants$temps, lwd =2,lty = 1, col = "darkcyan")
lines(( br.f12.briere_T.summary$mean- br.c.briere_T.summary$mean) ~ br.dtr12.sims.PostQuants$temps, lwd =2,lty=1, col = "#4a3fa8")
lines(rep(0, length( br.dtr9.sims.PostQuants$temps))~ br.dtr9.sims.PostQuants$temps, col ="black", lty=3, lwd =1.5)
#legend("topleft", legend = c("data: DTR 9C", "data: DTR 12C","RS: DTR 9C", "RS: DTR 12C"),lwd =2, lty = c(1,1,2,2), col = c("darkcyan", "#4a3fa8"), bty = 'n')
text(38,0.3, expression(italic(paste("T"[opt]))), cex =1.2)
mtext("c", side = 3, at = -5, cex = 1.8, line = 2)
par(mfrow = c(1,1))
dev.off()

#Lifespan RS estimates relative to each respective fluc temp
pdf("plots/br_RS_diff_to_observed.pdf", width = 8, height = 8)
par(mar=c(5,6,5,1))
par(mfrow = c(2,2))
plot(( br.dtr9.sims.PostQuants$mean-br.c.briere_T.summary$mean) ~ br.dtr9.sims.PostQuants$temps,xlim = c(0,45),ylim =c(-.4,.4),main = expression(paste("Bite rate (a)")),ylab = "Relative to fluctuating temperature", xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, col = "darkcyan",type="n")
abline(v= 29.7, col ="darkcyan", lty =1, lwd=2)
abline(v= 31.3, col ="#4a3fa8", lty =1, lwd=2)
abline(v= 33.8, col ="black", lty =1, lwd=2)
lines((br.dtr9.sims.PostQuants$mean- br.f9.briere_T.summary$mean) ~ br.dtr9.sims.PostQuants$temps, lwd =2,lty = 2, col = "darkcyan")
lines(( br.dtr12.sims.PostQuants$mean- br.f12.briere_T.summary$mean) ~ br.dtr12.sims.PostQuants$temps, lwd =2,lty=2, col = "#4a3fa8")
lines(rep(0, length( br.dtr9.sims.PostQuants$temps))~ br.dtr9.sims.PostQuants$temps, col ="black", lty=3, lwd =1.5)
text(25,0.3, expression(italic(paste("T"[opt]))), cex =1.2)
mtext("d", side = 3, at = -5, cex = 1.8, line = 2)
par(mfrow = c(1,1))
dev.off()

#Lifetime eggs
pdf("plots/le_RS_differential.pdf", width = 8, height = 8)
par(mar=c(5,6,5,1))
par(mfrow = c(2,2))
plot(( le.dtr9.sims.PostQuants$mean- le.c.briere_T.summary$mean) ~ le.dtr9.sims.PostQuants$temps,xlim = c(0,45),ylim =c(-220,220),main = expression(paste("Lifetime eggs (B)")),ylab = "Relative to constant temperature", xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, col = "darkcyan",type="n")
abline(v= 27.6, col =scales::alpha("black",1), lty =1, lwd=2)
lines(( le.dtr9.sims.PostQuants$mean-le.c.briere_T.summary$mean) ~ le.dtr9.sims.PostQuants$temps, lwd =2,lty = 2, col = "darkcyan")
lines(( le.dtr12.sims.PostQuants$mean-le.c.briere_T.summary$mean) ~ le.dtr12.sims.PostQuants$temps, lwd =2,lty=2, col = "#4a3fa8")
lines(( le.f9.briere_T.summary$mean-le.c.briere_T.summary$mean) ~ le.dtr9.sims.PostQuants$temps, lwd =2,lty = 1, col = "darkcyan")
lines(( le.f12.briere_T.summary$mean- le.c.briere_T.summary$mean) ~ le.dtr12.sims.PostQuants$temps, lwd =2,lty=1, col = "#4a3fa8")
lines(rep(0, length( le.dtr9.sims.PostQuants$temps))~ le.dtr9.sims.PostQuants$temps, col ="black", lty=3, lwd =1.5)
#legend("topleft", legend = c("data: DTR 9C", "data: DTR 12C","RS: DTR 9C", "RS: DTR 12C"),lwd =2, lty = c(1,1,2,2), col = c("darkcyan", "#4a3fa8"), bty = 'n')
text(23,170, expression(italic(paste("T"[opt]))), cex =1.2)
mtext("e", side = 3, at = -5, cex = 1.8, line = 2)
par(mfrow = c(1,1))
dev.off()

##Lifetime eggs -- fluctuating temperatures
pdf("plots/le_RS_diff_to_observed.pdf", width = 8, height = 8)
par(mar=c(5,6,5,1))
par(mfrow = c(2,2))
plot(( le.dtr9.sims.PostQuants$mean- le.c.briere_T.summary$mean) ~ le.dtr9.sims.PostQuants$temps,xlim = c(0,45),ylim =c(-220,220),main = expression(paste("Lifetime eggs (B)")),ylab = "Relative to fluctuating temperature", xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, col = "darkcyan",type="n")
abline(v= 26.4, col ="darkcyan", lty =1, lwd=2)
abline(v= 24.7, col ="#4a3fa8", lty =1, lwd=2)
abline(v= 27.6, col ="black", lty =1, lwd=2)
lines((le.dtr9.sims.PostQuants$mean- le.f9.briere_T.summary$mean) ~ le.dtr9.sims.PostQuants$temps, lwd =2,lty = 2, col = "darkcyan")
lines(( le.dtr12.sims.PostQuants$mean- le.f12.briere_T.summary$mean) ~ le.dtr12.sims.PostQuants$temps, lwd =2,lty=2, col = "#4a3fa8")
lines(rep(0, length( le.dtr9.sims.PostQuants$temps))~ le.dtr9.sims.PostQuants$temps, col ="black", lty=3, lwd =1.5)
text(20,170, expression(italic(paste("T"[opt]))), cex =1.2)
mtext("f", side = 3, at = -5, cex = 1.8, line = 2)
par(mfrow = c(1,1))
dev.off()

#Thermal Suitability
pdf("plots/S_RS_differential.pdf", width = 8, height = 8)
par(mar=c(5,6,5,1))
par(mfrow = c(2,2))
plot(( R0.RS.dtr9.out$mean/max(R0.RS.dtr9.out$mean)- R0.c.uni.out$mean/max(R0.c.uni.out$mean)) ~ Temp.xs,xlim = c(0,45),ylim =c(-0.7,0.7),main = expression(paste("Suitability (S)")),ylab = "Relative to constant temperature", xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, col = "darkcyan",type="n")
abline(v= 26.8, col =scales::alpha("black",1), lty =1, lwd=2)
lines(( R0.RS.dtr9.out$mean/max(R0.RS.dtr9.out$mean)- R0.c.uni.out$mean/max(R0.c.uni.out$mean)) ~ Temp.xs, lwd =2,lty = 2, col = "darkcyan")
lines(( R0.RS.dtr12.out$mean/max(R0.RS.dtr12.out$mean)- R0.c.uni.out$mean/max(R0.c.uni.out$mean)) ~ Temp.xs, lwd =2,lty = 2, col = "#4a3fa8")
lines(( R0.dtr9.uni.out$mean/max(R0.dtr9.uni.out$mean)- R0.c.uni.out$mean/max(R0.c.uni.out$mean)) ~ Temp.xs, lwd =2,lty = 1, col = "darkcyan")
lines(( R0.dtr12.uni.out$mean/max(R0.dtr12.uni.out$mean)- R0.c.uni.out$mean/max(R0.c.uni.out$mean)) ~ Temp.xs, lwd =2,lty = 1, col = "#4a3fa8")
lines(rep(0, length( R0.RS.dtr9.out$mean))~Temp.xs , col ="black", lty=3, lwd =1.5)
#legend("topleft", legend = c("data: DTR 9C", "data: DTR 12C","RS: DTR 9C", "RS: DTR 12C"),lwd =2, lty = c(1,1,2,2), col = c("darkcyan", "#4a3fa8"), bty = 'n')
text(31,0.6, expression(italic(paste("T"[opt]))), cex =1.2)
mtext("g", side = 3, at = -5, cex = 1.8, line = 2)
par(mfrow = c(1,1))
dev.off()

###fluctuating temps- suitability
pdf("plots/S_RS_diff_to_observed.pdf", width = 8, height = 8)
par(mar=c(5,6,5,1))
par(mfrow = c(2,2))
plot(( le.dtr9.sims.PostQuants$mean- le.c.briere_T.summary$mean) ~ le.dtr9.sims.PostQuants$temps,xlim = c(0,45),ylim =c(-0.7,0.7),main = expression(paste("Suitability (S)")),ylab = "Relative to fluctuating temperature", xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.5, col = "darkcyan",type="n")
abline(v= 25.8, col ="darkcyan", lty =1, lwd=2)
abline(v= 24.3, col ="#4a3fa8", lty =1, lwd=2)
abline(v= 26.8, col ="black", lty =1, lwd=2)
lines(( R0.RS.dtr9.out$mean/max(R0.RS.dtr9.out$mean)- R0.dtr9.uni.out$mean/max(R0.dtr9.uni.out$mean)) ~ Temp.xs, lwd =2,lty = 2, col = "darkcyan")
lines(( R0.RS.dtr12.out$mean/max(R0.RS.dtr12.out$mean)- R0.dtr12.uni.out$mean/max(R0.dtr12.uni.out$mean)) ~ Temp.xs, lwd =2,lty = 2, col = "#4a3fa8")
lines(rep(0, length( R0.RS.dtr9.out$mean))~ Temp.xs, col ="black", lty=3, lwd =1.5)
text(20,0.6, expression(italic(paste("T"[opt]))), cex =1.2)
mtext("h", side = 3, at = -5, cex = 1.8, line = 2)
par(mfrow = c(1,1))
dev.off()


###Use the code of the observed data for DTR 9 and 12 and then add the RS estimates on top of for both trait and direct S on.
####Main Figure 3a-  
# Plot thermal performance curves for observed data over complete RS estimates
pdf("plots/Main_Fig3.pdf", width = 8, height = 12)
par(mfrow = c(3,2))
plot.R0.poly3(Temp.xs, R0.c.uni.out, cols = "black", add = "false")
#median lines from RS on each trait
lines(Temp.xs, R0.RS.dtr9.out$median/max(R0.RS.dtr9.out$median), col='darkcyan', lty = 2, lwd = 2)
lines(Temp.xs, R0.RS.dtr12.out$median/max(R0.RS.dtr12.out$median), col='#4a3fa8', lty = 2, lwd = 2)

#mean line from RS directly on constant S
lines(R0.dtr9.sims.PostQuants$temps, R0.dtr9.sims.PostQuants$median/max(R0.dtr9.sims.PostQuants$median), col='darkcyan', lty = 3, lwd = 2)
lines(R0.dtr12.sims.PostQuants$temps, R0.dtr12.sims.PostQuants$median/max(R0.dtr12.sims.PostQuants$median), col='#4a3fa8', lty = 3, lwd = 2)
mtext("a", side = 3, at = -5, cex = 1.8, line = 2)
legend("topleft", legend = c("data: constant","RS (traits): DTR 9C", "RS (traits): DTR 12C","RS (S): DTR 9C","RS (S): DTR 12C"),lwd =2, lty = c(1,2,2,3,3), col = c("black","darkcyan", "#4a3fa8"), bty = 'n')

par(mfrow = c(1,1))
dev.off()