#Set working directory
mainDir = "C:/Users/Kerri/Desktop/Chapter2 InProgress"
setwd(mainDir)

# Load libraties for plotting traits
library(plyr) # Slices and dices data
library(plotrix) # For standard error function
library(coda)
library(dplyr)
library(tidyverse)

# load("saved posteriors/temps.Rdata")

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


#Plot  all  the temperature programs into a single 4-panel figurefor the Supplemental Info
library(ggplot2)

#panel A - constant temps
SI_Fig1a <- ggplot(temp.programs.2, aes(x= hr, y = X28.dtr12))+
  scale_x_continuous("Time (hr)", limits = c(0,24), breaks = seq(0,24,4))+
  scale_y_continuous( expression(paste("Temperature (",degree,"C)")), limits = c(10,41), breaks = seq(10,40,5))+
  geom_line(aes(y=X28.dtr12), colour = "NA")+ #tricks the graph to place hours on it
  geom_hline(aes(yintercept=16, colour = "16C"),size = 1.2)+
  geom_hline(aes(yintercept=20, colour = "20C"), size = 1.2)+
  geom_hline(aes(yintercept=24, colour = "24C"),size = 1.2)+
  geom_hline(aes(yintercept=28, colour = "28C"),size = 1.2)+
  geom_hline(aes(yintercept=32, colour = "32C"),size = 1.2)+
  geom_hline(aes(yintercept=36, colour = "36C"), size = 1.2)+
  scale_color_manual("", values = c("purple","royalblue","darkgreen","orange","red","black"))+
  theme_bw(base_size = 18)+
  guides(colour=guide_legend(nrow=1))+
  theme(legend.margin =margin(t = 0, unit='cm'),legend.direction = "horizontal", legend.position = c(0.5,.94), legend.title = element_text(size=9), legend.text = element_text(size = 9))


#panel B - dtr of 9
SI_Fig1b <- ggplot(temp.programs.2, aes(x= hr))+
  scale_x_continuous("Time (hr)", limits = c(0,24), breaks = seq(0,24,4))+
  scale_y_continuous( expression(paste("Temperature (",degree,"C)")), limits = c(10,41), breaks = seq(10,40,5))+
  geom_line(aes(y=X16.dtr9, colour = "16C"), linetype = "dashed", size = 1.2)+
  geom_line(aes(y=X20.dtr9, colour = "20C"), linetype = "dashed", size = 1.2)+
  geom_line(aes(y=X24.dtr9, colour = "24C"), linetype = "dashed", size = 1.2)+
  geom_line(aes(y=X28.dtr9, colour = "28C"), linetype = "dashed", size =1.2)+
  geom_line(aes(y=X32.dtr9, colour = "32C"), linetype = "dashed", size =1.2)+
  scale_color_manual("", values = c("purple","royalblue","darkgreen","orange","red","black"))+
  theme_bw(base_size = 18)+
  guides(colour=guide_legend(nrow=1))+
  theme(legend.margin =margin(t = 0, unit='cm'),legend.direction = "horizontal", legend.position = c(0.5,.94), legend.title = element_text(size=9), legend.text = element_text(size = 9))


#panel C - dtr of 12
SI_Fig1c <- ggplot(temp.programs.2, aes(x= hr))+
  scale_x_continuous("Time (hr)", limits = c(0,24), breaks = seq(0,24,4))+
  scale_y_continuous( expression(paste("Temperature (",degree,"C)")), limits = c(10,41), breaks = seq(10,40,5))+
  geom_line(aes(y=X16.dtr12, colour = "16C"), linetype = "dotted", size = 1.2)+
  geom_line(aes(y=X20.dtr12, colour = "20C"), linetype = "dotted", size = 1.2)+
  geom_line(aes(y=X24.dtr12, colour = "24C"), linetype = "dotted", size = 1.2)+
  geom_line(aes(y=X28.dtr12, colour = "28C"), linetype = "dotted", size = 1.2)+
  geom_line(aes(y=X32.dtr12, colour = "32C"), linetype = "dotted", size = 1.2)+
  scale_color_manual("", values = c("purple","royalblue","darkgreen","orange","red","black"))+
  theme_bw(base_size = 18)+
  guides(colour=guide_legend(nrow=1))+
  theme(legend.margin =margin(t = 0, unit='cm'),legend.direction = "horizontal", legend.position = c(0.5,.94), legend.title = element_text(size=9), legend.text = element_text(size = 9))


#panel D - overlay of all treatments
SI_Fig1d <- ggplot(temp.programs.2, aes(x= hr))+
  scale_x_continuous("Time (hr)", limits = c(0,24), breaks = seq(0,24,4))+
  scale_y_continuous( expression(paste("Temperature (",degree,"C)")), limits = c(10,41), breaks = seq(10,40,5))+
  geom_hline(aes(yintercept=16, colour = "16C", linetype = "constant"),size = 1.2)+
  geom_hline(aes(yintercept=20, colour = "20C", linetype = "constant"), size = 1.2)+
  geom_hline(aes(yintercept=24, colour = "24C", linetype = "constant"),size = 1.2)+
  geom_hline(aes(yintercept=28, colour = "28C", linetype = "constant"),size = 1.2)+
  geom_hline(aes(yintercept=32, colour = "32C", linetype = "constant"),size = 1.2)+
  geom_hline(aes(yintercept=36, colour = "36C", linetype = "constant"), size = 1.2)+
  geom_line(aes(y=X16.dtr9, colour = "16C", linetype = "DTR 09C"), size = 1.2)+
  geom_line(aes(y=X20.dtr9, colour = "20C", linetype = "DTR 09C"), size = 1.2)+
  geom_line(aes(y=X24.dtr9, colour = "24C", linetype = "DTR 09C"), size = 1.2)+
  geom_line(aes(y=X28.dtr9, colour = "28C", linetype = "DTR 09C"), size =1.2)+
  geom_line(aes(y=X32.dtr9, colour = "32C", linetype = "DTR 09C"), size =1.2)+
  geom_line(aes(y=X16.dtr12, colour = "16C", linetype = "DTR 12C"), size = 1.2)+
  geom_line(aes(y=X20.dtr12, colour = "20C", linetype = "DTR 12C"), size = 1.2)+
  geom_line(aes(y=X24.dtr12, colour = "24C", linetype = "DTR 12C"), size = 1.2)+
  geom_line(aes(y=X28.dtr12, colour = "28C", linetype = "DTR 12C"), size = 1.2)+
  geom_line(aes(y=X32.dtr12, colour = "32C", linetype = "DTR 12C"), size = 1.2)+
  scale_color_manual("", values = c("purple","royalblue","darkgreen","orange","red","black"))+
  scale_linetype_manual("", values = c("solid", "dashed","dotted"), drop = TRUE)+
  theme_bw(base_size = 18)+
  guides(linetype=guide_legend(nrow=1))+
  guides(colour = FALSE)+
  theme(legend.spacing.y = unit(-0.16, "cm"),legend.margin =margin(t = 0, unit='cm'),legend.direction = "horizontal", legend.position = c(0.5,.94), legend.title = element_text(size=9), legend.text = element_text(size = 9), legend.key = element_rect(colour = "transparent", fill = "white"))

library(cowplot)
library(ggpubr)

supp1 <- plot_grid(SI_Fig1a, SI_Fig1b, SI_Fig1c, SI_Fig1d, labels=c("a", "b", "c", "d"), ncol = 2, nrow = 2, align = "hv", label_size = 20)
supp1
save_plot("plots/supp1_tempprofiles.png", supp1, ncol = 2, nrow = 2, base_aspect_ratio = 1.4)


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
# For EEID Poster

# Export just the columns we need
temp.programs.subset <- data.frame(Time = seq(0, 23, 1),
								   X16.cons = rep(16, times = 24),
								   X20.cons = rep(20, times = 24),
								   X24.cons = rep(24, times = 24),
								   X28.cons = rep(28, times = 24),
								   X32.cons = rep(32, times = 24),
								   X36.cons = rep(36, times = 24),
								   X16.dtr9 = temp.programs.2$X16.dtr9,
								   X20.dtr9 = temp.programs.2$X20.dtr9,
								   X24.dtr9 = temp.programs.2$X24.dtr9,
								   X28.dtr9 = temp.programs.2$X28.dtr9,
								   X32.dtr9 = temp.programs.2$X32.dtr9,
								   X16.dtr12 = temp.programs.2$X16.dtr12,
								   X20.dtr12 = temp.programs.2$X20.dtr12,
								   X24.dtr12 = temp.programs.2$X24.dtr12,
								   X28.dtr12 = temp.programs.2$X28.dtr12,
								   X32.dtr12 = temp.programs.2$X32.dtr12)


# Change into long format
temp.programs.long <- temp.programs.subset %>%
	pivot_longer(!Time, names_to = "Treatment", values_to = "Temp") %>% 
	separate_wider_delim(Treatment, delim=".", names = c("Meantemp", "DTR")) %>%
	mutate(Meantemp=as.factor(sub("X", "", Meantemp))) %>% 
	mutate(DTR=as.factor(DTR))

# Reverse the order of of mean temp factor levels so they display with hottest on top
temp.programs.long$Meantemp <- factor(temp.programs.long$Meantemp, levels = rev(levels(temp.programs.long$Meantemp)))

# Make plot for EEID
EEID_Fig <- ggplot(temp.programs.long, aes(x = Time, y = Temp, colour = Meantemp, linetype = DTR)) +
	scale_x_continuous("Time (hour of day)", limits = c(0,24), breaks = seq(0,24,4)) +
	scale_y_continuous( expression(paste("Temperature (",degree,"C)")), limits = c(10,41), breaks = seq(10,40,5)) +
	geom_line(linewidth = 1.2, alpha = 0.8) +
	#scale_color_manual("Mean", values = c("#4A54B7","#944CB1","#C8419C","#EC4175","#FE565D","#FF7639")) +
	scale_color_manual("Mean", values = c("#e3ad1d", "#FF7639","#EC4175", "#944CB1", "#4A54B7", "#175377")) +
	#scale_color_viridis_d("Mean", option = "rocket") +
	scale_linetype_manual("DTR", values = c("solid", "twodash","dotted"), labels = c("0", "9", "12")) +
	theme_bw(base_size = 18) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(EEID_Fig)

# Make plot for publication
Pub_Fig <- ggplot(temp.programs.long, aes(x = Time, y = Temp, colour = Meantemp, linetype = DTR)) +
	scale_x_continuous("Time (hour of day)", limits = c(0,24), breaks = seq(0,24,4)) +
	scale_y_continuous( expression(paste("Temperature (",degree,"C)")), limits = c(10,41), breaks = seq(10,40,5)) +
	geom_line(linewidth = 1.2, alpha = 0.8) +
	#scale_color_manual("Mean", values = c("#4A54B7","#944CB1","#C8419C","#EC4175","#FE565D","#FF7639")) +
	scale_color_manual("Mean", values = c("#e3ad1d", "#FF7639","#EC4175", "#944CB1", "#4A54B7", "#175377")) +
	#scale_color_viridis_d("Mean", option = "rocket") +
	scale_linetype_manual("DTR", values = c("solid", "twodash","dotted"), labels = c("0", "9", "12")) +
	theme_bw(base_size = 18) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
	theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
	theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10))
print(Pub_Fig)

ggsave('figures/EEID_tempprofiles.pdf', EEID_Fig, width = 8.5, height = 4)
ggsave('figures/Publication_tempprofiles.pdf', EEID_Fig, width = 7, height = 4)
ggsave('figures/Publication_tempprofiles.pdf', Pub_Fig, width = 6, height = 3)
#save_plot("figures/EEID_tempprofiles.png", EEID_Fig, base_aspect_ratio = 2.2)
