#Test to see the difference  using the Bayesian model mean vs. median values on the visual of 
# R0 for a dtr of 9,12 with the emprical data.

#Used to determine if this is the reason for the discrepancy btw the TPdists value and the Figure in Figure 2.


#R0 equation
# constant to keep lifespan from being numerically zero
# assume minimum survival time is half an hour
ec.lf = 1/48

# constant to keep PDR from being numerically zero
# assume the maximum EIP is 100 days
ec.pdr = 1/100

# 2. M depends on lifetime fecundity; gamma TPC is substituted for exp^(-1/(lf*PDR))expression
R0.full = function(a, bc, lf, y, B, pEA, MDR){
  M = B * pEA * MDR * (lf+ec.lf)
  (a^2 * bc * y * M * (lf+ec.lf) )^0.5
}

#3. dtr 12 with mean values
#Define the function for each equation
  #Quadratic function; -c(T-Tmin)(T-Tmax)
  #Briere function; cT(T-Tmin)(Tmax-T)1/2
lf.q.mean <- function(x){-0.3191236  *(x-5.5729360)*(x-30.6481572)}
le.b.mean <- function(x){0.4801446 *x*(x- 12.8019914)*(29.0380737-x)^(1/2)}
br.b.mean <- function(x){1.514858e-04 *x*(x-1.468721e+00 )*(3.876619e+01-x)^(1/2)}
g.q.mean <- function(x){-0.002390491*(x-4.193786945)*(x- 42.059015808 )}
bc.q.mean <- function(x){-0.002747507* (x-8.767342129 )*(x-39.709208937)}
pea.q.mean <- function(x){-0.007996182*(x-  15.356924052)*(x- 37.078852667 )}
mdr.b.mean <- function(x){1.065409e-04*x*(x-1.335163e+01)*(3.596750e+01-x)^(1/2)}

#dtr 12 with median values
lf.q.med <- function(x){-0.266742*(x- 4.880488)*(x-30.618093)}
le.b.med <- function(x){0.4623845 *x*(x-13.0905227 )*(28.5405581 -x)^(1/2)}
br.b.med <- function(x){ 1.480739e-04*x*(x-1.041452e+00)*(3.867934e+01-x)^(1/2)}
g.q.med <- function(x){- 0.002238986 *(x-4.057460369 )*(x-42.251567958 )}
bc.q.med <- function(x){-0.002445606*(x-9.147378065)*(x-39.563056169)}
pea.q.med <- function(x){- 0.007984157 *(x-15.370718932)*(x-37.073399631)}
mdr.b.med <- function(x){1.064575e-04*x*(x-1.336170e+01)*(3.597729e+01-x)^(1/2)}

#plot together
R0.mean.values <- R0.full(br.b.mean(seq(0,50,0.1)),bc.q.mean(seq(0,50,0.1)),lf.q.mean(seq(0,50,0.1)),g.q.mean(seq(0,50,0.1)),le.b.mean(seq(0,50,0.1)),pea.q.mean(seq(0,50,0.1)),mdr.b.mean(seq(0,50,0.1)))
R0.median.values <- R0.full(br.b.med(seq(0,50,0.1)),bc.q.med(seq(0,50,0.1)),lf.q.med(seq(0,50,0.1)),g.q.med(seq(0,50,0.1)),le.b.med(seq(0,50,0.1)),pea.q.med(seq(0,50,0.1)),mdr.b.med(seq(0,50,0.1)))
#Replace NA with 0 for finding Tmin
R0.mean.values[is.na(R0.mean.values)]<-0
R0.median.values[is.na(R0.median.values)] <- 0
temps <- seq(0,50,0.1)
plot(R0.mean.values/max(R0.mean.values)~temps,na.rm = T, type = "n", ylab = "Suitability", xlab = "Temperature")
lines(R0.mean.values/max(R0.mean.values)~temps, col="black", lwd =2)
lines (R0.median.values/max(R0.median.values)~temps, col = "seagreen3", lwd =2)

#look at overall Tmin values
#step1: extract the index of the peak
index.peak <- which.max(R0.mean.values/max(R0.mean.values))
#step2: find index of last 0 temp before the peak
zero.indexies <- which(R0.mean.values == 0)
#first element that is above index.peak
above.index <- which(zero.indexies < index.peak)

#last element of that last should correspond to the Tmin
temps[above.index[length(above.index)]] #10.2

############Now the median list##########
#look at overall Tmin values
#step1: extract the index of the peak
index.peak <- which.max(R0.median.values/max(R0.median.values))
#step2: find index of last 0 temp before the peak
zero.indexies <- which(R0.median.values == 0)
#first element that is above index.peak
above.index <- which(zero.indexies < index.peak)

#last element of that last should correspond to the Tmin
temps[above.index[length(above.index)]]  #9.9


########Conclusion, either the R0 estimates are flawed or it is a rounding issue. look at the TP thresholds again


####Looking at why there is a performance boost at the thermal extremes for the quadratic fits.
#I suspect it is because we limit our functions to 0 for any negative values, which then yields predictions counter to
#what we would expect based on rate summation

#Use pea function as an example- I think the boost occurs where if the function went negative it is equal to the postive y-range portion
x.vals <- seq(0,50,0.1)
#round to circumvent matching errors
x.vals <- round(x.vals, 4)

pea.postive <- pea.q.mean(x.vals) #see index at where it starts to become negative (should be Tmin)
pea.postive[154]#last negative; 155 first positive
x.vals[154] #15.3C
pea.postive[371]#last postive; 372 first negative
x.vals[372] #37.1C
plot(x.vals, abs(pea.postive)) #good, find the max value in the middle
plot(x.vals, pea.postive)
#Then find the x value outside of the middle that matches on either side

#find max between 154 and 372
plot(seq(15.3,37.1,0.1),pea.q.mean(seq(15.3,37.1,0.1))) 
function.max <- max(pea.q.mean(seq(15.3,37.1,0.1)))
function.max #0.9432
#Now find lower and upper that matches
Linflection.index <- abs(pea.q.mean(seq(0,15.2,0.1)))
Linflection.index[109]
Linflection.index[110]
x.vals[109] #10.8
x.vals[110] #10.9

check.dtr9 <-15.3-4.5 #10.8
check.dtr12 <- 15.3-6 #9.3
#plot the RS graph over a smaller x-axis to see if the above numbers correspond

#Ok this represents where RS values start to be postive, 
#but not the point RS estimates go from being overpredictive to underpredictive....

#Following the same logic that position should be where temp fluctation does not cross the 0 point
inflection.dtr9 <- 15.3+4.5 #19.8
inflection.dtr12 <- 15.3+6 #21.6
##Not correct either....correct value is right around 16C....

###Look at the concavity of br without being restricted to 0
br.postive <- br.b.mean(seq(-20,50,0.1))
le.postive<-le.b.mean(seq(-20,50,0.1))
mdr.postive<-mdr.b.mean(seq(-20,50,0.1))
plot(seq(-20,50,0.1), br.postive) #convection increases;;;; would have more of a rescue effect
plot(seq(-20,50,0.1), le.postive)
plot(seq(-20,50,0.1), mdr.postive)
#but we don't see a rescue effect based on the experimental data with temp fluc
