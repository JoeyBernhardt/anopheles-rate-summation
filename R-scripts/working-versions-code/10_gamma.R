
##### Calculating gamma!

# # Further, the S(T) formulation uses the Gompertz function over daily adult 
# survival and the extrinsic incubation period (EIP, the inverse of the parasite
# development rate (PDR-1)) to calculate the proportion of mosquitoes surviving 
# the latency period (ϒ) as described in [19].

# ## To estimate ϒ, we fit a Gompertz distribution to survivorship
# data from each temperature and replicate. We then took the
# proportion of mosquitoes alive upon completion of the predicted
# extrinsic incubation period (PDR50(T)−1) of P. falciparum at each
# temperature. The amount of days to reach 50% of maximum
# infectiousness in a mosquito population is represented by
# PDR50(T)−1 [29]. This formulation allows us to account for agedependent
# mortality in the proportion of mosquitoes surviving
# the latency period (ϒ). We then compared the thermal responses
# of lifespan, biting rate and lifetime egg production for
# An. stephensi when these traits are directly observed (lf, a, B)
# versus estimated (lf*, a*, B*) from the data generated in this
# study, as well as if any observed differences alter the predicted
# thermal suitability for malaria transmission (R0).

data.constant <- read.csv("data-raw/constant.individual.trait.csv")
