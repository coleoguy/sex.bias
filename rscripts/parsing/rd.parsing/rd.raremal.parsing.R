# Julio Rincones-Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com
# P.I. Dr. Heath Blackmon
# coleoguy@gmail.com

# This script has been written to parse through the additional simulations carried out to determine
# at which rd value there is a "flip" in net feminization under the original conditions used for 
# our simulations of a XY sex determination system where the fate of sexually antagonistic
# variation is modeled using finite populations, biases in operational sex ratio, different selection
# coefficients (s = 0.1, 0.2, 0.5, and 0.9), genetic architectures (h = 0, 0.5, 1, and .99), and
# recombination distances.


# load the functions used to calculate the relative fitness of each sex.
source("../sex.bias/rscripts/functions/analysis.functions.R")

# load the results for instances where males are the rare sex and store to a vector. 
load("../sex.bias/rscripts/results/rd.rare.males.RData")
res.rare.male <- results

# remove the surplus "results" vector.
rm(results)

# For the sake of clarity, the RData structure takes the following form:
# res.rare.male$females500$males0.2$rd0$h0$s0.1


# Provide the formals for the functions in the loop used to parse the data
loc <- "sex"
h <- c(0, .5, 1, .99)
s <- .5

# save results to data structure in order to facilitate calculation of male and female relative
# fitness. We have 26 rd values, the first 25 are for sex chromosomes and so, this might require 2
# going through this script as many times as required to obtain results for the the desired number of
# rd values

# creates empty matrices where the results from the loop on line 72 will be stored.
# h = 0
results.fit.h0 <- matrix(, 3, 2)

colnames(results.fit.h0) <- c("500", "100")

rownames(results.fit.h0) <- c("OSR.2.h0", "OSR.1.h0", "OSR.05.h0")

# h = 0.5

results.fit.h0.5 <- matrix(, 3, 2)

colnames(results.fit.h0.5) <- c("500", "100")

rownames(results.fit.h0.5) <- c("OSR.2.h.5", "OSR.1.h.5", "OSR.05.h.5")

# h = 1

results.fit.h1 <- matrix(, 3, 2)

colnames(results.fit.h1) <- c("500", "100")

rownames(results.fit.h1) <- c("OSR.2.h1", "OSR.1.h1", "OSR.05.h1")

# h =  .99

results.fit.h99 <- matrix(, 3, 2)

colnames(results.fit.h99) <- c("500", "100")

rownames(results.fit.h99) <- c("OSR.2.h99", "OSR.1.h99", "OSR.05.h99")

# Cycle through the common sex population size (500, 100)
for(i in 1:2){
  # cycle through OSR bias in this order: 0.2, 0.1, 0.05
  for(j in 1:3){
    # Values of rd will need to be changed manually for every instance of interest
    # rd values are: 0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26,
    # 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, .5
    # Omit .5 as this is equivalent to autosomal loci.
    for(k in 1:25){
      if(names(res.rare.male[[i]][[j]])[k] == "rd0.26"){
        # cycle through genetic architectures. The selection coefficient is hard coded. If you want to choose a different
        # coefficient you need to change '$s0.5' to any of the values we have used for the simulations(0.1, 0.2, 0.5, 0.9).
        for(m in 1:4){
          if(m == 1){
            for(n in 1:4){
              fit <- getFitness(data = res.rare.male[[i]][[j]][[k]][[m]]$s0.5, loc, h[m], s)
              results.fit.h0[j, i] <- getDifference.fem(fit) 
            }
          }
          if(m == 2){
            for(n in 1:4){
              fit <- getFitness(data = res.rare.male[[i]][[j]][[k]][[m]]$s0.5, loc, h[m], s)
              results.fit.h0.5[j, i] <- getDifference.fem(fit) 
            }
          }
          if(m == 3){
            for(n in 1:4){
              fit <- getFitness(data = res.rare.male[[i]][[j]][[k]][[m]]$s0.5, loc, h[m], s)
              results.fit.h1[j, i] <- getDifference.fem(fit) 
            }
          }
          if(m == 4){
            for(n in 1:4){
              fit <- getFitness(data = res.rare.male[[i]][[j]][[k]][[m]]$s0.5, loc, h[m], s)
              results.fit.h99[j, i] <- getDifference.fem(fit) 
            }
          }
        }
      }
    }
  }
}

# write .csv files to store your results permanently if desired. Remember to change the file name
# to avoid overwriting previously saved results.

write.csv(results.fit.h0, file ="rd0.26.h0.rmal.sexchrom.csv")
write.csv(results.fit.h0.5, file ="rd0.26.h0.5.rmal.sexchrom.csv")
write.csv(results.fit.h1, file ="rd0.26.h1.rmal.sexchrom.csv")
write.csv(results.fit.h99, file ="rd0.26.h.99.rmal.sexchrom.csv")


