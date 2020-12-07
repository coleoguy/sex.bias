# Julio Rincones-Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com
# P.I. Dr. Heath Blackmon
# coleoguy@gmail.com

# This script has been written to parse through the additional simulations carried out to determine
# at which rd value there is a "flip" in net masculinization under the original conditions used for 
# our simulations of a XY sex determination system where the fate of sexually antagonistic
# variation is modeled using finite populations, biases in operational sex ratio, different selection
# coefficients (s = 0.1, 0.2, 0.5, and 0.9), genetic architectures (h = 0, 0.5, 1, and .99), and
# recombination distances.


# load the functions used to calculate the relative fitness of each sex.
source("../sex.bias/rscripts/functions/analysis.functions.R")

# load the results for instances where males are the rare sex and store to a vector. 
load("../sex.bias/rscripts/results/osr1.rd.model.RData")
res.nobias <- results

# remove the surplus data structures
rm(results)
rm(run.result)
rm(x)

# data structure is as follows: res.nobias$pop100$rd0.02$h0$s0.1

# Provide the formals for the functions in the loop used to parse the data
loc <- "sex"
h <- c(0, .5, 1, .99)
s <- .5

# creates empty matrices where the results from the loop on line 72 will be stored.
# h = 0
results.fit.h0 <- matrix(, length(res.nobias$pop100), 2) 
# We have two populations but both have the same rd values, hence why we are using only one
# population to fill the row names.

colnames(results.fit.h0) <- c("100", "500")

rownames(results.fit.h0) <- names(res.nobias$pop100)

# h = 0.5

results.fit.h0.5 <- matrix(, length(res.nobias$pop100), 2)

colnames(results.fit.h0.5) <- c("100", "500")

rownames(results.fit.h0.5) <-  names(res.nobias$pop100)

# h = 1

results.fit.h1 <- matrix(, length(res.nobias$pop100), 2)

colnames(results.fit.h1) <- c("100", "500")

rownames(results.fit.h1) <-  names(res.nobias$pop100)

# h =  .99

results.fit.h99 <- matrix(, length(res.nobias$pop100), 2)

colnames(results.fit.h99) <- c("100", "500")

rownames(results.fit.h99) <-  names(res.nobias$pop100)

# Cycle through the population size of one sex (100, 500)
for(i in 1:2){
  # Cycle through rd values. These values are:
  # 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26,
  # 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48
  for(j in 1:24){
        for(k in 1:4){
          if(k == 1){
            for(n in 1:4){
              fit <- getFitness(data = res.nobias[[i]][[j]][[k]]$s0.5, loc, h[k], s)
              results.fit.h0[j, i] <- getDifference.mal(fit) 
              # either 'getDifference' function should provide the same results in the absence of
              # OSR bias.
            }
          }
          if(k == 2){
            for(n in 1:4){
              fit <- getFitness(data = res.nobias[[i]][[j]][[k]]$s0.5, loc, h[k], s)
              results.fit.h0.5[j, i] <- getDifference.mal(fit) 
            }
          }
          if(k == 3){
            for(n in 1:4){
              fit <- getFitness(data = res.nobias[[i]][[j]][[k]]$s0.5, loc, h[k], s)
              results.fit.h1[j, i] <- getDifference.mal(fit) 
            }
          }
          if(k == 4){
            for(n in 1:4){
              fit <- getFitness(data = res.nobias[[i]][[j]][[k]]$s0.5, loc, h[k], s)
              results.fit.h99[j, i] <- getDifference.fem(fit) 
            }
          }
        }
      }
    }

rd <- c("0.02", "0.04", "0.06", "0.08", "0.1", "0.12", "0.14", "0.16", "0.18", "0.2", "0.22", 
        "0.24", "0.26", "0.28", "0.3", "0.32", "0.34", "0.36", "0.38", "0.4", "0.42", "0.44", 
        "0.46", "0.48")

results.fit.h0 <- cbind(results.fit.h0, rd)
results.fit.h0.5 <- cbind(results.fit.h0.5, rd)
results.fit.h1 <- cbind(results.fit.h1, rd)
# write .csv files to store your results permanently if desired. Remember to change the file name
# to avoid overwriting previously saved results.

write.csv(results.fit.h0, file ="h0.nobias.XYmodel.csv")
write.csv(results.fit.h0.5, file ="h0.5.nobias.XYmodel.csv")
write.csv(results.fit.h1, file ="h1.nobias.XYmodel.csv")

###