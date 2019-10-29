# Julio Rincones-Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com
# P.I. Dr. Heath Blackmon
# coleoguy@gmail.com

### this script is used for parsing the results of our simulations. Specifically to extract the results
### corresponding to both sex chromosomes.

#### I'm tempted to redo these to generate one tidyverse compliant data frame.
### The write.csv commands need changing to specify the directory where I want these saved (preferably 
### in the figures folder)


# load functions
source("../functions/analysis.functions.R")
load("../results/base.comp.model.RData")
res.base <- results
load("../results/rare.female.model.RData")
res.rare.fem <- results
load("../results/rare.male.model.RData")
res.rare.mal <- results
rm(results)

# providing the formals for the getFitness function
loc <- "sex"
h <- c(0, .5, 1, .99)
s <- .5

# create matrix for results when there is no bias

results.fit <- matrix(, 4, 4)

# insert column and row names
colnames(results.fit) <- c("50", "100", 
                           "500", "1000")

rownames(results.fit)<- c("OSR1.h0", "OSR1.h.5", "OSR1.h1", "OSR1.h99")



# Sexually antagonistic loci evaluation on sex chromosomes
# cycle through population sizes (50, 100, 500, 1000)
for(i in 1:4){
  # pull results when rd = 0.1 only
  for(j in 1:4){
    if(names(res.base[[i]])[j] == "rd0.1"){
      # cycle through genetic architectures but only when s = 0.5
      for(k in 1:4){ 
        s.coeff <- which(names(res.base[[i]][[j]][[k]]) == "s0.5")
        fit <- getFitness(data = res.base[[i]][[j]][[k]][[s.coeff]], loc, h[k], s)
        results.fit[k, i] <- getDifference.mal(fit)
      }
    }
  }
}


# create matrix for results when there is bias with genetic architectures h = 0, 0.5, 1, and  sex-
# specific dominance. Looking at the case when males are the common sex.

# h = 0
results.fit.biash0 <- matrix(, 6, 4)

colnames(results.fit.biash0) <- c("50", "100", 
                                  "500", "1000")

rownames(results.fit.biash0) <- c("OSR.8.h0", "OSR.6.h0", "OSR.4.h0", "OSR.2.h0", "OSR.1.h0", "OSR.05.h0")

# h = 0.5

results.fit.biash0.5 <- matrix(, 6, 4)

colnames(results.fit.biash0.5) <- c("50", "100", 
                                    "500", "1000")

rownames(results.fit.biash0.5) <- c("OSR.8.h.5", "OSR.6.h.5", "OSR.4.h.5", "OSR.2.h.5", "OSR.1.h.5", "OSR.05.h.5")

# h = 1

results.fit.biash1 <- matrix(, 6, 4)

colnames(results.fit.biash1) <- c("50", "100", 
                                  "500", "1000")

rownames(results.fit.biash1) <- c("OSR.8.h1", "OSR.6.h1", "OSR.4.h1", "OSR.2.h1", "OSR.1.h1", "OSR.05.h1")

# h =  .99

results.fit.biash99 <- matrix(, 6, 4)

colnames(results.fit.biash99) <- c("50", "100", 
                                   "500", "1000")

rownames(results.fit.biash99) <- c("OSR.8.h99", "OSR.6.h99", "OSR.4.h99", "OSR.2.h99", "OSR.1.h99", "OSR.05.h99")


# Evaluation of exually antagonistic loci on sex  chromosomes when there is bias
# structure of results is as follows: res.rare.fem$males50$females0.8$rd0.1$h0$s0.5 (loop will start iterating in this
# order in the case of sex chromosomes.)
# cycle through population sizes (50, 100, 500, 1000)  
for(i in 1:4){
  # cycle through OSR bias in this order: 0.8, 0.6, 0.4, 0.2, 0.1, 0.05
  for(j in 1:6){
     # pull results when rd = 0.1 only
     for(k in 1:4){
      if(names(res.rare.fem[[i]][[j]])[k] == "rd0.1"){
        # cycle through genetic architectures. The selection coefficient is hard coded. If you want to choose a different
        # coefficient you need to change '$s0.5' to any of the values we have used for the simulations(0.1, 0.2, 0.5, 0.9).
        for(m in 1:4){
          if(m == 1){
            for(n in 1:4){
              fit <- getFitness(data = res.rare.fem[[i]][[j]][[k]][[m]]$s0.5, loc, h[m], s)
              results.fit.biash0[j, i] <- getDifference.mal(fit) 
             }
          }
          if(m == 2){
            for(n in 1:4){
              fit <- getFitness(data = res.rare.fem[[i]][[j]][[k]][[m]]$s0.5, loc, h[m], s)
              results.fit.biash0.5[j, i] <- getDifference.mal(fit) 
             }
          }
          if(m == 3){
            for(n in 1:4){
               fit <- getFitness(data = res.rare.fem[[i]][[j]][[k]][[m]]$s0.5, loc, h[m], s)
               results.fit.biash1[j, i] <- getDifference.mal(fit) 
             }
           }
           if(m == 4){
             for(n in 1:4){
               fit <- getFitness(data = res.rare.fem[[i]][[j]][[k]][[m]]$s0.5, loc, h[m], s)
               results.fit.biash99[j, i] <- getDifference.mal(fit) 
             }
           }
         }
       }
     }
   }
 }
# Since the results of these  loops are one table for difference in relative fitness when there is no bias, and
# four tables for when there is bias (one for each genetic architecture), need to extract the relevant row from the
# "no bias" table for each genetic architecture, and then row bind the matching "bias" table.

results.fit <- t(results.fit)
results.nobias <- as.data.frame(results.fit)

results.nobias.h0 <- results.nobias$OSR1.h0
h0.rfem.sexchrom <- rbind(results.nobias.h0, results.fit.biash0)
rownames(h0.rfem.sexchrom) <- c(1, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05)

results.nobias.h0.5 <- results.nobias$OSR1.h.5
h0.5.rfem.sexchrom <- rbind(results.nobias.h0.5, results.fit.biash0.5)
rownames(h0.5.rfem.sexchrom) <- c(1, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05)

results.nobias.h1 <- results.nobias$OSR1.h1
h1.rfem.sexchrom <- rbind(results.nobias.h1, results.fit.biash1)
rownames(h1.rfem.sexchrom) <- c(1, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05)

results.nobias.h99 <- results.nobias$OSR1.h99
h99.rfem.sexchrom <- rbind(results.nobias.h99, results.fit.biash99)
rownames(h99.rfem.sexchrom) <- c(1, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05)

write.csv(h0.rfem.sexchrom, file ="h0.rfem.sexchrom.csv")
write.csv(h0.5.rfem.sexchrom, file ="h0.5.rfem.sexchrom.csv")
write.csv(h1.rfem.sexchrom, file ="h1.rfem.sexchrom.csv")
write.csv(h99.rfem.sexchrom, file ="h.99.rfem.sexchrom.csv")


##
# Now let's look at the case when females are the common sex.
# create matrix for results when there is bias with genetic architectures h = 0, 0.5, 1, and  sex-
# specific dominance. 

# create matrix for results when there is no bias

results.fit <- matrix(, 4, 4)

# insert column and row names
colnames(results.fit) <- c("50", "100", 
                           "500", "1000")

rownames(results.fit)<- c("OSR1.h0", "OSR1.h.5", "OSR1.h1", "OSR1.h99")


# Need to run through this loop once more since the getDifference function for the case where
# females are the common sex, is different.
# Sexually antagonistic loci evaluation on sex chromosomes
# cycle through population sizes (50, 100, 500, 1000)
for(i in 1:4){
  # pull results when rd = 0.1 only
  for(j in 1:4){
    if(names(res.base[[i]])[j] == "rd0.1"){
      # cycle through genetic architectures but only when s = 0.5
      for(k in 1:4){ 
        s.coeff <- which(names(res.base[[i]][[j]][[k]]) == "s0.5")
        fit <- getFitness(data = res.base[[i]][[j]][[k]][[s.coeff]], loc, h[k], s)
        results.fit[k, i] <- getDifference.fem(fit)
      }
    }
  }
}


# h = 0
results.fit.biash0 <- matrix(, 6, 4)

colnames(results.fit.biash0) <- c("50", "100", 
                                  "500", "1000")

rownames(results.fit.biash0) <- c("OSR.8.h0", "OSR.6.h0", "OSR.4.h0", "OSR.2.h0", "OSR.1.h0", "OSR.05.h0")

# h = 0.5

results.fit.biash0.5 <- matrix(, 6, 4)

colnames(results.fit.biash0.5) <- c("50", "100", 
                                    "500", "1000")

rownames(results.fit.biash0.5) <- c("OSR.8.h.5", "OSR.6.h.5", "OSR.4.h.5", "OSR.2.h.5", "OSR.1.h.5", "OSR.05.h.5")

# h = 1

results.fit.biash1 <- matrix(, 6, 4)

colnames(results.fit.biash1) <- c("50", "100", 
                                  "500", "1000")

rownames(results.fit.biash1) <- c("OSR.8.h1", "OSR.6.h1", "OSR.4.h1", "OSR.2.h1", "OSR.1.h1", "OSR.05.h1")

# h =  .99

results.fit.biash99 <- matrix(, 6, 4)

colnames(results.fit.biash99) <- c("50", "100", 
                                   "500", "1000")

rownames(results.fit.biash99) <- c("OSR.8.h99", "OSR.6.h99", "OSR.4.h99", "OSR.2.h99", "OSR.1.h99", "OSR.05.h99")


# Evaluation of exually antagonistic loci on sex  chromosomes when there is bias
# structure of results is as follows: res.rare.fem$males50$females0.8$rd0.1$h0$s0.5 (loop will start iterating in this
# order in the case of sex chromosomes.)
# cycle through population sizes (50, 100, 500, 1000)
for(i in 1:4){
  # cycle through OSR bias in this order: 0.8, 0.6, 0.4, 0.2, 0.1, 0.05
  for(j in 1:6){
    # pull results when rd = 0.5 only
    for(k in 1:4){
      if(names(res.rare.mal[[i]][[j]])[k] == "rd0.1"){
        # cycle through genetic architectures but only when s = 0.5
        for(m in 1:4){
          if(m == 1){
            for(n in 1:4){
              fit <- getFitness(data = res.rare.mal[[i]][[j]][[k]][[m]]$s0.5, loc, h[m], s)
              results.fit.biash0[j, i] <- getDifference.fem(fit) 
            }
          }
          if(m == 2){
            for(n in 1:4){
              fit <- getFitness(data = res.rare.mal[[i]][[j]][[k]][[m]]$s0.5, loc, h[m], s)
              results.fit.biash0.5[j, i] <- getDifference.fem(fit) 
            }
          }
          if(m == 3){
            for(n in 1:4){
              fit <- getFitness(data = res.rare.mal[[i]][[j]][[k]][[m]]$s0.5, loc, h[m], s)
              results.fit.biash1[j, i] <- getDifference.fem(fit) 
            }
          }
          if(m == 4){
            for(n in 1:4){
              fit <- getFitness(data = res.rare.fem[[i]][[j]][[k]][[m]]$s0.5, loc, h[m], s)
              results.fit.biash99[j, i] <- getDifference.fem(fit) 
            }
          }
        }
      }
    }
  }
}

# Since the results of these  loops are one table for difference in relative fitness when there is no bias, and
# four tables for when there is bias (one for each genetic architecture), need to extract the relevant row from the
# "no bias" table for each genetic architecture, and then row bind the matching "bias" table.

results.fit <- t(results.fit)
results.nobias <- as.data.frame(results.fit)

results.nobias.h0 <- results.nobias$OSR1.h0
h0.rmal.sexchrom <- rbind(results.nobias.h0, results.fit.biash0)
rownames(h0.rmal.sexchrom) <- c(1, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05)

results.nobias.h0.5 <- results.nobias$OSR1.h.5
h0.5.rmal.sexchrom <- rbind(results.nobias.h0.5, results.fit.biash0.5)
rownames(h0.5.rmal.sexchrom) <- c(1, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05)

results.nobias.h1 <- results.nobias$OSR1.h1
h1.rmal.sexchrom <- rbind(results.nobias.h1, results.fit.biash1)
rownames(h1.rmal.sexchrom) <- c(1, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05)

results.nobias.h99 <- results.nobias$OSR1.h99
h99.rmal.sexchrom <- rbind(results.nobias.h99, results.fit.biash99)
rownames(h99.rmal.sexchrom) <- c(1, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05)

write.csv(h0.rmal.sexchrom, file ="h0.rmal.sexchrom.csv")
write.csv(h0.5.rmal.sexchrom, file ="h0.5.rmal.sexchrom.csv")
write.csv(h1.rmal.sexchrom, file ="h1.rmal.sexchrom.csv")
write.csv(h99.rmal.sexchrom, file ="h.99.rmal.sexchrom.csv")
