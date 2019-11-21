# Julio Rincones-Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com
# P.I. Dr. Heath Blackmon
# coleoguy@gmail.com

### this script is used for parsing the results of our simulations. 
### Specifically to extract the results
### corresponding to autosomes when females are the common sex.

#### I'm tempted to redo these to generate one tidyverse compliant data frame.
### The write.csv commands need changing to specify the directory where 
### I want these saved (preferably in the figures folder)

source("../functions/analysis.functions.R")
load("../results/base.comp.model.RData")
res.base <- results
load("../results/rare.male.model.RData")
res.rare.mal <- results
rm(results)

# create matrix for results

results.fit <- matrix(, 4, 4)

# insert column and row names

col.names <- c("50", "100", 
               "500", "1000")

colnames(results.fit) <- col.names

row.names <- c("OSR1.h0", "OSR1.h.5", "OSR1.h1", "OSR1.h99")

rownames(results.fit)<- row.names

# providing the formals for the getFitness function
loc <- "auto"
h <- c(0, .5, 1, .99)
s <- .5

i <- j <- k <- 1



# Autosomal loci evaluation
# cycle through population sizes (50, 100, 500, 1000)
for(i in 1:4){
  # pull results when rd = 0.5 only
  for(j in 1:4){
    if(names(res.base[[i]])[j] == "rd0.5"){
      # cycle through genetic architectures but only when s = 0.5
      for(k in 1:4){ 
        s.coeff <- which(names(res.base[[i]][[j]][[k]]) == "s0.5")
        fit <- getFitness(data = res.base[[i]][[j]][[k]][[s.coeff]], 
                          loc, h[k], s)
        results.fit[k, i] <- getDifference.fem(fit)
      }
    }
  }
}


# create matrix for results

results.fit.biash0 <- matrix(, 6, 4)

# insert column and row names
colnames(results.fit.biash0) <- c("50", "100", 
                                  "500", "1000")

rownames(results.fit.biash0) <- 
  c("OSR.8.h0", "OSR.6.h0", "OSR.4.h0", "OSR.2.h0", "OSR.1.h0", "OSR.05.h0")

# create matrix for results

results.fit.biash0.5 <- matrix(, 6, 4)

# insert column and row names
colnames(results.fit.biash0.5) <- c("50", "100", 
                                    "500", "1000")

rownames(results.fit.biash0.5) <- 
  c("OSR.8.h.5", "OSR.6.h.5", "OSR.4.h.5", 
    "OSR.2.h.5", "OSR.1.h.5", "OSR.05.h.5")

# create matrix for results

results.fit.biash1 <- matrix(, 6, 4)

# insert column and row names
colnames(results.fit.biash1) <- c("50", "100", 
                                  "500", "1000")

rownames(results.fit.biash1) <- 
  c("OSR.8.h1", "OSR.6.h1", "OSR.4.h1", "OSR.2.h1", "OSR.1.h1", "OSR.05.h1")

# create matrix for results

results.fit.biash99 <- matrix(, 6, 4)

# insert column and row names
colnames(results.fit.biash99) <- c("50", "100", 
                                   "500", "1000")

rownames(results.fit.biash99) <- 
  c("OSR.8.h99", "OSR.6.h99", "OSR.4.h99", 
    "OSR.2.h99", "OSR.1.h99", "OSR.05.h99")

# providing the formals for the getFitness function
loc <- "auto"
h <- c(0, .5, 1, .99)
s <- .5


# Autosomal loci evaluation
# structure of results is as follows: 
# res.rare.fem$males50$females0.8$rd0.5$h0$s0.5 
# (loop will start iterating in this order in the case of autosomes.)
# cycle through population sizes (50, 100, 500, 1000)
for(i in 1:4){
  # cycle through OSR bias
  for(j in 1:6){
    # pull results when rd = 0.5 only
    for(k in 1:4){
      if(names(res.rare.mal[[i]][[j]])[k] == "rd0.5"){
        # cycle through genetic architectures but only when s = 0.5
        for(m in 1:4){
          if(m == 1){
            for(n in 1:4){
              fit <- getFitness(data = res.rare.mal[[i]][[j]][[k]][[m]]$s0.5, 
                                loc, h[m], s)
              results.fit.biash0[j, i] <- getDifference.fem(fit) 
            }
          }
          if(m == 2){
            for(n in 1:4){
              fit <- getFitness(data = res.rare.mal[[i]][[j]][[k]][[m]]$s0.5, 
                                loc, h[m], s)
              results.fit.biash0.5[j, i] <- getDifference.fem(fit) 
            }
          }
          if(m == 3){
            for(n in 1:4){
              fit <- getFitness(data = res.rare.mal[[i]][[j]][[k]][[m]]$s0.5, 
                                loc, h[m], s)
              results.fit.biash1[j, i] <- getDifference.fem(fit) 
            }
          }
          if(m == 4){
            for(n in 1:4){
              fit <- getFitness(data = res.rare.mal[[i]][[j]][[k]][[m]]$s0.5, 
                                loc, h[m], s)
              results.fit.biash99[j, i] <- getDifference.fem(fit) 
            }
          }
        }
      }
    }
  }
}

row.names.osr <- 
  c("OSR.1", "OSR0.8", "OSR0.6", "OSR0.4", "OSR0.2", "OSR0.1", "OSR0.05")

results.fit <- t(results.fit)
results.nobias <- as.data.frame(results.fit)

results.nobias.h0 <- results.nobias$OSR1.h0
h0.autosomes <- rbind(results.nobias.h0, results.fit.biash0)
rownames(h0.autosomes) <- row.names.osr

results.nobias.h.5 <- results.nobias$OSR1.h.5
h.5.autosomes <- rbind(results.nobias.h.5, results.fit.biash0.5)
rownames(h.5.autosomes) <- row.names.osr

results.nobias.h1 <- results.nobias$OSR1.h1
h1.autosomes <- rbind(results.nobias.h1, results.fit.biash1)
rownames(h1.autosomes) <- row.names.osr

results.nobias.h99 <- results.nobias$OSR1.h99
h99.autosomes <- rbind(results.nobias.h99, results.fit.biash99)
rownames(h99.autosomes) <- row.names.osr

write.csv(h0.autosomes, file ="h0.autosomes.csv")
write.csv(h.5.autosomes, file ="h0.5.autosomes.csv")
write.csv(h1.autosomes, file ="h1.autosomes.csv")
write.csv(h99.autosomes, file ="h.99.autosomes.csv")


