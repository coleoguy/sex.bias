# Julio Rincones-Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com
# P.I. Dr. Heath Blackmon
# coleoguy@gmail.com


# This script has functions for the analysis of results
# The intention being to assess the impact of allele frequency 
# change in our simulations and how this
# leads to feminization and masculinization of the genome.


getFitness <- function(data, loc, h, s){
  results <- matrix(,nrow(data),2)
  colnames(results) <- c("wbar.m", "wbar.f")
  for(i in 1:nrow(data)){
    if(loc=="auto"){
      p <- data[i, 3]
      a11 <- p ^ 2
      a12 <- 2 * p * (1 - p)
      a22 <- (1 - p) ^ 2
      wm.a11 <- 1 + s
      wm.a12 <- 1 + h * s
      wm.a22 <- 1
      wf.a11 <- 1 / (1 + s)
      wf.a12 <- 1 / (1 + h * s)
      wf.a22 <- 1
      # standardise by dividing the mean fitness by the highest 
      # fitness each sex can have yielding relative fitnesses.
      wmean.m <- (a11 * wm.a11 + a12 * wm.a12 + a22 * wm.a22)/wm.a11
      wmean.f <- (a11 * wf.a11 + a12 * wf.a12 + a22 * wf.a22)/wf.a22
    }
    if(loc=="sex"){
      px <- data[i, 1]
      py <- data[i, 2]
      m.a11 <- px * py
      m.a12 <- px * (1 - py) + py * (1 - px)
      m.a22 <- (1 - px) * (1 - py)
      f.a11 <- px ^ 2
      f.a12 <- 2 * px * (1 - px)
      f.a22 <- (1 - px) ^ 2
      wm.a11 <- 1 + s
      wm.a12 <- 1 + h * s
      wm.a22 <- 1
      wf.a11 <- 1 / (1 + s)
      wf.a12 <- 1 / (1 + h * s)
      wf.a22 <- 1
      # standardise by dividing the mean fitness by the highest fitness 
      # each sex can have yielding relative fitnesses.
      wmean.m <- (m.a11 * wm.a11 + m.a12 * wm.a12 + m.a22 * wm.a22)/wm.a11
      wmean.f <- (f.a11 * wf.a11 + f.a12 * wf.a12 + f.a22 * wf.a22)/wf.a22
    }
    results[i, 1:2] <- c(wmean.m, wmean.f)
  }
  return(as.data.frame(results))
}
  
# Function to get the fitness difference when males are the common sex
getDifference.mal <- function(w.beta){
  mean(fit$wbar.m - fit$wbar.f)
}

# Function to get the fitness difference when females are the common sex

getDifference.fem <- function(w.beta){
  mean(fit$wbar.f - fit$wbar.m)
}

