## Julio Rincones Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com
# P.I. Dr. Heath Blackmon

# This script runs simulations for instances where there is no OSR bias but with recombination
# distance values ranging between 0.02 and 0.48. The results of which we parse to determine
# the impact of rd values and genetic architecture on the fate of sexually antagonistic variation

# first we load our functions
source("../functions/functions.R")

# here we set up all the variables that we will explore in our simulations below
# we have chosen to focus on intermediate values for numbers of females and number of males.

females <- c(100, 500)
males <- c(100, 500)

# recombination distance
rd <- c(0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26,
        0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48)

# dominance factor of allele 1
h <- c(0, .5, 1, 99) # 99 indicates to run SSD model

# strength of sexual selection (i.e. the selection coefficient)
s <- c(.1,.2,.5,.9)

# number of trials to run
iter <- 500

# allele 1 is fit for males
# allele 2 is fit for females

# this allows the program below to use multiple CPUs
library(doMC)
registerDoMC(4)

# this is really just for testing usually these are incremented
# in the loops described below
i <- j <- k <- m <- 1

# this loop sets up the structure of the results
# it creates a tree of nested lists with the terminal
# lists being matrices 3 by iter the columns of the matrices
# are X, Y, and A giving the frequency of alleles that 
# are present on those chromosomes
results <- list()
for(i in 1:length(females)){
  results[[i]] <- list()
  names(results)[i] <- paste("pop", max(females[i], males[i]), sep = "")[1]
  for(j in 1:length(rd)){
    results[[i]][[j]] <- list()
    names(results[[i]])[j] <- paste("rd", rd[j], sep = "")
    for(k in 1:length(h)){
      results[[i]][[j]][[k]] <- list()
      names(results[[i]][[j]])[k] <- paste("h", h[k], sep = "")
      for(m in 1:length(s)){
        results[[i]][[j]][[k]][[m]] <- list()
        names(results[[i]][[j]][[k]])[m] <- paste("s", s[m], sep = "")
      }
    }
  }
}

# this runs the base comparison conditions for our sim project
# i cycles through pop sizes
# j cycles through recombination distances
# k cycles through dominance factors
# m cycles through strengths of selection
# n cycles through iterations/trials
for(i in 1:length(females)){
  cat("running pop size:", females[i], "\n")
  for(j in 1:length(rd)){
    cat("running rd:", rd[j], "\n")
    for(k in 1:length(h)){
      cat("running h:", h[k], "\n")
      for(m in 1:length(s)){
        cat("running s:", s[m], "\n")
        run.result <- matrix(,iter,3)
        colnames(run.result) <- c("X", "Y", "A")
        x <- foreach(n=1:iter, .combine = "rbind") %dopar%{
          # here we make the initial population with all equal
          # frequencies of each genotype in males and females
          pop <- makeGenomes(females, males, 
                             freqs = c(rep(females[i]/4, 4),
                                       rep(males[i]/4, 4)))
          # this sets up the contain for a given sim
          resultY <- resultX <- resultA <- c()
          # this will be true unless an allele has fixed as the SAL
          segregating <- T
          # this is just a counter
          p <- 1
          # this while loop will run till something fixes or 
          # until we reach 2000 generations
          while(segregating){
            print(p)
            # this gets the allele frequencies we are interested in
            resultA[p] <- GetFreq(pop, chrom="A", allele = 1, females=females[i], males=males[i])
            resultY[p] <- GetFreq(pop, chrom="Y", allele = 1, females=females[i], males=males[i])
            resultX[p] <- GetFreq(pop, chrom="X", allele = 2, females=females[i], males=males[i])
            # this runs a generation of the simulation
            pop <- Generation(pop, females=females[i], males=males[i], rd=rd[j], h=h[k], s=s[m])
            # this checks to see if something is fixed in the pop
            if(resultY[p] %in% c(1,0) & resultX[p] %in% c(1,0)){
              segregating <- F
              A <- resultA[p]
              X <- resultX[p]
              Y <- resultY[p]
            }
            # this checks to see if we have run 1000 generations
            if(p > 1000){
              segregating <- F
              A <- resultA[p]
              X <- resultX[p]
              Y <- resultY[p]
            }
            # this increments our counter
            p <- p + 1
          }
          # this puts the final result of our sim together
          c(X, Y, A)
        }
        # we have now left the parallel loop and the variable
        # x contains all of our results
        run.result[1:iter, 1:3] <- x
        results[[i]][[j]][[k]][[m]] <- run.result
      }
    }
  }
}

# save as osr1.rd.model.RData
