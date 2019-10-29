# Julio Rincones Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com
# P.I. Dr. Heath Blackmon
# coleoguy@gmail.com

# This script simulates a subset of conditions in our XY model where bias ranges from 0.2 to 0.05
# and we test different recombination distances to determine more accurately the inflection point at
# which significant changes in the fitness of one sex takes place.

# More precisely, in these simulations females are the common sex and males are the rare sex.

# first we load our functions
source("../functions/functions.R")

# here we set up all the variables that we will explore in our 
# simulations below
# numbers of females
females <- c(500, 100)
# levels of bias to explore
bias <- c(.2, .1, 0.05)
# recombination distance
rd <- c(0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26,
        0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, .5)
# dominance factor of allele 1
h <- c(0, .5, 1, 99) # 99 indicates to run SSD model
# selection strengths
s <- c(.1,.2,.5,.9)
# number of trials to run
iter <- 1000

# allele 1 is fit for males
# allele 2 is fit for females

# this allows the program to be run on multiple CPUs
library(doMC)
registerDoMC(2)

# these are used mainly for trouble shooting
# normally these are iterated in the loops below
ii<- i <- j <- k <- m <- 1

# this loop sets up the structure of the results
# it creates a tree of nested lists with the terminal
# lists being matrices 3 by iter the columns of the matrices
# are X, Y, and A giving the frequency of alleles that 
# are present on those chromosomes
results <- list()
for(i in 1:length(females)){
  results[[i]] <- list()
  names(results)[i] <- paste("females", females[i], sep = "")[1]
  for(ii in 1:length(bias)){
    results[[i]][[ii]] <- list()
    names(results[[i]])[ii] <- paste("males", bias[ii], sep = "")[1]
    for(j in 1:length(rd)){
      results[[i]][[ii]][[j]] <- list()
      names(results[[i]][[ii]])[j] <- paste("rd", rd[j], sep = "")
      for(k in 1:length(h)){
        results[[i]][[ii]][[j]][[k]] <- list()
        names(results[[i]][[ii]][[j]])[k] <- paste("h", h[k], sep = "")
        for(m in 1:length(s)){
          results[[i]][[ii]][[j]][[k]][[m]] <- list()
          names(results[[i]][[ii]][[j]][[k]])[m] <- paste("s", s[m], sep = "")
        }
      }
    }
  }
}

# this runs the rare male conditions for our sim project
# i cycles through pop sizes
# ii cycles through levels of bias
# j cycles through recombination distances
# k cycles through dominance factors
# m cycles through strengths of selection
# n cycles through iterations/trials
for(i in 1:length(females)){
  cat("running pop size:", females[i], "\n")
  for(ii in 1:length(bias)){
    cat("running bias:", bias[ii], "\n")
    for(j in 1:length(rd)){
      cat("running rd:", rd[j], "\n")
      for(k in 1:length(h)){
        cat("running h:", h[k], "\n")
        for(m in 1:length(s)){
          cat("running s:", s[m], "\n")
          run.result <- matrix(,iter,3)
          colnames(run.result) <- c("X", "Y", "A")
          x <- foreach(n=1:iter, .combine = "rbind") %dopar%{
            # these next two lines handle figuring out how many males
            # we have based on female count and current bias level
            males.each <- round(females[i]*bias[ii]/4)
            males <- males.each*4
            # here we make the initial population with all equal
            # frequencies of each genotype in males and females
            pop <- makeGenomes(females, males, 
                               freqs = c(rep(females[i]/4, 4),
                                         rep(males.each, 4)))
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
              resultA[p] <- GetFreq(pop, chrom="A", allele = 1, females=females[i], males=males)
              resultY[p] <- GetFreq(pop, chrom="Y", allele = 1, females=females[i], males=males)
              resultX[p] <- GetFreq(pop, chrom="X", allele = 2, females=females[i], males=males)
              # this runs a generation of the simulation
              pop <- Generation(pop, females=females[i], males=males, rd=rd[j], h=h[k], s=s[m])
              # this checks to see if something is fixed in the pop
              if(resultY[p] %in% c(1,0) & resultX[p] %in% c(1,0)){
                segregating <- F
                A <- resultA[p]
                X <- resultX[p]
                Y <- resultY[p]
              }
              # this checks to see if we have run 2000 generations
              if(p > 500){
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
          results[[i]][[ii]][[j]][[k]][[m]] <- run.result
        }
      }
    }
  }
}


##
save(results, file = "../results/rd.rare.male.RData")


