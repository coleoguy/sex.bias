# Heath Blackmon
# coleoguy@gmail.com
# simulate baseline populations

# first we load our functions
source("../sex.bias/rscripts/functions/functions.R")

# here we set up all the variables that we will explore in our 
# simulations below
# numbers of females
females <- c(50, 100, 500, 1000)
males <- c(50, 100, 500, 1000)
# recombination distance
rd <- c(0, .1, .2, .5)
# dominance factor of allele 1
h <- c(0, .5, 1, 99) # 99 indicates to run SSD model
# strength of sexual selection
s <- c(.1,.2,.5,.9)
# number of trials to run
iter <- 500
# allele 1 is fit for males
# allele 2 is fit for females

# this allows the program below to use multiple cores
library(doMC)
registerDoMC(7)

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
rm(list=ls()[-18])
# save as base.comp.model.RData
