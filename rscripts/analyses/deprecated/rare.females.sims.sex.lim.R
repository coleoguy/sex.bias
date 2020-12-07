# Heath Blackmon
# coleoguy@gmail.com
# TODO add mutation function such that 1/1000 alleles suffers mutation


# first we load our functions
source("../functions/functions.sex.lim.R")

# here we set up all the variables that we will explore in our 
# simulations below
# numbers of males
males <- c(50, 100, 500, 1000)
# levels of bias to explore
bias <- c(1, .8,.6,.4,.2,.1,.05)
# dominance factor of allele 1
h <- c(.5) #h <- c(0, .5, 1)
# selection strengths
s <- -c(0.025, 0.05, .1) 
# number of trials to run
iter <- 1000
gen.counter <- c()
# allele 1 is fit for males
# allele 2 is fit for females

# this allows the program to be run on multiple CPUs
library(doMC)
registerDoMC(6)

# these are used mainly for trouble shooting
# normally these are iterated in the loops below
# ii<- i <- j <- k <- m <- 1

# this sets up the structure of the results
results <- as.data.frame(matrix(NA,0,9))
colnames(results) <- c("X","Y","A","males","OSR", "rd","h","s", "gens")

# this runs the rare male conditions for our sim project
# i cycles through pop sizes
# ii cycles through levels of bias
# j cycles through recombination distances
# k cycles through dominance factors
# m cycles through strengths of selection
# n cycles through iterations/trials
result.counter <- 1
for(i in 1:length(males)){
  cat("running pop size:", males[i], "\n")
  for(ii in 1:length(bias)){
    cat("running bias:", bias[ii], "\n")
    for(k in 1:length(h)){
      cat("running h:", h[k], "\n")
      for(m in 1:length(s)){
        cat("running s:", s[m], "\n")
        x <- foreach(n=1:iter, .combine = "rbind") %dopar%{
          # these next two lines handle figuring out how many males
          # we have based on female count and current bias level
          females <- round(males[i]*bias[ii])
          # here we make the initial population with all equal
          # frequencies of each genotype in males and females
          pop <- makeGenomes(females, males, 
                             freqs = c(c(females,0,0,0),
                                       c(males[i],0,0,0)))
          # this sets up the contain for a given sim
          resultA <- c()
          # this will be true unless an allele has fixed as the SAL
          segregating <- T
          # this is just a counter
          c.ite <- 1
          # this while loop will run till something fixes or 
          # until we reach 2000 generations
          while(segregating){
            #print(p)
            # this gets the allele frequencies we are interested in
            resultA[c.ite] <- GetFreq(pop, chrom="A", allele = 1, 
                                      females=females, males=males[i])
            # this runs a generation of the simulation
            pop <- Generation(pop, females=females, males=males[i], 
                              rd=.5, h=h[k], s=s[m],mut=mut.rate)
            if(c.ite > 1000){
              segregating <- F
              A <- resultA[c.ite]
            }
            # this increments our counter
            c.ite <- c.ite + 1
          }
          # this puts the final result of our sim together
          c(A, c.ite)
        }
        # we have now left the parallel loop and the variable
        # x contains all of our results
        results[result.counter:(result.counter+iter-1), 1] <- NA
        results[result.counter:(result.counter+iter-1), 2] <- NA
        results[result.counter:(result.counter+iter-1), 3] <- 1 - as.vector(x[,3])
        results[result.counter:(result.counter+iter-1), 4] <- males[i]
        results[result.counter:(result.counter+iter-1), 5] <- bias[ii]
        results[result.counter:(result.counter+iter-1), 6] <- .5
        results[result.counter:(result.counter+iter-1), 7] <- h[k]
        results[result.counter:(result.counter+iter-1), 8] <- s[m]
        results[result.counter:(result.counter+iter-1), 9] <- as.vector(x[,4])
        result.counter <- result.counter+iter
      }
    }
  }
}


