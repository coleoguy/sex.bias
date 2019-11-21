# Julio Rincones Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com
# P.I. Dr. Heath Blackmon
# coleoguy@gmail.com

# This script simulates instances in our haplodiploidy model where males 
# are the common sex 

library(doMC)
source("functions.HD.R")
# simulations for rare females

# allele 0 good for females
# allele 1 good for males
# pop = list of two vectors containing genotype counts
# s = selection coefficient
# h = dominance coefficient for 0 allele

# repeat
comm.sex <- c(1000, 500, 100, 50)
osr <- c(1, .8, .6, .4, .2, .1,.05)
s <- c(0.1, 0.2, 0.5, 0.9)
h <- c(0.0, 0.5, 1.0)
replicates <- 1000
max.gens <- 500

# build the structure of the results as a tree of lists with names
# that describe the parameter values the terminal lists in this 
# will be tables with rows equal to number iterations and columns
# equal to different dominance (h) values.
results <- list()
for(i in 1:length(comm.sex)){
  results[[i]] <- list()
  names(results)[i] <- paste("comm.sex", comm.sex[i], sep="")
  for(j in 1:length(osr)){
    results[[i]][[j]] <- list()
    names(results[[i]])[j] <- paste("osr", osr[j], sep="")
    for(k in 1:length(s)){
      results[[i]][[j]][[k]] <- list()
      names(results[[i]][[j]])[k] <- paste("s", s[k], sep="")
    }
  }
}


# these nested loops will test each scenario pairing different 
# parameters as appropriate
for(i in 1:length(comm.sex)){
  # just prints to let us know progress
  print(paste("working on common sex =", comm.sex[i]))
  for(j in 1:length(osr)){
    
    # here we calculate number of males and females in the pop
    females <- round(comm.sex[i]*osr[j])
    males <- comm.sex[i]
    
    # here we calculate the total number of chromosomes
    chroms <- males + 2 * females
    
    for(k in 1:length(s)){
      
      # this sets up the tables that will hold a set of results
      cur.result <- as.data.frame(matrix(NA,replicates, length(h)))
      colnames(cur.result) <- c("h=0", "h=0.5", "h=1")
      
      for(m in 1:length(h)){
        # how many cores to run on
        registerDoMC(2)
        x <- foreach (iter = 1:replicates, .combine = "c") %dopar% {

          # sets up the initial population
          pop <- GetInitialPop(females, males)
          
          # this flag and counter will allow us to break out of the
          # while loop if an allele fixes or we have run as long
          # as we want to allow it to run
          segregating <- T
          counter <- 1
          
          while(segregating){
            # one generation of selection
            pop <- Generation(pop, s=s[k], h=h[m], females, males)
            
            # calculate the frequency of 0 allele. i.e. female beneficial
            fre <- (pop$pop.mal[1] + pop$pop.fem[1]*2 + pop$pop.fem[2]) / chroms
            
            # test whether we have met stopping conditions
            if(fre == 0 | fre == 1 | counter == max.gens){
              segregating <- F
            }
            counter <- counter + 1
          }
          # store the result of a single simulation
          fre
        }
        cur.result[, m] <- x
        }
      # store a table of results for 112,000 simulations
      results[[i]][[j]][[k]] <- cur.result
    }
    }
  }





