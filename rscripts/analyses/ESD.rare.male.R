# Julio Rincones Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com
# P.I. Dr. Heath Blackmon
# coleoguy@gmail.com

library(doMC)
source("../functions/functions.esd.R")
# simulations for rare males

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
max.gens <- 100


results <- as.data.frame(matrix(NA,0,6))
colnames(results) <- c("freq0", "OSR", "sex.com", "num.com", "h", "s")

# these nested loops will test each scenario pairing different 
# parameters as appropriate
for(i in 1:length(comm.sex)){
  # just prints to let us know progress
  print(paste("working on common sex =", comm.sex[i]))
  for(j in 1:length(osr)){
    
    # here we calculate number of males and females in the pop
    females <- comm.sex[i]
    males <- round(comm.sex[i]*osr[j])
    
    # here we calculate the total number of chromosomes
    chroms <- 2 * males + 2 * females
    
    for(k in 1:length(s)){
      
      for(m in 1:length(h)){
        # how many cores to run on
        registerDoMC(5)
        x <- foreach (iter = 1:replicates, .combine = "c") %dopar% {
          
          # sets up the initial population
          pop <- GetInitialPop(females, males)
          fre <- GetFreq(pop,allele=0,males,females)
          # this flag and counter will allow us to break out of the
          # while loop if an allele fixes or we have run as long
          # as we want to allow it to run
          segregating <- T
          counter <- 1
          
          while(segregating){
            # one generation of selection
            pop <- Generation(pop, s=s[k], h=h[m], females, males)
            
            # calculate the frequency of 0 allele
            old.fre <- fre
            fre <- GetFreq(pop,allele=0,males,females)
            
            # test whether we have met stopping conditions
            if(fre == 0 | fre == 1 | counter == max.gens | round(old.fre, digits=5) == round(fre, digits=5)){
              segregating <- F
            }
            counter <- counter + 1
          }
          fre
        }    
        # x = results from one setting
        if(males > females) sex.com <- "male"
        if(males < females) sex.com <- "female"
        if(males == females) sex.com <- "equal"
        temp.result <- data.frame(x,
                                  rep(osr[j], length(x)),
                                  rep(sex.com, length(x)),
                                  rep(comm.sex[i], length(x)),
                                  rep(h[m], length(x)),
                                  rep(s[k], length(x)))
        colnames(temp.result) <- colnames(results)
        results <- rbind(results, temp.result)
      }
    }
  }
}





