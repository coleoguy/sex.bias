# Julio Rincones Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com
# P.I. Dr. Heath Blackmon
# coleoguy@gmail.com

# This is code to look at the impact of operational
# sex ratio bias and its impact on sexual antagonism
# in haplodiploid species

# single locus biallelic model with locus having
# an allele that benefits females and an allele
# that benefits males


# allele 0 good for females
# allele 1 good for males

# sets up initial genomes
GetInitialPop <- function(females, males){
  # females = number of females in the population
  # males = number of males in the population
  
  # create a vector for female genotypes
  pop.fem <- vector(length = 3)
  names(pop.fem) <- c("00", "01", "11")
  
  # create a vector for male genotypes
  pop.mal <- vector(length = 2)
  names(pop.mal) <- c("0","1")

  # determine the largest value that can go into
  # each of the possible genotypes
  pop.fem[1:3] <- floor(females/3)
  pop.mal[1:2] <- floor(males/2)
  
  # determine any "leftover" individuals when
  # we try to split the popsize into genotypes
  extra.fem <- females%%3
  extra.mal <- males%%2
  
  # add in those "leftover" individuals to 
  # genotypes randomly
  if(extra.fem==1){
    pick <- sample(1:3, 1)
    pop.fem[pick] <- pop.fem[pick] + 1
  }
  if(extra.fem==2){
    pick <- sample(1:3, 2)
    pop.fem[pick] <- pop.fem[pick] + 1
  }
  if(extra.mal==1){
    pick <- sample(1:2, 1)
    pop.mal[pick] <- pop.mal[pick] + 1
  }
  
  result <- list(pop.fem, pop.mal)
  names(result) <- c("pop.fem","pop.mal")
  return(result)
}

# run a generation
Generation <- function(pop, s=.5,h=.5, females, males){
  # assess fitness
  GetFitness <- function(pop, s, h){
    # pop = list of two vectors containing genotype counts
    # s = selection coefficient
    # h = dominance coefficient for 0 allele
    
    # absolute fitnesses
    fits.fem <- c(1+s, 1+h*s, 1)
    fits.mal <- c(1/(1+s), 1)
    
    # get mean fitness
    mean.fit.fem <- sum(pop$pop.fem * fits.fem) / sum(pop$pop.fem)
    mean.fit.mal <- sum(pop$pop.mal * fits.mal) / sum(pop$pop.mal)
    
    # get relative fitness
    rel.fit.fem <- fits.fem/mean.fit.fem
    rel.fit.mal <- fits.mal/mean.fit.mal
    
    results <- list(rel.fit.fem, rel.fit.mal)
    names(results) <- c("rel.fit.fem", "rel.fit.mal")
    return(results)
  }
  
  # makes gametes
  GetGametes <- function(fits, females, males, pop){
    # sample female genotypes based on fitness
    fem.gens <- c(rep("00", pop$pop.fem[1]),
                  rep("01", pop$pop.fem[2]),
                  rep("11", pop$pop.fem[3]))
    fem.fits <- c(rep(fits$rel.fit.fem[1], pop$pop.fem[1]),
                  rep(fits$rel.fit.fem[2], pop$pop.fem[2]),
                  rep(fits$rel.fit.fem[3], pop$pop.fem[3]))
        fems <- sample(fem.gens, prob = fem.fits, 
                   replace = T, 
                   size = females)
        
    # sample male genotypes based on fitness
        mal.gens <- c(rep("0", pop$pop.mal[1]),
                      rep("1", pop$pop.mal[2]))
        mal.fits <- c(rep(fits$rel.fit.mal[1], pop$pop.mal[1]),
                      rep(fits$rel.fit.mal[2], pop$pop.mal[2]))
        mals <- sample(mal.gens, prob = mal.fits, 
                   replace = T, size = males)
    
    # record number of each haplotype that will be produced
    # by the selected genotypes
    eggs <- round(c(sum(fems == "00") + sum(fems == "01") * 0.5,
                    sum(fems == "11") + sum(fems == "01") * 0.5))
    sperm <- c(sum(mals == "0"),
               sum(mals == "1"))
    
    # format the output
    gametes <- rbind(eggs,sperm)
    colnames(gametes) <- c("0","1")
    gametes <- t(gametes)
    gametes <- as.data.frame(gametes)
    return(gametes)
  }
  
  # fertilization of next generation
  GetNextGen <- function(pop, gametes, females, males){
    # sample sperm haplotypes
    new.mal.hap <- sample(0:1, size = males, 
                          prob = gametes$eggs, 
                          replace = T)
    # sample egg haplotypes
    new.fem.gen <- paste(sample(0:1, size = females,
                                prob = gametes$eggs,
                                replace = T),
                         sample(0:1, size = females, 
                                prob = gametes$sperm, 
                                replace = T),
                         sep = "")
    # record the fertilized genotypes and haplotypes
    pop$pop.fem[1] <- sum(new.fem.gen == "00")
    pop$pop.fem[2] <- sum(new.fem.gen == "01") + 
      sum(new.fem.gen == "10") 
    pop$pop.fem[3] <- sum(new.fem.gen == "11")
    pop$pop.mal[1] <- sum(new.mal.hap == "0")
    pop$pop.mal[2] <- sum(new.mal.hap == "1")
    return(pop)
  }
  
  
  fits <- GetFitness(pop, s = s, h = h)
  gametes <- GetGametes(fits, females, males, pop)
  pop <- GetNextGen(pop, gametes, females, males)
  return(pop)
}
