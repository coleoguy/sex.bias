# Julio Rincones Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com
# P.I. Dr. Heath Blackmon
# coleoguy@gmail.com
# simulating populations with sex bias
# and sexual antagonism and ESD


# sets up initial genomes
GetInitialPop <- function(females, males){
  # females = number of females in the population
  # males = number of males in the population
  
  # create a vector for female genotypes
  pop.fem <- vector(length = 3)
  names(pop.fem) <- c("00", "01", "11")
  
  # create a vector for male genotypes
  pop.mal <- vector(length = 3)
  names(pop.mal) <- c("00", "01", "11")
  
  # determine the largest value that can go into
  # each of the possible genotypes
  pop.fem[1:3] <- floor(females/3)
  pop.mal[1:3] <- floor(males/3)
  
  # determine any "leftover" individuals when
  # we try to split the popsize into genotypes
  extra.fem <- females%%3
  extra.mal <- males%%3
  
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
    pick <- sample(1:3, 1)
    pop.mal[pick] <- pop.mal[pick] + 1
  }
  if(extra.mal==2){
    pick <- sample(1:3, 2)
    pop.mal[pick] <- pop.mal[pick] + 1
  }
  result <- c(pop.fem[1:2], 0, pop.fem[3], 
              pop.mal[1:2], 0, pop.mal[3])
  return(result)
}



makeGenomes <- function(females, males, freqs=NULL){
  population <- rep(0, 8)
  if(!is.null(freqs)){
    population <- freqs
  }else{
    print("supply frequencies")
  }
  names(population) <- c("fem.11", "fem.12", "fem.21", "fem.22",
                         "mal.11", "mal.12", "mal.21", "mal.22")
  return(population)
}

measureFit <- function(pop, h, s){
  if(h==0){
    fit <- c(1,   1+s,   1+s,   1+s,
             1+s, 1, 1, 1)
    return(fit)
  }
  if(h==.5){
    fit <- c(1, 1+h*s, 1+h*s, 1+s,
             1+s,   1+h*s, 1+h*s, 1)
    return(fit)
  }
  if(h==1){
    fit <- c(1, 1,   1,   1+s,
             1+s,   1+s, 1+s, 1)
    return(fit)
  }
  if(h==99){
    fit <- c(1, 1+s, 1+s, 1+s,
             1+s, 1+s, 1+s, 1)
    return(fit)
  }
}

GetParentsGeno <- function(pop, fit, females, males){
  x.mom <- c(rep(1, pop[1]), rep(2, pop[2]), rep(3, pop[3]), rep(4, pop[4]))
  x.dad <- c(rep(5, pop[5]), rep(6, pop[6]), rep(7, pop[7]), rep(8, pop[8]))
  fit.mom <- c(rep(fit[1], pop[1]), rep(fit[2], pop[2]), 
               rep(fit[3], pop[3]), rep(fit[4], pop[4]))
  fit.dad <- c(rep(fit[5], pop[5]), rep(fit[6], pop[6]), 
               rep(fit[7], pop[7]), rep(fit[8], pop[8]))
  # this effectively performs viability selection
  # so we sample females as moms based on fitness
  mom.genomes <- sample(x = x.mom, 
                        size = females,
                        replace = T,
                        prob=fit.mom)
  # THIS IS WHERE THE PROBLEM IS WE ARE SAMPLING BASED ON FITNESS BUT EQUAL NUMBERS OF EACH
  # then we sample males as dads based on fitness
  dad.genomes <- sample(x.dad, 
                        size = males,
                        replace = T,
                        prob=fit.dad)
  # we can simplify our coding a bit by taking into consideration that
  # there are 4 genotypes of moms and 4 genotypes of dads lets make a
  # table that allows us to draw sperm and eggs based on the genotype 
  # distribution of parents and we will account for recombination in 
  # this step as well.
  mom.geno <- rep(0, 4)
  names(mom.geno) <- c("X1X1", "X1X2", "X2X1", "X2X2")
  mom.geno[1] <- sum(mom.genomes==1)
  mom.geno[2] <- sum(mom.genomes==2)
  mom.geno[3] <- sum(mom.genomes==3)
  mom.geno[4] <- sum(mom.genomes==4)
  dad.geno <- rep(0, 4)
  names(dad.geno) <- c("X1Y1", "X1Y2", "X2Y1", "X2Y2")
  dad.geno[1] <- sum(dad.genomes==5)
  dad.geno[2] <- sum(dad.genomes==6)
  dad.geno[3] <- sum(dad.genomes==7)
  dad.geno[4] <- sum(dad.genomes==8)
  
  parents <- list(mom.geno, dad.geno)
  names(parents) <- c("moms", "dads")
  return(parents)
}

makeEggs <- function(mom.geno){
  eggs <- rep(0, 2)
  names(eggs) <- c("X1", "X2")
  eggs[1] <- mom.geno[1] +
    0.5 * mom.geno[2] +
    0.5 * mom.geno[3]
  eggs[2] <- mom.geno[4] + 
    0.5 * mom.geno[2] +
    0.5 * mom.geno[3]
  return(eggs)
}

makeSperm <- function(dad.geno, rd){
  sperm <- rep(0, 2)
  names(sperm) <- c("1", "2")
  sperm[1] <- dad.geno[1] +
    0.5 * dad.geno[2] +
    0.5 * dad.geno[3] 
  sperm[2] <- dad.geno[4] + 
    0.5 * dad.geno[2] +
    0.5 * dad.geno[3] 
  return(sperm)
}

makeNewPop <- function(pop, eggs, sperm, females, males){
  x  <-  paste(sample(1:2, size = females, replace=T, prob = eggs),
               sample(1:2, size = females, replace=T, prob = sperm))
  pop[1] <- sum(x == "1 1")
  pop[2] <- sum(x == "1 2")
  pop[3] <- sum(x == "2 1")
  pop[4] <- sum(x == "2 2")
  x <- paste(sample(1:2, size = males, replace=T, prob = eggs),
             sample(1:2, size = males, replace=T, prob = sperm))
  pop[5] <- sum(x == "1 1")
  pop[6] <- sum(x == "1 2")
  pop[7] <- sum(x == "2 1")
  pop[8] <- sum(x == "2 2")
  return(pop)
}

Generation <- function(pop, females, males, h, s){
  fit <- measureFit(pop, h, s)
  parents <- GetParentsGeno(pop, fit, females, males)
  eggs <- makeEggs(parents$moms)
  sperm <- makeSperm(parents$dads, rd)
  pop <- makeNewPop(pop, eggs, sperm, females, males)
  return(pop)
}

GetFreq <- function(pop, allele, males, females){
    zeros <- (pop[1]*2 + pop[2] + pop[3] +
             pop[5]*2 + pop[6] + pop[7]) / (2 * sum(pop))
  if(allele == 0) return(zeros)
  if(allele == 1) return(1 - zeros)
}
