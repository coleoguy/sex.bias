# Julio Rincones Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com
# P.I. Dr. Heath Blackmon
# coleoguy@gmail.com
# simulating populations with sex bias
# and sexual antagonism

makeGenomes <- function(females, males, freqs=NULL){
  population <- rep(0, 8)
  if(!is.null(freqs)){
    population <- freqs
  }else{
    print("supply frequencies")
  }
  names(population) <- c("fem.X1X1", "fem.X1X2", "fem.X2X1", "fem.X2X2",
                         "mal.X1Y1", "mal.X1Y2", "mal.X2Y1", "mal.X2Y2")
  return(population)
}


measureFit <- function(pop, h, s){
  if(h!=99){
    fit <- c(1,   1+h*s,   1+h*s,   1+s,
             1, 1/(1+h*s), 1/(1+h*s), 1/(1+s))
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
  # THIS IS WHERE THE PROBLEM IS WE ARE SAMPLING BASED ON FITNESS BUT EQUAL 
  # NUMBERS OF EACH then we sample males as dads based on fitness
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
  sperm <- rep(0, 4)
  names(sperm) <- c("X1", "X2", "Y1", "Y2")
  sperm[1] <- 0.5 * dad.geno[1] +
    0.5 * dad.geno[2] * (1 - rd)  +
    0.5 * dad.geno[3] * rd
  sperm[2] <- 0.5 * dad.geno[4] + 
    0.5 * dad.geno[2] * rd +
    0.5 * dad.geno[3] * (1 - rd)
  sperm[3] <- 0.5 * dad.geno[1] + 
    0.5 * dad.geno[2] * rd +
    0.5 * dad.geno[3] * (1 - rd)
  sperm[4] <- 0.5 * dad.geno[4] + 
    0.5 * dad.geno[2] * (1 - rd)  +
    0.5 * dad.geno[3] * rd
  return(sperm)
}

makeNewPop <- function(pop, eggs, sperm, females, males){
  x  <-  paste(sample(1:2, size = females, replace=T, prob = eggs),
               sample(1:2, size = females, replace=T, prob = sperm[1:2]))
  pop[1] <- sum(x == "1 1")
  pop[2] <- sum(x == "1 2")
  pop[3] <- sum(x == "2 1")
  pop[4] <- sum(x == "2 2")
  x <- paste(sample(1:2, size = males, replace=T, prob = eggs),
             sample(1:2, size = males, replace=T, prob = sperm[3:4]))
  pop[5] <- sum(x == "1 1")
  pop[6] <- sum(x == "1 2")
  pop[7] <- sum(x == "2 1")
  pop[8] <- sum(x == "2 2")
  return(pop)
}

Generation <- function(pop, females, males, rd, h, s){
  fit <- measureFit(pop, h, s)
  parents <- GetParentsGeno(pop, fit, females, males)
  eggs <- makeEggs(parents$moms)
  sperm <- makeSperm(parents$dads, rd)
  pop <- makeNewPop(pop, eggs, sperm, females, males)
  return(pop)
}

GetFreq <- function(pop, chrom, allele, males, females){
  ones <- twos <- 0
  if(chrom == "Y"){
    if(allele == 1) ones <- (pop[5] + pop[7]) / males
    if(allele == 2) twos <- (pop[6] + pop[8]) / males
  }
  if(chrom == "A"){
    ones <- (pop[1]*2 + pop[2] + pop[3] +
             pop[5]*2 + pop[6] + pop[7]) / (2 * sum(pop))
    twos <- 1 - ones
  }
  if(chrom == "X"){
    if(allele == 1){
      ones <- (pop[1] * 2 + pop[2] + pop[3] + pop[5] + pop[6]) / 
        (males + females * 2)
    }
    if(allele == 2){
      twos <- (pop[4] * 2 + pop[2] + pop[3] + pop[7] + pop[8]) / 
        (males + females * 2)
    }
  }
  if(allele == 1) return(ones)
  if(allele == 2) return(twos)
}
