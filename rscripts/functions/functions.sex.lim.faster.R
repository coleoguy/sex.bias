# Heath Blackmon
# coleoguy@gmail.com
# simulating populations with sex bias
# and deleterious alleles
# 1/1000 alleles mutates to sex deleterious


males = 100
females = 20
s = .3
iter = 50
gens = 100
cores = 14

runSim <- function(males = 100, females = 20, s = .3,
                   iter = 50000, gens = 1000, cores = 14){
  registerDoMC(14)
  x <- foreach(n=1:iter, .combine = "rbind") %dopar%{
    resultA <- c()
    pop <- makeGenomes(females = females,
                       males = males,
                       freqs = c(c(females, 0, 0, 0),
                                 c(males, 0, 0, 0)))
    for(k in 1:gens){
      resultA[k] <- GetFreq(pop, chrom = "A", allele = 2,
                            females = females, males = males)
      pop <- Generation(pop, females = females, males = males,
                        rd = .5, h = .5, s = s)
    }
    # this puts the final result of our sim together
    resultA
  }
  return(x)
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
measureFit <- function(h, s){
  # mutation is deleterious in females
  fit <- c(1,   1-h*s, 1-h*s,   1-s,
           1,   1,     1,       1)
  return(fit)
}
GetParentsGeno <- function(pop, fit, females, males){
  # this effectively performs viability selection
  # so we sample females as moms based on fitness
  mom.genomes <- sample(x = 1:4,
                        size = females,
                        replace = T,
                        prob=pop[1:4]*fit[1:4])
  dad.genomes <- sample(5:8,
                        size = males,
                        replace = T,
                        prob=pop[5:8]*fit[5:8])
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
Generation <- function(pop, females, males, rd, h, s){
  fit <- measureFit(h, s)
  parents <- GetParentsGeno(pop, fit, females, males)
  eggs <- makeEggs(parents$moms)
  sperm <- makeSperm(parents$dads, rd)
  pop <- makeNewPop(pop, eggs, sperm, females, males)
  muts <- sum(rbinom(2*sum(pop), 1, .0001))
  if(muts > 0){
    z <- which(pop > 0)
    z <- z[!z %in% c(4, 8)]
    mut.possible <- c()
    # possible genotypes to mutate
    if(length(z)>0){
      for(q in z){
        mut.possible <- c(rep(q, pop[q]), mut.possible)
      }
      # check to make sure it isnt already fixed
      # draw genotypes to mutate
      hits <- sample(mut.possible, muts, replace=T) # DID THIS FIX IT
      # remove ind from old genotype
      pop[1] <- pop[1] - sum(hits==1)
      pop[2] <- pop[2] - sum(hits==2)
      pop[3] <- pop[3] - sum(hits==3)
      pop[5] <- pop[5] - sum(hits==5)
      pop[6] <- pop[6] - sum(hits==6)
      pop[7] <- pop[7] - sum(hits==7)
      # add to new genotypes
      pop[2] <- pop[2] + sum(hits==1)
      pop[4] <- pop[4] + sum(hits%in%c(2,3))
      pop[6] <- pop[6] + sum(hits==5)
      pop[8] <- pop[8] + sum(hits%in%c(6,7))
    }
  }
  #if(min(pop)<0) print("AHHH")
  pop[pop<0]<-0
  return(pop)
}

