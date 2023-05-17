# simulating populations with sex bias
# and sexual antagonism

# rd = .5
# males = 1000
# females = 1000
# h = 0
# s = .9
# 
# pop <- makeGenomes(males, females,
#                    freqs = c(rep(males/4, 4),
#                              rep(females/4, 4)))


Generation <- function(pop, females, males, rd, h, s){
  fit <- measureFit(h, s)
  parents <- GetParentsGeno(pop, fit, females, males)
  eggs <- makeEggs(parents$moms)
  sperm <- makeSperm(parents$dads, rd)
  pop <- makeNewPop(pop, eggs, sperm, females, males)
  return(pop)
}

measureFit <- function(h, s){
  if(h!=99){
    #hmal <- h
    #hfem <- 1-h
    fit <- c(1,   1+((1-h)*s),   1+((1-h)*s), (1+s),
             (1+s), (1+h*s),   (1+h*s),   1)
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
  dad.genomes <- sample(x.dad,
                        size = males,
                        replace = T,
                        prob=fit.dad)
  mom.geno <- rep(0, 4)
  names(mom.geno) <- c("X1X1", "X1X2", "X2X1", "X2X2")
  mom.geno[1] <- sum(mom.genomes==1)
  mom.geno[2] <- sum(mom.genomes==2)
  mom.geno[3] <- sum(mom.genomes==3)
  mom.geno[4] <- sum(mom.genomes==4)
  mom.geno
  
  dad.geno <- rep(0, 4)
  names(dad.geno) <- c("X1Y1", "X1Y2", "X2Y1", "X2Y2")
  dad.geno[1] <- sum(dad.genomes==5)
  dad.geno[2] <- sum(dad.genomes==6)
  dad.geno[3] <- sum(dad.genomes==7)
  dad.geno[4] <- sum(dad.genomes==8)
  dad.geno
  
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
  #make new females
  x  <-  paste(sample(1:2, size = females, replace=T, prob = eggs),
               sample(1:2, size = females, replace=T, prob = sperm[1:2]))
  pop[1] <- sum(x == "1 1")
  pop[2] <- sum(x == "1 2")
  pop[3] <- sum(x == "2 1")
  pop[4] <- sum(x == "2 2")
  #make new males
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
    if(allele == 2) twos <- (pop[6] + pop[8]) / males #Only males in calc bc only males can have Y and 1 in each male
  }
  if(chrom == "A"){
    ones <- (pop[1]*2 + pop[2] + pop[3] +
             pop[5]*2 + pop[6] + pop[7]) / (2 * sum(pop)) # *2 because 2 gametes per indiv
    twos <- 1 - ones
  }
  if(chrom == "X"){
    if(allele == 1){
      ones <- (pop[1] * 2 + pop[2] + pop[3] + pop[5] + pop[6]) /
        (males + females * 2) # mult by 2 because 2 gametes per indiv, both males and females can have 
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

getFits <- function(x){
  x <- as.numeric(x)
  rd <- x[6]
  s <- x[8]
  h <- x[7]
  A <- x[3]
  X <- x[1]
  Y <- x[2]
  # autosome
  if(rd == 0.5){
    # male benefit allele frequency
    p <- A #p is allele 1, 1 is beneficial for males
    q <- 1 - p #q is allele 2, 2 is beneficial for females
    # todo add fitness to below
    wMal <- ((p^2 * (1+s))  + (2*p*q * (1+(h*s))) + (q^2 * 1)) 
    
    wFem <- (p^2 * (1)) + (2*p*q * (1+((1-h)*s))) + (q^2 * (1+s))
    wdiff <- wMal - wFem #Male vs female
  }
  # sex chromosome
  if(rd != 0.5){
    # males # X Y
    p2 <- (1-X) * Y #p is allele 1, ben for males #x chr table data is for female ben. allele, subtract from 1 for male ben. allele freq
    q2 <- X * (1-Y) #q is allele 2, ben for females
    # todo add fitness to below
    
    #s <- x[8]
    #hmal <- 
    wMal <- ((p2 * (1+s)) + ((X*Y +  (1-X)*(1-Y)) * (1 + (h*s))) + (q2 * 1)) #
    
    # females # X X
    p <- X #p is allele 2, beneficial for females #x chr table data is for female ben. allele
    q <- 1 - p #q is allele 1, ben for males
    # todo add fitness to below
    
    wFem <- (p^2 * (1+s)) + (2*p*q * (1+((1-h)*s))) + (q^2 * 1)
    wdiff <- wMal - wFem #Male vs female
  }
  fits <- c(wMal, wFem, wdiff)
  names(fits) <- c("male","female", "wdiff")
  return(fits)
  
}


getProb <- function(x, y){
  #Taking each row/simulation that matches the scenario to extract how many 1s (fixed) and 0s (deleted) were for each scenario
  
  #y[,4] is equivalent to the commonsex column, labeled as males or females in input files
  temp <- y[y[,4] == x[i,4] & y$OSR == x[i,5] & y$rd == x[i,6] & y$h == x[i,7] & y$s == x[i,8],]
  
  #X chr fixation (allele 2 fixation, benefits female)
  Xfix <- sum(temp$X == 1)/length(temp$X)
  Xdel <- sum(temp$X == 0)/length(temp$X)

  #Y chromosome fixation (allele 1 fixation, benefits males)
  Yfit <- sum(temp$Y == 1)/length(temp$Y)
  Ydel <- sum(temp$Y == 0)/length(temp$Y)

  #Autosome fixation (allele 1 fixation, benefits male)
  Afit <- sum(temp$A == 1)/length(temp$A)
  Adel <- sum(temp$A == 0)/length(temp$A)

  fixes <- c(Xfix, Xdel, Yfit, Ydel, Afit, Adel)
  names(fixes) <- c("Xfix","Xdel","Yfix","Ydel","Afix","Adel")
  return(fixes)
}






