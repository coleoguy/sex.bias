# Julio Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com
# Advisor: Heath Blackmon
# coleoguy@gmail.com
# Script to parse XY model data. The aim being to store all the allele frequencies 
# in long format along with the common sex, osr, h, rd, and corresponding selection 
# coefficient values associated with each frequency.

load("../results/base.comp.model.RData")
no.bias <- results
rm(results)

# allele 1 is fit for males
# allele 2 is fit for females

# First let's parse the instances where OSR=1 for the autosomal loci
# structure of the data = no.bias[["pop50"]][["rd0.1"]][["h0"]][["s0.1"]]
# Column 1 = Freq of A2 on the X chromosome
# Column 2 = Freq of A1 on the Y chromosome
# Column 3 = Freq of A1 on Autosomes

result.no.bias <- as.data.frame(matrix(,0,5))
colnames(result.no.bias) <- c("freq2", "common.sex", "rd", "h", "s")
com.nums <- c(50,100,500,1000)
rd.levels <- c(0,0.1,0.2,0.5)
h.factor <- c("h0", "h0.5", "h1", "h99")
s.coeff <- c("s0.1", "s0.2", "s0.5", "s0.9")
osr <- 1

i <- j <- k <- m <- 4
for(i in 1:length(no.bias)){ # iterate through popsize
  print(i)
  for(j in 1:length(no.bias[[i]])){  #iterate through rd
    for(k in 1:length(no.bias[[i]][[j]])){ # iterate through h
      for(m in 1:length(no.bias[[i]][[j]][[k]])){ # iterate through s
        # Results of these loops saved in vector freq1 
        freq2 <- no.bias[[i]][[j]][[k]][[m]][,3] # store the frequency 
        comm <- rep(com.nums[i], 500)
        rd <- rep(rd.levels[j], 500)
        h <- rep(h.factor[k], 500)
        s <- rep(s.coeff[m], 500)
        foo <- cbind(freq2,comm, osr, rd,h,s)
        result.no.bias <- rbind(foo, result.no.bias)
      }
    }
  }
}


# Then parse when there is OSR bias
load("../results/rare.male.model.RData")
rare.mal <- results
rm(results)

result <- as.data.frame(matrix(,0,6))
colnames(result) <- c("freq2", "comm", "osr", "rd", "h", "s")
com.nums <- c(50,100,500,1000)
osr.level <- c(.8,.6,.4,.2,.1,.05)
rd.levels <- c(0,0.1,0.2,0.5)

for(i in 1:length(rare.mal)){ # iterate through popsize
  print(i)
  for(j in 1:length(rare.mal[[i]])){  #iterate through osr
    for(k in 1:length(rare.mal[[i]][[j]])){ # iterate through rd
      for(m in 1:length(rare.mal[[i]][[j]][[k]])){ # iterate through h
        for(n in 1:length(rare.mal[[i]][[j]][[k]][[m]])){# iterate through s
          # Results of these loops saved in vector x 
          freq2 <- rare.mal[[i]][[j]][[k]][[m]][[n]][,3] # store the frequency 
          comm <- rep(com.nums[i], 1000)
          osr <- rep(osr.level[j], 1000)
          rd <- rep(rd.levels[k],1000)
          h <- rep(names(rare.mal[[i]][[j]][[k]])[m], 1000)
          s <- rep(names(rare.mal[[i]][[j]][[k]][[m]])[n], 1000)
          foo <- cbind(freq2,comm,osr,rd,h,s)
          result <- rbind(foo, result)
        }
      }
    }
  }
}

autosomal <- rbind(result, result.no.bias)

write.csv(autosomal, file="Autosome.rare.male.csv")

rm(list = ls())

### Now parse the frequency of the A1  allele on the Y chromosome (when females are the rare sex).

# for rare males(i.e females are the common sex) 
# these are the frequencies that were saved with  each iteration
# resultA[p] <- GetFreq(pop, chrom="A", allele = 1, females=females[i], males=males)
# i.e The frequency  of the male beneficial allele on the Autosome
# resultY[p] <- GetFreq(pop, chrom="Y", allele = 1, females=females[i], males=males) 
# i.e The frequency of the male beneficial allele on the Y chromosome
# resultX[p] <- GetFreq(pop, chrom="X", allele = 2, females=females[i], males=males)
# i.e.The frequency  of the female beneficial allele on the X chromosome.

# column 1 of the results corresponds to the X chromosome
# column 2 to the Y chromosome
# column 3 to the autosome.

##### both simulations for rare males or rare females follow the same scheme

load("../results/base.comp.model.RData")
no.bias <- results
rm(results)

# allele 1 is fit for males
# allele 2 is fit for females

# First let's parse the instances where OSR=1 for the autosomal loci
# structure of the data = no.bias[["pop50"]][["rd0.1"]][["h0"]][["s0.1"]]
# Column 1 = Freq of A2 on the X chromosome
# Column 2 = Freq of A1 on the Y chromosome
# Column 3 = Freq of A1 on Autosomes

result.no.bias <- as.data.frame(matrix(,0,5))
colnames(result.no.bias) <- c("freq1", "common.sex", "rd", "h", "s")
com.nums <- c(50,100,500,1000)
rd.levels <- c(0,0.1,0.2,0.5)
h.factor <- c("h0", "h0.5", "h1", "h99")
s.coeff <- c("s0.1", "s0.2", "s0.5", "s0.9")
osr <- 1

i <- j <- k <- m <- 4
for(i in 1:length(no.bias)){ # iterate through popsize
  print(i)
  for(j in 1:length(no.bias[[i]])){  #iterate through rd
    for(k in 1:length(no.bias[[i]][[j]])){ # iterate through h
      for(m in 1:length(no.bias[[i]][[j]][[k]])){ # iterate through s
        # Results of these loops saved in vector freq1 
        freq1 <- no.bias[[i]][[j]][[k]][[m]][,2] # store the frequency 
        comm <- rep(com.nums[i], 500)
        rd <- rep(rd.levels[j], 500)
        h <- rep(h.factor[k], 500)
        s <- rep(s.coeff[m], 500)
        foo <- cbind(freq1,comm, osr, rd,h,s)
        result.no.bias <- rbind(foo, result.no.bias)
      }
    }
  }
}


load("../results/rare.female.model.RData")
rare.fem <- results
rm(results)

# Convert the list of results above as data frame and assign column names
result <- as.data.frame(matrix(,0,6))
colnames(result) <- c("freq1", "comm", "osr", "rd", "h", "s")
com.nums <- c(50,100,500,1000)
osr.level <- c(.8,.6,.4,.2,.1,.05)
rd.levels <- c(0,0.1,0.2,0.5)
# nested for loops to iterate through the results
# cols are X,  y, a
for(i in 1:length(rare.fem)){ # iterate through popsize
  print(i)
  for(j in 1:length(rare.fem[[i]])){  #iterate through osr
    for(k in 1:length(rare.fem[[i]][[j]])){ # iterate through rd
      for(m in 1:length(rare.fem[[i]][[j]][[k]])){ # iterate through h
        for(n in 1:length(rare.fem[[i]][[j]][[k]][[m]])){# iterate through s
          # Results of these loops saved in vector x 
          freq1 <- rare.fem[[i]][[j]][[k]][[m]][[n]][,2] # store the frequency 
          comm <- rep(com.nums[i], 1000)
          osr <- rep(osr.level[j], 1000)
          rd <- rep(rd.levels[k],1000)
          h <- rep(names(rare.fem[[i]][[j]][[k]])[m], 1000)
          s <- rep(names(rare.fem[[i]][[j]][[k]][[m]])[n], 1000)
          foo <- cbind(freq1,comm,osr,rd,h,s)
          result <- rbind(foo, result)
        }
      }
    }
  }
}

ychrom <- rbind(result, result.no.bias)

write.csv(ychrom, file="Ychrom.rare.female.csv")

rm(list = ls())

### Parse through the results for the A2 allele on the X chromosome (when males are the rare sex)

load("../results/base.comp.model.RData")
no.bias <- results
rm(results)

# allele 1 is fit for males
# allele 2 is fit for females

# First let's parse the instances where OSR=1 for the autosomal loci
# structure of the data = no.bias[["pop50"]][["rd0.1"]][["h0"]][["s0.1"]]
# Column 1 = Freq of A2 on the X chromosome
# Column 2 = Freq of A1 on the Y chromosome
# Column 3 = Freq of A1 on Autosomes

result.no.bias <- as.data.frame(matrix(,0,5))
colnames(result.no.bias) <- c("freq2", "common.sex", "rd", "h", "s")
com.nums <- c(50,100,500,1000)
rd.levels <- c(0,0.1,0.2,0.5)
h.factor <- c("h0", "h0.5", "h1", "h99")
s.coeff <- c("s0.1", "s0.2", "s0.5", "s0.9")
osr <- 1

i <- j <- k <- m <- 4
for(i in 1:length(no.bias)){ # iterate through popsize
  print(i)
  for(j in 1:length(no.bias[[i]])){  #iterate through rd
    for(k in 1:length(no.bias[[i]][[j]])){ # iterate through h
      for(m in 1:length(no.bias[[i]][[j]][[k]])){ # iterate through s
        # Results of these loops saved in vector freq1 
        freq2 <- no.bias[[i]][[j]][[k]][[m]][,1] # store the frequency 
        comm <- rep(com.nums[i], 500)
        rd <- rep(rd.levels[j], 500)
        h <- rep(h.factor[k], 500)
        s <- rep(s.coeff[m], 500)
        foo <- cbind(freq2, comm, osr, rd,h,s)
        result.no.bias <- rbind(foo, result.no.bias)
      }
    }
  }
}


load("../results/rare.male.model.RData")
rare.mal <- results
rm(results)


result <- as.data.frame(matrix(,0,6))
colnames(result) <- c("freq2", "comm", "osr", "rd", "h", "s")
com.nums <- c(50,100,500,1000)
osr.level <- c(.8,.6,.4,.2,.1,.05)
rd.levels <- c(0,0.1,0.2,0.5)
# nested for loops to iterate through the results
# cols are X,  Y, A
for(i in 1:length(rare.mal)){ # iterate through popsize
  print(i)
  for(j in 1:length(rare.mal[[i]])){  #iterate through osr
    for(k in 1:length(rare.mal[[i]][[j]])){ # iterate through rd
      for(m in 1:length(rare.mal[[i]][[j]][[k]])){ # iterate through h
        for(n in 1:length(rare.mal[[i]][[j]][[k]][[m]])){# iterate through s
          # Results of these loops saved in vector x 
          freq2 <- rare.mal[[i]][[j]][[k]][[m]][[n]][,1] # store the frequency 
          comm <- rep(com.nums[i], 1000)
          osr <- rep(osr.level[j], 1000)
          rd <- rep(rd.levels[k],1000)
          h <- rep(names(rare.mal[[i]][[j]][[k]])[m], 1000)
          s <- rep(names(rare.mal[[i]][[j]][[k]][[m]])[n], 1000)
          foo <- cbind(freq2,comm,osr,rd,h,s)
          result <- rbind(foo, result)
        }
      }
    }
  }
}

xchrom <- rbind(result, result.no.bias)

write.csv(xchrom,file="Xchrom.rare.male.csv")

rm(list = ls())


