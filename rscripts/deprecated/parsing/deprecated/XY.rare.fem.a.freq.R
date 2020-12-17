# Julio Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com
# Advisor: Heath Blackmon
# coleoguy@gmail.com
# Script to parse XY model data. The aim being to store all the allele frequencies 
# in long format along with the common sex, osr, h, rd, and corresponding selection 
# coefficient values associated with each frequency.

load("../results/rare.female.model.RData")
rare.fem <- results
rm(results)

# allele 1 is fit for males
# allele 2 is fit for females

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


write.csv(result, file="Ychrom.rare.female.csv")

rm(list = ls())

load("../results/rare.female.model.RData")
rare.fem <- results
rm(results)

result <- as.data.frame(matrix(,0,6))
colnames(result) <- c("freq1", "comm", "osr", "rd", "h", "s")
com.nums <- c(50,100,500,1000)
osr.level <- c(.8,.6,.4,.2,.1,.05)
rd.levels <- c(0,0.1,0.2,0.5)

result <- as.data.frame(matrix(,0,6))
colnames(result) <- c("freq1", "comm", "osr", "rd", "h", "s")
com.nums <- c(50,100,500,1000)
osr.level <- c(.8,.6,.4,.2,.1,.05)
rd.levels <- c(0,0.1,0.2,0.5)
# nested for loops to iterate through the results
# cols are X,  Y, A
for(i in 1:length(rare.fem)){ # iterate through popsize
  print(i)
  for(j in 1:length(rare.fem[[i]])){  #iterate through osr
    for(k in 1:length(rare.fem[[i]][[j]])){ # iterate through rd
      for(m in 1:length(rare.fem[[i]][[j]][[k]])){ # iterate through h
        for(n in 1:length(rare.fem[[i]][[j]][[k]][[m]])){# iterate through s
          # Results of these loops saved in vector x 
          freq1 <- rare.fem[[i]][[j]][[k]][[m]][[n]][,1] # store the frequency 
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

write.csv(result, file="Xchrom.rare.female.csv")

rm(list = ls())

load("..results/rare.female.model.RData")
rare.fem <- results
rm(results)

result <- as.data.frame(matrix(,0,6))
colnames(result) <- c("freq1", "comm", "osr", "rd", "h", "s")
com.nums <- c(50,100,500,1000)
osr.level <- c(.8,.6,.4,.2,.1,.05)
rd.levels <- c(0,0.1,0.2,0.5)

for(i in 1:length(rare.fem)){ # iterate through popsize
  print(i)
  for(j in 1:length(rare.fem[[i]])){  #iterate through osr
    for(k in 1:length(rare.fem[[i]][[j]])){ # iterate through rd
      for(m in 1:length(rare.fem[[i]][[j]][[k]])){ # iterate through h
        for(n in 1:length(rare.fem[[i]][[j]][[k]][[m]])){# iterate through s
          # Results of these loops saved in vector x 
          freq1 <- rare.fem[[i]][[j]][[k]][[m]][[n]][,3] # store the frequency 
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

write.csv(result, file="Autosome.rare.female.csv")
