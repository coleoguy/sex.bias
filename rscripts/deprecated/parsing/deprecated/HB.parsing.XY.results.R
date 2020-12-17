# This script build single R object to hold everything

load("~/Desktop/Dropbox/gitrepos/sex.bias/results/rare.male.model.RData")
final.results <- matrix(NA, 3000000, 11)
colnames(final.results) <- c("freqX", "freqY", "freqA", "model", "com.sex",
                       "com.num", "rare.num", "s",
                       "r", "h.mal.ben","OSR")
final.results <- as.data.frame(final.results)


# col 1 is X
# col 2 is Y
# col 3 is A
counter <- 1
for(i in 1:length(results)){
  print(counter)
  for(j in 1:length(results$females50)){
    for(k in 1:length(results$females50$males0.8)){
      for(m in 1:length(results$females50$males0.8$rd0)){
        for(n in 1:length(results$females50$males0.8$rd0$h0)){
          foo <- results[[i]][[j]][[k]][[m]][[n]]
          # pull autosome data
          if(names(results[[i]][[j]])[k] == "rd0.5"){
            comnum <- as.numeric(strsplit(names(results)[i],"s")[[1]][2])
            osr <- as.numeric(strsplit(names(results[[i]])[j], "s")[[1]][2])
            rar.num <- osr*comnum
            rd <- as.numeric(strsplit(names(results[[i]][[j]])[k], "d")[[1]][2])
            h <- as.numeric(strsplit(names(results[[i]][[j]][[k]])[m], "h")[[1]][2])
            s <- as.numeric(strsplit(names(results[[i]][[j]][[k]][[m]])[n], "s")[[1]][2])
            final.results[counter:(counter + nrow(foo) - 1), 1] <- foo[,1]
            final.results[counter:(counter + nrow(foo) - 1), 2] <- foo[,2]
            final.results[counter:(counter + nrow(foo) - 1), 3] <- foo[,3]
            final.results[counter:(counter + nrow(foo) - 1), 4] <- "XY"
            final.results[counter:(counter + nrow(foo) - 1), 5] <- "females"
            final.results[counter:(counter + nrow(foo) - 1), 6] <- comnum
            final.results[counter:(counter + nrow(foo) - 1), 7] <- rar.num
            final.results[counter:(counter + nrow(foo) - 1), 8] <- s
            final.results[counter:(counter + nrow(foo) - 1), 9] <- rd
            final.results[counter:(counter + nrow(foo) - 1), 10] <- h
            final.results[counter:(counter + nrow(foo) - 1), 11] <- osr
            counter <- counter + nrow(foo)
          }
          # pull sex chrom data
          if(names(results[[i]][[j]])[k] != "rd0.5"){
            comnum <- as.numeric(strsplit(names(results)[i],"s")[[1]][2])
            osr <- as.numeric(strsplit(names(results[[i]])[j], "s")[[1]][2])
            rar.num <- osr*comnum
            rd <- as.numeric(strsplit(names(results[[i]][[j]])[k], "d")[[1]][2])
            h <- as.numeric(strsplit(names(results[[i]][[j]][[k]])[m], "h")[[1]][2])
            s <- as.numeric(strsplit(names(results[[i]][[j]][[k]][[m]])[n], "s")[[1]][2])
            final.results[counter:(counter + nrow(foo) - 1), 1] <- foo[,1]
            final.results[counter:(counter + nrow(foo) - 1), 2] <- foo[,2]
            final.results[counter:(counter + nrow(foo) - 1), 3] <- foo[,3]
            final.results[counter:(counter + nrow(foo) - 1), 4] <- "XY"
            final.results[counter:(counter + nrow(foo) - 1), 5] <- "females"
            final.results[counter:(counter + nrow(foo) - 1), 6] <- comnum
            final.results[counter:(counter + nrow(foo) - 1), 7] <- rar.num
            final.results[counter:(counter + nrow(foo) - 1), 8] <- s
            final.results[counter:(counter + nrow(foo) - 1), 9] <- rd
            final.results[counter:(counter + nrow(foo) - 1), 10] <- h
            final.results[counter:(counter + nrow(foo) - 1), 11] <- osr
            counter <- counter + nrow(foo)
          }
        }
      }
    }
  }
}
rm(list=ls()[!ls() %in% c("counter", "final.results")])

load("~/Desktop/Dropbox/gitrepos/sex.bias/results/rare.female.model.RData")


# col 1 is X
# col 2 is Y
# col 3 is A
for(i in 1:length(results)){
  print(counter)
  for(j in 1:length(results$males50)){
    for(k in 1:length(results$males50$females0.8)){
      for(m in 1:length(results$males50$females0.8$rd0)){
        for(n in 1:length(results$males50$females0.8$rd0$h0)){
          foo <- results[[i]][[j]][[k]][[m]][[n]]
          # autosomes are same for both sexes so we can skip rd0.5
          # pull sex chrom data
          if(names(results[[i]][[j]])[k] != "rd0.5"){
            comnum <- as.numeric(strsplit(names(results)[i],"s")[[1]][2])
            osr <- as.numeric(strsplit(names(results[[i]])[j], "s")[[1]][2])
            rar.num <- osr*comnum
            rd <- as.numeric(strsplit(names(results[[i]][[j]])[k], "d")[[1]][2])
            h <- as.numeric(strsplit(names(results[[i]][[j]][[k]])[m], "h")[[1]][2])
            s <- as.numeric(strsplit(names(results[[i]][[j]][[k]][[m]])[n], "s")[[1]][2])
            final.results[counter:(counter + nrow(foo) - 1), 1] <- foo[,1]
            final.results[counter:(counter + nrow(foo) - 1), 2] <- foo[,2]
            final.results[counter:(counter + nrow(foo) - 1), 3] <- foo[,3]
            final.results[counter:(counter + nrow(foo) - 1), 4] <- "XY"
            final.results[counter:(counter + nrow(foo) - 1), 5] <- "males"
            final.results[counter:(counter + nrow(foo) - 1), 6] <- comnum
            final.results[counter:(counter + nrow(foo) - 1), 7] <- rar.num
            final.results[counter:(counter + nrow(foo) - 1), 8] <- s
            final.results[counter:(counter + nrow(foo) - 1), 9] <- rd
            final.results[counter:(counter + nrow(foo) - 1), 10] <- h
            final.results[counter:(counter + nrow(foo) - 1), 11] <- osr
            counter <- counter + nrow(foo)
          }
        }
      }
    }
  }
}

rm(list=ls()[!ls() %in% c("counter", "final.results")])



load("~/Desktop/Dropbox/gitrepos/sex.bias/results/base.comp.model.RData")

# col 1 is X
# col 2 is Y
# col 3 is A
for(i in 1:length(results)){
  print(counter)
  for(j in 1:length(results$pop50)){
    for(k in 1:length(results$pop50$rd0)){
      for(m in 1:length(results$pop50$rd0$h0)){
          foo <- results[[i]][[j]][[k]][[m]]
          # autosomes are same for both sexes so we can skip rd0.5
          # pull sex chrom data
            comnum <- rar.num <- as.numeric(strsplit(names(results)[i],"p")[[1]][3])
            rd <- as.numeric(strsplit(names(results[[i]])[j], "d")[[1]][2])
            h <- as.numeric(strsplit(names(results[[i]][[j]])[k], "h")[[1]][2])
            s <- as.numeric(strsplit(names(results[[i]][[j]][[k]])[m], "s")[[1]][2])
            final.results[counter:(counter + nrow(foo) - 1), 1] <- foo[,1]
            final.results[counter:(counter + nrow(foo) - 1), 2] <- foo[,2]
            final.results[counter:(counter + nrow(foo) - 1), 3] <- foo[,3]
            final.results[counter:(counter + nrow(foo) - 1), 4] <- "XY"
            final.results[counter:(counter + nrow(foo) - 1), 5] <- "equal"
            final.results[counter:(counter + nrow(foo) - 1), 6] <- comnum
            final.results[counter:(counter + nrow(foo) - 1), 7] <- rar.num
            final.results[counter:(counter + nrow(foo) - 1), 8] <- s
            final.results[counter:(counter + nrow(foo) - 1), 9] <- rd
            final.results[counter:(counter + nrow(foo) - 1), 10] <- h
            final.results[counter:(counter + nrow(foo) - 1), 11] <- 1
            counter <- counter + nrow(foo)
          }
        }
      }
    }

rm(list=ls()[!ls() %in% c("counter", "final.results")])

results <- final.results[1:(counter-1),]

library(ggraptR)
ggraptR(final.results)
