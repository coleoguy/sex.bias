# parsing HD results
results.f <- as.data.frame(matrix(NA,0,11))
colnames(results.f) <- c("freqX", "freqY", "freqA", "model",
                         "com.sex", "com.num", "rare.num",
                         "s", "r", "h.mal.ben", "OSR")

load("~/Desktop/Dropbox/gitrepos/sex.bias/results/deprecated/hd.rare.males.RData")

counter <- 1
for(i in 1:length(results)){ # itterate through four levels of common sex
  for(j in 1:length(results$females1000)){
    for(k in 1:length(results$females1000$osr1)){
        cur.res <- results[[i]][[j]][[k]]
        com.num <- as.numeric(strsplit(names(results)[i], "s")[[1]][2])
        osr <-  as.numeric(strsplit(names(results[[i]])[j], "r")[[1]][2])
        rare.num <- com.num * osr
        s <- as.numeric(strsplit(names(results[[i]][[j]])[k], "s")[[1]][2])
        hvals <- c(0, .5, 1)
        for(n in 1:ncol(cur.res)){
          results.f[counter:(counter + nrow(cur.res) -1), 1] <- NA
          results.f[counter:(counter + nrow(cur.res) -1), 2] <- NA
          results.f[counter:(counter + nrow(cur.res) -1), 3] <- cur.res[, n]
          results.f[counter:(counter + nrow(cur.res) -1), 4] <- "HD"
          results.f[counter:(counter + nrow(cur.res) -1), 5] <- "female"
          results.f[counter:(counter + nrow(cur.res) -1), 6] <- com.num
          results.f[counter:(counter + nrow(cur.res) -1), 7] <- rare.num
          results.f[counter:(counter + nrow(cur.res) -1), 8] <- s
          results.f[counter:(counter + nrow(cur.res) -1), 9] <- NA
          results.f[counter:(counter + nrow(cur.res) -1), 10] <- hvals[n]
          results.f[counter:(counter + nrow(cur.res) -1), 11] <- osr
          counter <- counter + nrow(cur.res)
        }
      }
    }
  }
}












load("~/Desktop/Dropbox/gitrepos/sex.bias/results/deprecated/hd.rare.females.RData")

for(i in 1:length(results)){ # itterate through four levels of common sex
  for(j in 1:length(results$males1000)){
    for(k in 1:length(results$males1000$osr1)){
      cur.res <- results[[i]][[j]][[k]]
      com.num <- as.numeric(strsplit(names(results)[i], "s")[[1]][2])
      osr <-  as.numeric(strsplit(names(results[[i]])[j], "r")[[1]][2])
      rare.num <- com.num * osr
      s <- as.numeric(strsplit(names(results[[i]][[j]])[k], "s")[[1]][2])
      hvals <- c(0, .5, 1)
      for(n in 1:ncol(cur.res)){
        results.f[counter:(counter + nrow(cur.res) -1), 1] <- NA
        results.f[counter:(counter + nrow(cur.res) -1), 2] <- NA
        results.f[counter:(counter + nrow(cur.res) -1), 3] <- cur.res[, n]
        results.f[counter:(counter + nrow(cur.res) -1), 4] <- "HD"
        results.f[counter:(counter + nrow(cur.res) -1), 5] <- "male"
        results.f[counter:(counter + nrow(cur.res) -1), 6] <- com.num
        results.f[counter:(counter + nrow(cur.res) -1), 7] <- rare.num
        results.f[counter:(counter + nrow(cur.res) -1), 8] <- s
        results.f[counter:(counter + nrow(cur.res) -1), 9] <- NA
        results.f[counter:(counter + nrow(cur.res) -1), 10] <- hvals[n]
        results.f[counter:(counter + nrow(cur.res) -1), 11] <- osr
        counter <- counter + nrow(cur.res)
      }
    }
  }
}

library(ggraptR) 
ggraptR(results.f)
