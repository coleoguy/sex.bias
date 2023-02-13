library(ggplot2)
library(viridis)
library(ggraptR)

# set up matrix to store results and fill commsex, osr, Ne, and log.Ne columns
xchrom.comm.osr <- matrix(,28,5)
colnames(xchrom.comm.osr) <- c("freq", "commsex", "osr", "Ne", "log.Ne")
commsex <- rep(c(1000, 500, 100, 50), each = 7)
xchrom.comm.osr[, 2] <- commsex
osr <- rep(c(1, .8, .6, .4, .2, .1, .05))
xchrom.comm.osr[, 3] <- osr

Ne <- read_csv("Desktop/Dropbox/sex.bias/rscripts/results/Ne.allele.freqs/var.ne.xchrom.csv")[c(5:8), c(4:10)]

Ne.xchrom <- c()
for(i in 1:nrow(Ne)){
  Ne.xchrom <- c(Ne.xchrom, unlist(Ne[i,1:7]))
}
xchrom.comm.osr[,4] <- Ne.xchrom

log.Ne <- log10(Ne.xchrom)

xchrom.comm.osr <- cbind(xchrom.comm.osr, log.Ne)[,-c(5)]

# X chromosome allele frequencies
load("~/Desktop/Dropbox/sex.bias/rscripts/results/rare.male.model.RData")
res.rare.mal <- results

# first lets load the equal sex data
load("~/Desktop/Dropbox/sex.bias/rscripts/results/base.comp.model.RData")
res.base <- results

# lets remove the extra copy of the results
rm(results)

# Freq when there is no bias (uses the base.comp.model.RData file)
freq.additive.s05.osr1 <- res.base$pop1000$rd0.1$h0.5$s0.5[,1]
comm1000.osr1 <- mean(freq.additive.s05.osr1)
     
# First at 1000 comm sex
freq.additive.s05.osr.8 <- res.rare.mal$females1000$males0.8$rd0.1$h0.5$s0.5[,1]
comm1000.osr.8 <- mean(freq.additive.s05.osr.8)
     
freq.additive.s05.osr.6 <- res.rare.mal$females1000$males0.6$rd0.1$h0.5$s0.5[,1]
comm1000.osr.6 <- mean(freq.additive.s05.osr.6)
     
freq.additive.s05.osr.4 <- res.rare.mal$females1000$males0.4$rd0.1$h0.5$s0.5[,1]
comm1000.osr.4 <- mean(freq.additive.s05.osr.4)

freq.additive.s05.osr.2 <- res.rare.mal$females1000$males0.2$rd0.1$h0.5$s0.5[,1]
comm1000.osr.2 <- mean(freq.additive.s05.osr.2)

freq.additive.s05.osr.1 <- res.rare.mal$females1000$males0.1$rd0.1$h0.5$s0.5[,1]
comm1000.osr.1 <- mean(freq.additive.s05.osr.1)

freq.additive.s05.osr.05 <- res.rare.mal$females1000$males0.05$rd0.1$h0.5$s0.5[,1]
comm1000.osr.05 <- mean(freq.additive.s05.osr.05)
     
xchrom.comm1000 <- rbind(comm1000.osr1, comm1000.osr.8, comm1000.osr.6, comm1000.osr.4, comm1000.osr.2, 
                         comm1000.osr.1, comm1000.osr.05)

# Then 500

freq.additive.s05.osr1 <- res.base$pop500$rd0.1$h0.5$s0.5[,1]
comm500.osr1 <- mean(freq.additive.s05.osr1)

freq.additive.s05.osr.8 <- res.rare.mal$females500$males0.8$rd0.1$h0.5$s0.5[,1]
comm500.osr.8 <- mean(freq.additive.s05.osr.8)

freq.additive.s05.osr.6 <- res.rare.mal$females500$males0.6$rd0.1$h0.5$s0.5[,1]
comm500.osr.6 <- mean(freq.additive.s05.osr.6)

freq.additive.s05.osr.4 <- res.rare.mal$females500$males0.4$rd0.1$h0.5$s0.5[,1]
comm500.osr.4 <- mean(freq.additive.s05.osr.4)

freq.additive.s05.osr.2 <- res.rare.mal$females500$males0.2$rd0.1$h0.5$s0.5[,1]
comm500.osr.2 <- mean(freq.additive.s05.osr.2)

freq.additive.s05.osr.1 <- res.rare.mal$females500$males0.1$rd0.1$h0.5$s0.5[,1]
comm500.osr.1 <- mean(freq.additive.s05.osr.1)

freq.additive.s05.osr.05 <- res.rare.mal$females500$males0.05$rd0.1$h0.5$s0.5[,1]
comm500.osr.05 <- mean(freq.additive.s05.osr.05)

xchrom.comm500 <- rbind(comm500.osr1, comm500.osr.8, comm500.osr.6, comm500.osr.4, comm500.osr.2, 
                         comm500.osr.1, comm500.osr.05)

# Then 100
     
freq.additive.s05.osr1 <- res.base$pop100$rd0.1$h0.5$s0.5[,1]
comm100.osr1 <- mean(freq.additive.s05.osr1)

freq.additive.s05.osr.8 <- res.rare.mal$females100$males0.8$rd0.1$h0.5$s0.5[,1]
comm100.osr.8 <- mean(freq.additive.s05.osr.8)

freq.additive.s05.osr.6 <- res.rare.mal$females100$males0.6$rd0.1$h0.5$s0.5[,1]
comm100.osr.6 <- mean(freq.additive.s05.osr.6)

freq.additive.s05.osr.4 <- res.rare.mal$females100$males0.4$rd0.1$h0.5$s0.5[,1]
comm100.osr.4 <- mean(freq.additive.s05.osr.4)

freq.additive.s05.osr.2 <- res.rare.mal$females100$males0.2$rd0.1$h0.5$s0.5[,1]
comm100.osr.2 <- mean(freq.additive.s05.osr.2)

freq.additive.s05.osr.1 <- res.rare.mal$females100$males0.1$rd0.1$h0.5$s0.5[,1]
comm100.osr.1 <- mean(freq.additive.s05.osr.1)

freq.additive.s05.osr.05 <- res.rare.mal$females100$males0.05$rd0.1$h0.5$s0.5[,1]
comm100.osr.05 <- mean(freq.additive.s05.osr.05)

xchrom.comm100 <- rbind(comm100.osr1, comm100.osr.8, comm100.osr.6, comm100.osr.4, comm100.osr.2, 
                        comm100.osr.1, comm100.osr.05)

# Then 50

freq.additive.s05.osr1 <- res.base$pop50$rd0.1$h0.5$s0.5[,1]
comm50.osr1 <- mean(freq.additive.s05.osr1)

freq.additive.s05.osr.8 <- res.rare.mal$females50$males0.8$rd0.1$h0.5$s0.5[,1]
comm50.osr.8 <- mean(freq.additive.s05.osr.8)

freq.additive.s05.osr.6 <- res.rare.mal$females50$males0.6$rd0.1$h0.5$s0.5[,1]
comm50.osr.6 <- mean(freq.additive.s05.osr.6)

freq.additive.s05.osr.4 <- res.rare.mal$females50$males0.4$rd0.1$h0.5$s0.5[,1]
comm50.osr.4 <- mean(freq.additive.s05.osr.4)

freq.additive.s05.osr.2 <- res.rare.mal$females50$males0.2$rd0.1$h0.5$s0.5[,1]
comm50.osr.2 <- mean(freq.additive.s05.osr.2)

freq.additive.s05.osr.1 <- res.rare.mal$females50$males0.1$rd0.1$h0.5$s0.5[,1]
comm50.osr.1 <- mean(freq.additive.s05.osr.1)

freq.additive.s05.osr.05 <- res.rare.mal$females50$males0.05$rd0.1$h0.5$s0.5[,1]
comm50.osr.05 <- mean(freq.additive.s05.osr.05)

xchrom.comm50 <- rbind(comm50.osr1, comm50.osr.8, comm50.osr.6, comm50.osr.4, comm50.osr.2, 
                        comm50.osr.1, comm50.osr.05)
  
xchrom.freq <- rbind(xchrom.comm1000, xchrom.comm500, xchrom.comm100, xchrom.comm50)

xchrom.comm.osr[, 1] <- xchrom.freq
     
xchrom.freq <- as.data.frame(xchrom.comm.osr, stringsAsFactors = F)
  
ggplot(data = xchrom.freq, aes(x = log.Ne, y = as.numeric(xchrom.freq$freq))) +
       geom_point(aes(col=as.factor(osr)) , stat="identity", position="jitter", alpha=0.8, size=4) +
       scale_x_continuous(breaks =c(0,1,2,3), labels=c(c(1,10,100,1000))) +
       xlab("Variance effective population size") +
       ylab("Mean allele frequency") +
       guides(fill=guide_legend(title="OSR")) +
       scale_colour_viridis_d("OSR") +
       ggtitle("Mean frequency of the allele benefitting the common sex (females) on X chromosomes.  
Genetic architecture is additive. s= 0.5")

     
write.csv(xchrom.freq, file= "freq.xchrom.additive.s05.csv")
