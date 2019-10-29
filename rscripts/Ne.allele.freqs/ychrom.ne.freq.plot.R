library(ggplot2)
library(viridis)
library(ggraptR)

# set up matrix to store results and fill commsex, osr, Ne, and log.Ne columns
ychrom.comm.osr <- matrix(,28,5)
colnames(ychrom.comm.osr) <- c("freq", "commsex", "osr", "Ne", "log.Ne")
commsex <- rep(c(1000, 500, 100, 50), each = 7)
ychrom.comm.osr[, 2] <- commsex
osr <- rep(c(1, .8, .6, .4, .2, .1, .05))
ychrom.comm.osr[, 3] <- osr

Ne <- rep(c(1000, 500, 100, 50), each = 7)
ychrom.comm.osr[, 4] <- Ne

log.Ne <- log10(Ne)

ychrom.comm.osr <- cbind(ychrom.comm.osr, log.Ne)[,-c(5)]

# X chromosome allele frequencies
load("../results/rare.female.model.RData")
res.rare.fem <- results

# first lets load the equal sex data
load("..results/base.comp.model.RData")
res.base <- results

# lets remove the extra copy of the results
rm(results)

# Freq when there is no bias (uses the base.comp.model.RData file)
freq.additive.s05.osr1 <- res.base$pop1000$rd0.1$h0.5$s0.5[,2]
comm1000.osr1 <- mean(freq.additive.s05.osr1)

# First at 1000 comm sex
freq.additive.s05.osr.8 <- res.rare.fem$males1000$females0.8$rd0.1$h0.5$s0.5[,2]
comm1000.osr.8 <- mean(freq.additive.s05.osr.8)

freq.additive.s05.osr.6 <- res.rare.fem$males1000$females0.6$rd0.1$h0.5$s0.5[,2]
comm1000.osr.6 <- mean(freq.additive.s05.osr.6)

freq.additive.s05.osr.4 <- res.rare.fem$males1000$females0.4$rd0.1$h0.5$s0.5[,2]
comm1000.osr.4 <- mean(freq.additive.s05.osr.4)

freq.additive.s05.osr.2 <- res.rare.fem$males1000$females0.2$rd0.1$h0.5$s0.5[,2]
comm1000.osr.2 <- mean(freq.additive.s05.osr.2)

freq.additive.s05.osr.1 <- res.rare.fem$males1000$females0.1$rd0.1$h0.5$s0.5[,2]
comm1000.osr.1 <- mean(freq.additive.s05.osr.1)

freq.additive.s05.osr.05 <- res.rare.fem$males1000$females0.05$rd0.1$h0.5$s0.5[,2]
comm1000.osr.05 <- mean(freq.additive.s05.osr.05)

ychrom.comm1000 <- rbind(comm1000.osr1, comm1000.osr.8, comm1000.osr.6, comm1000.osr.4, comm1000.osr.2, 
                         comm1000.osr.1, comm1000.osr.05)

# Then 500

freq.additive.s05.osr1 <- res.base$pop500$rd0.1$h0.5$s0.5[,2]
comm500.osr1 <- mean(freq.additive.s05.osr1)

freq.additive.s05.osr.8 <- res.rare.fem$males500$females0.8$rd0.1$h0.5$s0.5[,2]
comm500.osr.8 <- mean(freq.additive.s05.osr.8)

freq.additive.s05.osr.6 <- res.rare.fem$males500$females0.6$rd0.1$h0.5$s0.5[,2]
comm500.osr.6 <- mean(freq.additive.s05.osr.6)

freq.additive.s05.osr.4 <- res.rare.fem$males500$females0.4$rd0.1$h0.5$s0.5[,2]
comm500.osr.4 <- mean(freq.additive.s05.osr.4)

freq.additive.s05.osr.2 <- res.rare.fem$males500$females0.2$rd0.1$h0.5$s0.5[,2]
comm500.osr.2 <- mean(freq.additive.s05.osr.2)

freq.additive.s05.osr.1 <- res.rare.fem$males500$females0.1$rd0.1$h0.5$s0.5[,2]
comm500.osr.1 <- mean(freq.additive.s05.osr.1)

freq.additive.s05.osr.05 <- res.rare.fem$males500$females0.05$rd0.1$h0.5$s0.5[,2]
comm500.osr.05 <- mean(freq.additive.s05.osr.05)

ychrom.comm500 <- rbind(comm500.osr1, comm500.osr.8, comm500.osr.6, comm500.osr.4, comm500.osr.2, 
                        comm500.osr.1, comm500.osr.05)

# Then 100

freq.additive.s05.osr1 <- res.base$pop100$rd0.1$h0.5$s0.5[,2]
comm100.osr1 <- mean(freq.additive.s05.osr1)

freq.additive.s05.osr.8 <- res.rare.fem$males100$females0.8$rd0.1$h0.5$s0.5[,2]
comm100.osr.8 <- mean(freq.additive.s05.osr.8)

freq.additive.s05.osr.6 <- res.rare.fem$males100$females0.6$rd0.1$h0.5$s0.5[,2]
comm100.osr.6 <- mean(freq.additive.s05.osr.6)

freq.additive.s05.osr.4 <- res.rare.fem$males100$females0.4$rd0.1$h0.5$s0.5[,2]
comm100.osr.4 <- mean(freq.additive.s05.osr.4)

freq.additive.s05.osr.2 <- res.rare.fem$males100$females0.2$rd0.1$h0.5$s0.5[,2]
comm100.osr.2 <- mean(freq.additive.s05.osr.2)

freq.additive.s05.osr.1 <- res.rare.fem$males100$females0.1$rd0.1$h0.5$s0.5[,2]
comm100.osr.1 <- mean(freq.additive.s05.osr.1)

freq.additive.s05.osr.05 <- res.rare.fem$males100$females0.05$rd0.1$h0.5$s0.5[,2]
comm100.osr.05 <- mean(freq.additive.s05.osr.05)

ychrom.comm100 <- rbind(comm100.osr1, comm100.osr.8, comm100.osr.6, comm100.osr.4, comm100.osr.2, 
                        comm100.osr.1, comm100.osr.05)

# Then 50

freq.additive.s05.osr1 <- res.base$pop50$rd0.1$h0.5$s0.5[,2]
comm50.osr1 <- mean(freq.additive.s05.osr1)

freq.additive.s05.osr.8 <- res.rare.fem$males50$females0.8$rd0.1$h0.5$s0.5[,2]
comm50.osr.8 <- mean(freq.additive.s05.osr.8)

freq.additive.s05.osr.6 <- res.rare.fem$males50$females0.6$rd0.1$h0.5$s0.5[,2]
comm50.osr.6 <- mean(freq.additive.s05.osr.6)

freq.additive.s05.osr.4 <- res.rare.fem$males50$females0.4$rd0.1$h0.5$s0.5[,2]
comm50.osr.4 <- mean(freq.additive.s05.osr.4)

freq.additive.s05.osr.2 <- res.rare.fem$males50$females0.2$rd0.1$h0.5$s0.5[,2]
comm50.osr.2 <- mean(freq.additive.s05.osr.2)

freq.additive.s05.osr.1 <- res.rare.fem$males50$females0.1$rd0.1$h0.5$s0.5[,2]
comm50.osr.1 <- mean(freq.additive.s05.osr.1)

freq.additive.s05.osr.05 <- res.rare.fem$males50$females0.05$rd0.1$h0.5$s0.5[,2]
comm50.osr.05 <- mean(freq.additive.s05.osr.05)

ychrom.comm50 <- rbind(comm50.osr1, comm50.osr.8, comm50.osr.6, comm50.osr.4, comm50.osr.2, 
                       comm50.osr.1, comm50.osr.05)

ychrom.freq <- rbind(ychrom.comm1000, ychrom.comm500, ychrom.comm100, ychrom.comm50)

ychrom.comm.osr[, 1] <- ychrom.freq

ychrom.freq <- as.data.frame(ychrom.comm.osr, stringsAsFactors = F)

ggplot(data = ychrom.freq, aes(x = log.Ne, y = as.numeric(ychrom.freq$freq))) +
  geom_point(aes(col=as.factor(osr)) , stat="identity", position="jitter", alpha=0.8, size=4) +
  scale_x_continuous(breaks =c(0,1,2,3), labels=c(c(1,10,100,1000))) +
  xlab("Variance effective population size") +
  ylab("Mean allele frequency") +
  guides(fill=guide_legend(title="OSR")) +
  scale_colour_viridis_d("OSR") +
  ggtitle("Mean frequency of the allele benefitting the common sex (males) on the Y chromosome.  
Genetic architecture is additive. s= 0.5")


write.csv(ychrom.freq, file= "freq.ychrom.additive.s05.csv")

