library(ggplot2)
library(viridis)
library(ggraptR)

# Freq when there is no OSR bias
load("../results/base.comp.model.RData")
auto.freq <- results

# Freq when there is bias
load("../results/rare.male.model.RData")
auto.freq.osr <- results


# Freq when there is no bias (uses the base.comp.model.RData file)
freq.additive.s05.osr1 <- auto.freq$pop1000$rd0.5$h0.5$s0.5
auto.1000.osr1 <- freq.additive.s05.osr1[, -c(1,2)]
comm1000.osr1 <- mean(auto.1000.osr1)

# First at 1000 comm sex
freq.additive.osr.8 <- auto.freq.osr$females1000$males0.8$rd0.5$h0.5$s0.5
auto.1000.osr.8 <- freq.additive.osr.8[, -c(1,2)]
comm1000.osr.8 <- mean(auto.1000.osr.8)

freq.additive.osr.6 <- auto.freq.osr$females1000$males0.6$rd0.5$h0.5$s0.5
auto.1000.osr.6 <- freq.additive.osr.6[, -c(1,2)]
comm1000.osr.6 <- mean(auto.1000.osr.6)

freq.additive.osr.4 <- auto.freq.osr$females1000$males0.4$rd0.5$h0.5$s0.5
auto.1000.osr.4 <- freq.additive.osr.4[, -c(1,2)]
comm1000.osr.4 <- mean(auto.1000.osr.4)

freq.additive.osr.2 <- auto.freq.osr$females1000$males0.2$rd0.5$h0.5$s0.5
auto.1000.osr.2 <- freq.additive.osr.2[, -c(1,2)]
comm1000.osr.2 <- mean(auto.1000.osr.2)

freq.additive.osr.1 <- auto.freq.osr$females1000$males0.1$rd0.5$h0.5$s0.5
auto.1000.osr.1 <- freq.additive.osr.1[, -c(1,2)]
comm1000.osr.1 <- mean(auto.1000.osr.1)

freq.additive.osr.05 <- auto.freq.osr$females1000$males0.05$rd0.5$h0.5$s0.5
auto.1000.osr.05 <- freq.additive.osr.05[, -c(1,2)]
comm1000.osr.05 <- mean(auto.1000.osr.05)

autosomes.comm1000 <- rbind(comm1000.osr1, comm1000.osr.8, comm1000.osr.6, comm1000.osr.4, comm1000.osr.2, comm1000.osr.1, comm1000.osr.05)

# Then 500

freq.additive.s05.osr1 <- auto.freq$pop500$rd0.5$h0.5$s0.5
auto.500.osr1 <- freq.additive.s05.osr1[, -c(1,2)]
comm500.osr1 <- mean(auto.500.osr1)

freq.additive.osr.8 <- auto.freq.osr$females500$males0.8$rd0.5$h0.5$s0.5
auto.500.osr.8 <- freq.additive.osr.8[, -c(1,2)]
comm500.osr.8 <- mean(auto.500.osr.8)

freq.additive.osr.6 <- auto.freq.osr$females500$males0.6$rd0.5$h0.5$s0.5
auto.500.osr.6 <- freq.additive.osr.6[, -c(1,2)]
comm500.osr.6 <- mean(auto.500.osr.6)

freq.additive.osr.4 <- auto.freq.osr$females500$males0.4$rd0.5$h0.5$s0.5
auto.500.osr.4 <- freq.additive.osr.4[, -c(1,2)]
comm500.osr.4 <- mean(auto.500.osr.4)

freq.additive.osr.2 <- auto.freq.osr$females500$males0.2$rd0.5$h0.5$s0.5
auto.500.osr.2 <- freq.additive.osr.2[, -c(1,2)]
comm500.osr.2 <- mean(auto.500.osr.2)

freq.additive.osr.1 <- auto.freq.osr$females500$males0.1$rd0.5$h0.5$s0.5
auto.500.osr.1 <- freq.additive.osr.1[, -c(1,2)]
comm500.osr.1 <- mean(auto.500.osr.1)

freq.additive.osr.05 <- auto.freq.osr$females500$males0.05$rd0.5$h0.5$s0.5
auto.500.osr.05 <- freq.additive.osr.05[, -c(1,2)]
comm500.osr.05 <- mean(auto.500.osr.05)

autosomes.comm500 <- rbind(comm500.osr1, comm500.osr.8, comm500.osr.6, comm500.osr.4, comm500.osr.2, comm500.osr.1, comm500.osr.05)

# Then 100

freq.additive.s05.osr1 <- auto.freq$pop100$rd0.5$h0.5$s0.5
auto.100.osr1 <- freq.additive.s05.osr1[, -c(1,2)]
comm100.osr1 <- mean(auto.100.osr1)

freq.additive.osr.8 <- auto.freq.osr$females100$males0.8$rd0.5$h0.5$s0.5
auto.100.osr.8 <- freq.additive.osr.8[, -c(1,2)]
comm100.osr.8 <- mean(auto.100.osr.8)

freq.additive.osr.6 <- auto.freq.osr$females100$males0.6$rd0.5$h0.5$s0.5
auto.100.osr.6 <- freq.additive.osr.6[, -c(1,2)]
comm100.osr.6 <- mean(auto.100.osr.6)

freq.additive.osr.4 <- auto.freq.osr$females100$males0.4$rd0.5$h0.5$s0.5
auto.100.osr.4 <- freq.additive.osr.4[, -c(1,2)]
comm100.osr.4 <- mean(auto.100.osr.4)

freq.additive.osr.2 <- auto.freq.osr$females100$males0.2$rd0.5$h0.5$s0.5
auto.100.osr.2 <- freq.additive.osr.2[, -c(1,2)]
comm100.osr.2 <- mean(auto.100.osr.2)

freq.additive.osr.1 <- auto.freq.osr$females100$males0.1$rd0.5$h0.5$s0.5
auto.100.osr.1 <- freq.additive.osr.1[, -c(1,2)]
comm100.osr.1 <- mean(auto.100.osr.1)

freq.additive.osr.05 <- auto.freq.osr$females100$males0.05$rd0.5$h0.5$s0.5
auto.100.osr.05 <- freq.additive.osr.05[, -c(1,2)]
comm100.osr.05 <- mean(auto.100.osr.05)

autosomes.comm100 <- rbind(comm100.osr1, comm100.osr.8, comm100.osr.6, comm100.osr.4, comm100.osr.2, comm100.osr.1, comm100.osr.05)

# Then 50
freq.additive.s05.osr1 <- auto.freq$pop50$rd0.5$h0.5$s0.5
auto.50.osr1 <- freq.additive.s05.osr1[, -c(1,2)]
comm50.osr1 <- mean(auto.50.osr1)

freq.additive.osr.8 <- auto.freq.osr$females50$males0.8$rd0.5$h0.5$s0.5
auto.50.osr.8 <- freq.additive.osr.8[, -c(1,2)]
comm50.osr.8 <- mean(auto.50.osr.8)

freq.additive.osr.6 <- auto.freq.osr$females50$males0.6$rd0.5$h0.5$s0.5
auto.50.osr.6 <- freq.additive.osr.6[, -c(1,2)]
comm50.osr.6 <- mean(auto.50.osr.6)

freq.additive.osr.4 <- auto.freq.osr$females50$males0.4$rd0.5$h0.5$s0.5
auto.50.osr.4 <- freq.additive.osr.4[, -c(1,2)]
comm50.osr.4 <- mean(auto.50.osr.4)

freq.additive.osr.2 <- auto.freq.osr$females50$males0.2$rd0.5$h0.5$s0.5
auto.50.osr.2 <- freq.additive.osr.2[, -c(1,2)]
comm50.osr.2 <- mean(auto.50.osr.2)

freq.additive.osr.1 <- auto.freq.osr$females50$males0.1$rd0.5$h0.5$s0.5
auto.50.osr.1 <- freq.additive.osr.1[, -c(1,2)]
comm50.osr.1 <- mean(auto.50.osr.1)

freq.additive.osr.05 <- auto.freq.osr$females50$males0.05$rd0.5$h0.5$s0.5
auto.50.osr.05 <- freq.additive.osr.05[, -c(1,2)]
comm50.osr.05 <- mean(auto.50.osr.05)

autosomes.comm50 <- rbind(comm50.osr1, comm50.osr.8, comm50.osr.6, comm50.osr.4, comm50.osr.2, comm50.osr.1, comm50.osr.05)

autosome.freq <- rbind(autosomes.comm1000, autosomes.comm500, autosomes.comm100, autosomes.comm50)

auto.comm.osr <- matrix(,28,4)

colnames(auto.comm.osr) <- c("freq", "commsex", "osr", "Ne")
auto.comm.osr[, 1] <- autosome.freq[,1]

commsex <- rep(c(1000, 500, 100, 50), each = 7)
auto.comm.osr[, 2] <- commsex

osr <- rep(c("osr1", "osr.8", "osr.6", "osr.4", "osr.2", "osr.1", "osr.05"), each = 4)

Ne <- var_ne_autosomes

log.Ne <- log10(Ne.autosomes)

Ne.autosomes <- c()
for(i in 1:nrow(Ne)){
  Ne.autosomes <- c(Ne.autosomes, unlist(Ne[i,4:10]))
}
auto.comm.osr[,4] <- Ne.autosomes

auto.comm.osr <- cbind(auto.comm.osr, log.Ne)

autosome.freq <- as.data.frame(auto.comm.osr, stringsAsFactors = F)
class(autosome.freq$Ne) <- "numeric"

autosome.freq$osr[autosome.freq$osr == "osr.05"] <- .05
autosome.freq$osr[autosome.freq$osr == "osr.1"] <- .1
autosome.freq$osr[autosome.freq$osr == "osr.2"] <- .2
autosome.freq$osr[autosome.freq$osr == "osr.4"] <- .4
autosome.freq$osr[autosome.freq$osr == "osr.6"] <- .6
autosome.freq$osr[autosome.freq$osr == "osr.8"] <- .8
autosome.freq$osr[autosome.freq$osr == "osr1"] <- 1

ggplot(data = autosome.freq, aes(x = log.Ne, y = 1-as.numeric(autosome.freq$freq)))+
  geom_point(aes(col=as.factor(osr)) , stat="identity", position="jitter", alpha=0.8, size=4) +
  scale_x_discrete(breaks =c(0,1,2,3), labels=c(c(1,10,100,1000))) +
  xlab("Variance effective population size") +
  ylab("Mean allele frequency") +
  guides(fill=guide_legend(title="OSR")) +
  scale_colour_viridis_d("OSR") +
  ggtitle("Mean frequency of the allele benefitting the common sex (females) on Autosomes.  
Genetic architecture is additive. s= 0.5")

write.csv(autosome.freq, file= "freq.auto.additive.s05.csv")


