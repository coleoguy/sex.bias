library(ggplot2)
library(viridis)
library(ggraptR)

xchrom.freq <- read.csv("freq.xchrom.additive.s05.csv")

ggplot(data = xchrom.freq, aes(x = log.Ne, y = as.numeric(xchrom.freq$freq))) +
  geom_point(aes(col=as.factor(osr)) , stat="identity", position="jitter", alpha=0.8, size=4) +
  scale_x_continuous(breaks =c(0,1,2,3), labels=c(c(1,10,100,1000))) +
  xlab("Variance effective population size") +
  ylab("Allele frequency of the female beneficial allele") +
  guides(fill=guide_legend(title="OSR")) +
  scale_colour_viridis_d("OSR") +
  ggtitle("Allele frequency on X chromosomes. h = 0.5. s= 0.5")


ychrom.freq <- read.csv("freq.ychrom.additive.s05.csv")

ggplot(data = ychrom.freq, aes(x = log.Ne, y = as.numeric(ychrom.freq$freq))) +
  geom_point(aes(col=as.factor(osr)) , stat="identity", alpha=0.8, size=4) +
  scale_x_continuous(breaks =c(0,1,2,3), labels=c(c(1,10,100,1000))) +
  xlab("Variance effective population size") +
  ylab("Mean allele frequency") +
  guides(fill=guide_legend(title="OSR")) +
  scale_colour_viridis_d("OSR") +
  ggtitle("Mean frequency of the allele benefitting the common\nsex (males) on the Y chromosome. h = 0.5. s= 0.5")

  autosome.freq <- read.csv("freq.auto.additive.s05.csv")
  
ggplot(data = autosome.freq, aes(x = log.Ne, y = 1-as.numeric(autosome.freq$freq)))+
  geom_point(aes(col=as.factor(osr)) , stat="identity", alpha=0.8, size=4) +
  scale_x_continuous(breaks =c(0,1,2,3), labels=c(c(1,10,100,1000))) +
  xlab("Variance effective population size") +
  ylab("Allele frequency of the female beneficial allele") +
  guides(fill=guide_legend(title="OSR")) +
  scale_colour_viridis_d("OSR") +
  ggtitle("Allele frequency on Autosomes. h = 0.5. s= 0.5")
