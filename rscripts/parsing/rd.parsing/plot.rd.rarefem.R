# Julio Rincones-Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.comn

# experimenting with plotting the melted results of the different genetic architectures,
# fitnesses and OSR where we suspect rd values will show an interesting change in the direction of
# net masculinization towards net feminization.

library(ggplot2)
library(viridis)

pop500.melted <- read.csv("../../sex.bias/rscripts/parsing/rd.parsing/rare.females/pop500.melted.csv")


ggplot(pop500.melted, aes(y= value, x= rd)) + 
  geom_point(aes(colour= as.factor(variable)), stat="identity", alpha=0.4, size=3) + 
  theme_classic() + theme(text=element_text(family="sans", face="plain", color="#000000", size=15, 
                                         hjust=0.5, vjust=0.5)) + scale_size(range=c(1, 3)) + 
  guides(colour=guide_legend(title="OSR value 
and genetic 
architecture")) + 
  xlab("Recombination distance (rd)") + 
  ylab("Net masculinization (males/females)")


ggplot(pop500.melted, aes(x=rd, y=value, col=variable)) + geom_line() + xlab('rd') +
  ylab('fitness ratio')


pop100.melted <- read.csv("../../sex.bias/rscripts/parsing/rd.parsing/rare.females/pop100.melted.csv")


ggplot(pop100.melted, aes(y= value, x= rd)) + 
  geom_point(aes(colour= as.factor(variable)), stat="identity", alpha=0.4, size=3) + 
  theme_grey() + theme(text=element_text(family="sans", face="plain", color="#000000", size=15, 
                                         hjust=0.5, vjust=0.5)) + scale_size(range=c(1, 3)) + 
  guides(colour=guide_legend(title="OSR value 
and genetic 
architecture")) + 
  xlab("Recombination distance (rd)") + 
  ylab("Net masculinization (males/females)")
