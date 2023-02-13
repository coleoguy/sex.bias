# Julio Rincones-Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com

# set working directory  to sex.bias
# Merging then melting the parsed results for plotting. Pop size 100. Males are the common sex.

pop100.h0 <- read.csv("../../sex.bias/rscripts/parsing/rd.parsing/rare.females/pop100.h0.csv")
pop100.h05 <- read.csv("../../sex.bias/rscripts/parsing/rd.parsing/rare.females/pop100.h05.csv")
pop100.h1 <- read.csv("../../sex.bias/rscripts/parsing/rd.parsing/rare.females/pop100.h1.csv")

pop100.merged = merge(pop100.h0, pop100.h05, by="rd")
pop100.merged.h = merge(pop100.merged, pop100.h1, by = "rd")
head(pop100.merged.h)


pop100.melted <- reshape2::melt(pop100.merged.h, id.var='rd')
write.csv(pop100.melted, file = "pop100.melted.csv")


#Test
library(ggraptR)
ggraptR(pop100.melted)
ggplot(pop100.melted, aes(x=rd, y=value, col=variable)) + geom_line() + xlab('rd') +
  ylab('fitness ratio')

