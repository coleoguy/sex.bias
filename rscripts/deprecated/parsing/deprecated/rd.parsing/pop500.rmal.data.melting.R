# Julio Rincones-Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com

# Merging then melting the parsed results for plotting. Pop size 500. Femles are the common sex

pop500.h0 <- read.csv("../sex.bias/rscripts/parsing/rd.parsing/rare.males/pop500.rmal.h0.csv")
pop500.h05 <- read.csv("../sex.bias/rscripts/parsing/rd.parsing/rare.males/pop500.rmal.h05.csv")
pop500.h1 <- read.csv("../sex.bias/rscripts/parsing/rd.parsing/rare.males/pop500.rmal.h1.csv")

pop500.merged = merge(pop500.h0, pop500.h05, by="rd")
pop500.merged.h = merge(pop500.merged, pop500.h1, by = "rd")
head(pop500.merged.h)


pop500.melted <- reshape2::melt(pop500.merged.h, id.var='rd')
head(pop500.melted)
write.csv(pop500.melted, file = "pop500.melted.csv")
