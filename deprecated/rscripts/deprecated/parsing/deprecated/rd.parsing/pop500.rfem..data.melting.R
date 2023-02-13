# Julio Rincones-Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com

# Merging then melting the parsed results for plotting. Pop size 500. Males are the common sex

pop500.h0 <- read.csv("../sex.bias/rscripts/parsing/rd.parsing/rare.females/pop500.h0.csv")
pop500.h05 <- read.csv("../sex.bias/rscripts/parsing/rd.parsing/rare.females/pop500.h05.csv")
pop500.h1 <- read.csv("../sex.bias/rscripts/parsing/rd.parsing/rare.females/pop500.h1.csv")

pop500.merged = merge(pop500.h0, pop500.h05, by="rd")
pop500.merged.h = merge(pop500.merged, pop500.h1, by = "rd")
head(pop500.merged.h)


pop500.melted <- reshape2::melt(pop.500, id.var='rd')
head(pop500.melted)
write.csv(pop500.melted, file = "pop500.melted.csv")

# library(ggraptR)
# ggraptR(pop500.melted)
# ggplot(pop500.melted, aes(x=rd, y=value, col=variable)) + geom_line() + xlab('rd') +
#   ylab('fitness ratio')
# 
# p = ggplot() + 
#   geom_line(data = pop500.h0, aes(x = rd, y = OSR.2), color = "blue") +
#   geom_line(data = pop500.h0, aes(x  = rd, y = OSR.1), color = "blue") +
#   geom_line(data = pop500.h0, aes(x =  rd, y = OSR.05), color = "blue") +
#   geom_line(data = pop500.h05, aes(x = rd, y = OSR.2), color = "black") +
#   geom_line(data = pop500.h05, aes(x  = rd, y = OSR.1), color = "black") +
#   geom_line(data = pop500.h05, aes(x =  rd, y = OSR.05), color = "black") +
#   geom_line(data = pop500.h1, aes(x = rd, y = OSR.2), color = "orange") +
#   geom_line(data = pop500.h1, aes(x  = rd, y = OSR.1), color = "orange") +
#   geom_line(data = pop500.h1, aes(x =  rd, y = OSR.05), color = "orange") +
#   xlab('rd') +
#   ylab('fitness ratio') +
# print(p)
