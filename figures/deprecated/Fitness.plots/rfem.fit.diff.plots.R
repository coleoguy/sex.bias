# Julio Rincones Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com
# P.I. Dr. Heath Blackmon
# coleoguy@gmail.com

# Script to plot fitness differences between males and females for the XY model. The aim is to facet 
# genetic architecture 

library(ggplot2)
library(ggraptR)
library(viridis)

#set your working directory to sex.bias/rscripts
# load  the individual files for pop size 100
pop.100.h0 <- read.csv("../../sex.bias/rscripts/parsing/rd.parsing/rare.females/pop100.h0.csv")
pop.100.h05 <- read.csv("../../sex.bias/rscripts/parsing/rd.parsing/rare.females/pop100.h05.csv")
pop.100.h1 <- read.csv("../../sex.bias/rscripts/parsing/rd.parsing/rare.females/pop100.h1.csv")

# merge them
pop100.merged = merge(pop.100.h0, pop.100.h05, by="rd")
pop100.merged = merge(pop100.merged, pop.100.h1, by = "rd")

# melt them

pop100.melted.test <- reshape2::melt(pop100.merged, id.var='rd')

# need to do add the osr column and a genetic architecture column then the variable column can be replaced 
# to look like the one loaded on line 45

# save them

write.csv(pop100.melted.test, file = "../Fitness.plots/XY.model/cvs.rare.fem/pop100.rfem.melted.test.csv")

# the file needs to be amended so there is a column for rd values, one column for genetic 
# architecture, one column for fitness difference, and one column for OSR values corresponding
# to each genetic architecture. I have done this manually but the values can be saved as a vector
# and inserted into the data frame.

# I have reloaded the .csv file to plot in ggraptR  and then adjust the  colour paletter to be
# colourblind friendly.

test <- pop100.melted.test

# if loading the data
test <- read.csv("figures/Fitness.plots/XY.model/rd.plots/pop100.rfem.melted.test.csv", 
                 row.names = 1, as.is = T, header = T, check.names = F)

pop100 <- ggplot(test, aes(y= Fitness.difference, x= rd)) + 
  geom_point(aes(colour= architecture), stat="identity", position="identity", alpha=1, size=2.5) + 
  scale_color_viridis_d(begin = 0, end = 1)+
  facet_grid(OSR ~ .) + theme_light() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=15, 
                          hjust=0.5, vjust=0.5)) + scale_size(range=c(1, 3)) + 
  ggtitle("Fitness difference when males are the common sex. 
Males = 100. OSR = 0.2, 0.1, and 0.05") + 
  xlab("recombination distance (rd)") + ylab("Fitness difference") 

ggsave(filename = "xy.fitdiff.rd.rfem.100.pdf",
       plot= pop100, width=8, height=7, units="in", path = "../Figures/rd.plots")
 

# Now load  the individual files for pop size 500
pop.500.h0 <- read.csv("../sex.bias/rscripts/parsing/rd.parsing/rare.females/pop500.h0.csv")
pop.500.h05 <- read.csv("../sex.bias/rscripts/parsing/rd.parsing/rare.females/pop500.h05.csv")
pop.500.h1 <- read.csv("../sex.bias/rscripts/parsing/rd.parsing/rare.females/pop500.h1.csv")

# merge them
pop500.merged = merge(pop.500.h0, pop.500.h05, by="rd")
pop500.merged = merge(pop500.merged, pop.500.h1, by = "rd")

# melt them

pop500.melted.test <- reshape2::melt(pop500.merged, id.var='rd')

# save them

write.csv(pop500.melted.test, file = "pop500.rfem.melted.test.csv")

# the file needs to be amended so there is a column for rd values, one column for genetic 
# architecture, one column for fitness difference, and one column for OSR values corresponding
# to each genetic architecture. I have done this manually but the values can be saved as a vector
# and inserted into the data frame.

# I have reloaded the .csv file to plot in ggraptR  and then adjust the  colour palette to be
# colourblind friendly.

test <- read.csv("../sex.bias/pop500.rfem.melted.test.csv", 
                 row.names = 1, as.is = T, header = T, check.names = F)

ggraptR(test)

ggplot(test, aes(y= Fitness.difference, x= rd)) + 
  geom_point(aes(colour= architecture), stat="identity", position="identity", alpha=1, size=3) + 
  scale_color_viridis_d(begin = 0, end = 1)+
  facet_grid(OSR ~ .) + theme_grey() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=15, 
                          hjust=0.5, vjust=0.5)) + scale_size(range=c(1, 3)) + 
  ggtitle("Fitness difference when males are the common sex. 
          Males = 500. OSR = 0.2, 0.1, and 0.05") + 
  xlab("recombination distance (rd)") + ylab("Fitness difference") 

