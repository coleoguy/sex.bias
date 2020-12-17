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

# load  the individual files for pop size 100
no.bias.h0 <- read.csv("../sex.bias/rscripts/parsing/rd.parsing/nobias/h0.nobias.XYmodel.csv", 
                       row.names = 1, as.is = T, header = T, check.names = F)
no.bias.h05 <- read.csv("../sex.bias/rscripts/parsing/rd.parsing/nobias/h0.5.nobias.XYmodel.csv", 
                        row.names = 1, as.is = T, header = T, check.names = F)
no.bias.h1 <- read.csv("../sex.bias/rscripts/parsing/rd.parsing/nobias/h1.nobias.XYmodel.csv", 
                       row.names = 1, as.is = T, header = T, check.names = F)



# merge them
merged.dat1 = merge(no.bias.h0, no.bias.h05, by = "rd")
merged.dat = merge(merged.dat1, no.bias.h1, by = "rd")
rm(merged.dat1)

merged.col.names <- c("rd", "100.h0", "500.h0", "100.h0.5", "500.h0.5", "100.h1", "500.h1")
colnames(merged.dat) <- merged.col.names

# melt them
melted.dat <- reshape2::melt(merged.dat, id.var='rd')

# save them
write.csv(melted.dat, file = "no.bias.melted.data.csv")


# the file needs to be amended so there is a column for rd values, one column for genetic 
# architecture, one column for fitness difference, and one column for OSR values corresponding
# to each genetic architecture. I have done this manually but the values can be saved as a vector
# and inserted into the data frame.

# I have reloaded the .csv file to plot in ggraptR  and then adjust the  colour palette to be
# colourblind friendly.

ggraptR(melted.dat)

test <- read.csv("../sex.bias/figures/Fitness.plots/XY.model/rd.plots/no.bias.melted.data.csv", 
                 row.names = 1, as.is = T, header = T, check.names = F)

ggplot(test, aes(y= value, x= rd)) + 
  geom_point(aes(colour= variable), stat="identity", position="identity", alpha= 0.5, size=3) + 
  scale_color_viridis_d(begin = 0, end = 1)+
  theme_grey() + theme(text=element_text(family="sans", face="plain", color="#000000", size=15, 
                                         hjust=0.5, vjust=0.5)) + scale_size(range=c(1, 3)) + 
  ggtitle("Fitness difference when OSR is 1. 
  Population size = 500 and 100. s = 0.5") + xlab("recombination distance (rd)") + 
  ylab("Fitness difference")


