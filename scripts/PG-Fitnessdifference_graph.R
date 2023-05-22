
setwd("C:/Users/pdglenn/OneDrive/Github/sex.bias/results")
#Load libraries
library(gridExtra)
library(grid)
library(ggpubr)
library(ggplot2)
library(viridis)
library(vioplot)
library(cowplot)

# first we load our functions
source("../scripts/functions.R")


#### Rare Male analysis, Common Female, allele 2 is beneficial for common sex (female) ####

mdat<-read.csv("rare.male.250.iter.csv", as.is=T, header = TRUE)
simple.mdat <- aggregate(x = mdat,
                         by = list(mdat$females, mdat$OSR, mdat$rd, mdat$h, mdat$s),
                         FUN = mean)[,-c(1:5)]


simple.mdat$Wdiff <- simple.mdat$malW <- simple.mdat$femW <- NA

for(i in 1:nrow(simple.mdat)){
  fits <- getFits(simple.mdat[i,])
  simple.mdat[i,11] <- fits[1]
  simple.mdat[i,10] <- fits[2]
  simple.mdat[i,12] <- fits[3]
  
}

simple.mdat$delA <- simple.mdat$fixA <-simple.mdat$delY <-simple.mdat$fixY <- simple.mdat$delX <-simple.mdat$fixX <-  NA

for(i in 1:nrow(simple.mdat)){
  fixes <- getProb(simple.mdat, mdat)
  simple.mdat[i,13] <- fixes[1]
  simple.mdat[i,14] <- fixes[2]
  simple.mdat[i,15] <- fixes[3]
  simple.mdat[i,16] <- fixes[4]
  simple.mdat[i,17] <- fixes[5]
  simple.mdat[i,18] <- fixes[6]
}


#Sex chromosome

simple.msex<-simple.mdat[simple.mdat$rd == 0.2,] #Sex-linked

#Fitness

#Exclude h = 99
mfit.sex <- simple.msex[simple.msex$h != 99,]

mfit.sex$females<-as.factor(mfit.sex$females)
mfit.sex$s<-as.factor(mfit.sex$s)
mfit.sex$h<-as.factor(mfit.sex$h)
str(mfit.sex)

#h = 0.5

#mfit.sex.5h <- mfit.sex[mfit.sex$h == 0.5,]


#### Female 

#rare female, common male

fdat<-read.csv("rare.female.250.iter.csv", as.is=T, header = TRUE)
simple.fdat <- aggregate(x = fdat,
                         by = list(fdat$males, fdat$OSR, fdat$rd, fdat$h, fdat$s),
                         FUN = mean)[,-c(1:5)]


simple.fdat$Wdiff <- simple.fdat$malW <- simple.fdat$femW <- NA

for(i in 1:nrow(simple.fdat)){
  fits <- getFits(simple.fdat[i,])
  simple.fdat[i,11] <- fits[1]
  simple.fdat[i,10] <- fits[2]
  simple.fdat[i,12] <- fits[3]
  
}


simple.fdat$delA <- simple.fdat$fixA <-simple.fdat$delY <-simple.fdat$fixY <- simple.fdat$delX <-simple.fdat$fixX <-  NA

for(i in 1:nrow(simple.fdat)){
  fixes <- getProb(simple.fdat, fdat)
  simple.fdat[i,13] <- fixes[1]
  simple.fdat[i,14] <- fixes[2]
  simple.fdat[i,15] <- fixes[3]
  simple.fdat[i,16] <- fixes[4]
  simple.fdat[i,17] <- fixes[5]
  simple.fdat[i,18] <- fixes[6]
}



simple.fsex<-simple.fdat[simple.fdat$rd == 0.2,] #Sex-linked

#remove 99
ffit.sex <- simple.fsex[simple.fsex$h != 99,]

ffit.sex$males<-as.factor(ffit.sex$males)
ffit.sex$s<-as.factor(ffit.sex$s)
ffit.sex$h<-as.factor(ffit.sex$h)
#str(ffit.sex)

#ffit.sex.5h <- ffit.sex[ffit.sex$h == 0.5,]

#Combine the two together since everything is in the same order
SexFit <- mfit.sex[,c(4:8,12)]
colnames(SexFit)[1] <- "common"
colnames(SexFit)[6] <- "M_Wdiff"
SexFit[7] <- ffit.sex[12]
colnames(SexFit)[7] <- "F_Wdiff"


#### Autosome ####

### Male


simple.maut<-simple.mdat[simple.mdat$rd == 0.5,] #Autosome


#Fitness

#Exclude h = 99
mfit.aut <- simple.maut[simple.maut$h != 99,]

mfit.aut$females<-as.factor(mfit.aut$females)
mfit.aut$s<-as.factor(mfit.aut$s)
mfit.aut$h<-as.factor(mfit.aut$h)
#str(mfit.aut)

#h = 0.5

#mfit.aut.5h <- mfit.aut[mfit.aut$h == 0.5,]

### Female
simple.faut<-simple.fdat[simple.fdat$rd == 0.5,] #Autosome

#eclude h =99
ffit.aut <- simple.faut[simple.faut$h != 99,]

ffit.aut$males<-as.factor(ffit.aut$males)
ffit.aut$s<-as.factor(ffit.aut$s)
ffit.aut$h<-as.factor(ffit.aut$h)
#str(ffit.aut)

#ffit.aut.5h <- ffit.aut[ffit.aut$h == 0.5,]


#Combine the two together since everything is in the same order
AutFit <- mfit.aut[,c(4:8,12)]
colnames(AutFit)[1] <- "common"
colnames(AutFit)[6] <- "M_Wdiff"
AutFit[7] <- ffit.aut[12]
colnames(AutFit)[7] <- "F_Wdiff"


#Graphs


my_labels = c("1.0","0.75","0.50","0.25","0.00","0.25","0.50","0.75","1.0")
my_ticks = c(-1.0,-0.75,-0.50,-0.25,0,0.25,0.5, 0.75,1.0)

SexFit_h <- SexFit[SexFit$h == 0.5,]

SexFitnessComp <- ggplot(SexFit_h, aes(x = OSR))+
  facet_grid(h~common)+
  geom_line(aes(y=M_Wdiff, colour =s))+
  geom_point(aes(y=M_Wdiff, shape = s, colour = s))+
  geom_point(aes(y=F_Wdiff, shape = s, colour = s))+
  geom_line(aes(y=F_Wdiff, colour = s))+
  scale_y_continuous(limits = c(-1,1), breaks = my_ticks, labels=my_labels,expand = c(0,0))+
  theme(axis.line.y = element_blank(), plot.margin = unit(c(1,1,1,1), "lines"))+
  ylab(expression("Fitness difference (common - rare)\n"))+
  ggtitle("Fitness, rd = 0.2, h = 0.5")
SexFitnessComp <- ggdraw() + draw_plot(SexFitnessComp) +
  draw_label("Common Male", x = 0.03, y = 0.65, angle = 90, size =10) +
  draw_label("Common Female", x = 0.03, y = 0.30, angle = 90, size =10)
SexFitnessComp


AutFit_h <- AutFit[AutFit$h == 0.5,]

AutFitnessComp <- ggplot(AutFit_h, aes(x = OSR))+
  facet_grid(h~common)+
  geom_line(aes(y=M_Wdiff, colour =s,))+
  geom_point(aes(y=M_Wdiff, shape = s, colour = s))+
  geom_point(aes(y=F_Wdiff, shape = s, colour = s))+
  geom_line(aes(y=F_Wdiff, colour = s))+
  scale_y_continuous(limits = c(-1,1.0), breaks = my_ticks, labels=my_labels,expand = c(0,0))+
  theme(axis.line.y = element_blank(), plot.margin = unit(c(1,1,1,1), "lines"))+
  ylab(expression("Fitness difference (common - rare)\n"))+
  ggtitle("Fitness, rd = 0.5, h = 0.5") 
AutFitnessComp <- ggdraw() + draw_plot(AutFitnessComp) +
  draw_label("Common Male", x = 0.03, y = 0.65, angle = 90, size =10) +
  draw_label("Common Female", x = 0.03, y = 0.30, angle = 90, size =10)
AutFitnessComp


### Fixation 

# X
# set up for allele 2 which is beneficial for females
# delA is allele 1 (male)
# fixA is allele 2 (female)

AutFix <- mfit.aut[,c(4:8,13:18)]
colnames(AutFix)[1] <- "common"
colnames(AutFix)[c(6:11)] <- c("M_FfixX","M_FdelX", "M_FfixY","M_FdelY","M_FfixA","M_FdelA")
AutFix[c(12:17)] <- ffit.aut[c(13:18)]
colnames(AutFix)[c(12:17)] <- c("F_MfixX","F_MdelX", "F_MfixY","F_MdelY","F_MfixA","F_MdelA")

AutFix_h <- AutFix[AutFix$h == 0.5,]

AutFixCompA <- ggplot(AutFix, aes(x = OSR))+
  facet_grid(h~common)+
  geom_line(aes(y=(M_FdelA), colour =s,))+ #common female (allele 2)
  geom_point(aes(y=(M_FdelA), shape = "Common Female", colour = s))+
  geom_point(aes(y=(F_MfixA), shape = "Common Male", colour = s))+ #common male, allele 1 
  geom_line(aes(y=(F_MfixA), colour = s))+
  scale_y_continuous(limits = c(-.1,1.0), breaks = my_ticks, labels=my_labels,expand = c(0,0))+
  theme(axis.line.y = element_blank(), plot.margin = unit(c(1,1,1,1), "lines"))+
  ylab(expression("Frequeny of beneficial A \n"))+
  ggtitle("Prob of beneficial A Fixed, rd = 0.5, h = 0.5") 
AutFixCompA


SexFix <- mfit.sex[,c(4:8,13:18)]
colnames(SexFix)[1] <- "common"
colnames(SexFix)[c(6:11)] <- c("M_FfixX","M_FdelX", "M_FfixY","M_FdelY","M_FfixA","M_FdelA")
SexFix[c(12:17)] <- ffit.aut[c(13:18)]
colnames(SexFix)[c(12:17)] <- c("F_MfixX","F_MdelX", "F_MfixY","F_MdelY","F_MfixA","F_MdelA")

#X chromosome
SexFixCompX <- ggplot(SexFix, aes(x = OSR))+
  facet_grid(h~common)+
  geom_line(aes(y=(M_FfixX), colour =s,))+ #common female (allele 2) X - allele 2 is ben (fixed)
  geom_point(aes(y=(M_FfixX), shape = "Common Female", colour = s))+
  geom_point(aes(y=(F_MdelX), shape = "Common Male", colour = s))+ #common male, allele 1 - delX
  geom_line(aes(y=(F_MdelX), colour = s))+
  scale_y_continuous(limits = c(-.1,1.0), breaks = my_ticks, labels=my_labels,expand = c(0,0))+
  theme(axis.line.y = element_blank(), plot.margin = unit(c(1,1,1,1), "lines"))+
  ylab(expression("Frequeny of beneficial X \n"))+
  ggtitle("Prob of beneficial X Fixed, rd = 0.5, h = 0.5") 
SexFixCompX

#The two sex chromosome probabilities of fixation are the same in both male and female coomon/rare
#Why?
#Male has one X and Y
#Probability of fixing Y is solely dependent on the male so it seems reasonable that the probability of
#Fixing the beneficial Y depends on drift and whether or not ben Y is recessive or dominant
#Check the original graphs

#Y chromosome
SexFixCompY <- ggplot(SexFix, aes(x = OSR))+
  facet_grid(h~common)+
  geom_line(aes(y=(M_FdelY), colour =s,))+ #common female (allele 2) - allele 2 is ben (delY)
  geom_point(aes(y=(M_FdelY), shape = "Common Female", colour = s))+
  geom_point(aes(y=(F_MfixY), shape = "Common Male", colour = s))+ #common male, allele 1 - fixY
  geom_line(aes(y=(F_MfixY), colour = s))+
  scale_y_continuous(limits = c(-.1,1.0), breaks = my_ticks, labels=my_labels,expand = c(0,0))+
  theme(axis.line.y = element_blank(), plot.margin = unit(c(1,1,1,1), "lines"))+
  ylab(expression("Frequeny of beneficial Y \n"))+
  ggtitle("Prob of beneficial Y Fixed, rd = 0.5, h = 0.5") 
SexFixCompY



SexFixComp <- ggplot(SexFix, aes(x = OSR))+
  facet_grid(h~common)+
  geom_line(aes(y=(M_FfixX), colour =s,))+ #common female (allele 2) X - fixX 
  geom_point(aes(y=(M_FfixX), shape = "Common Female", colour = s))+
  geom_point(aes(y=(F_MfixY), shape = "Common Male", colour = s))+ #common male, allele 1 - fixY
  geom_line(aes(y=(F_MfixY), colour = s))+
  scale_y_continuous(limits = c(-.1,1.0), breaks = my_ticks, labels=my_labels,expand = c(0,0))+
  theme(axis.line.y = element_blank(), plot.margin = unit(c(1,1,1,1), "lines"))+
  ylab(expression("Frequeny of beneficial Chromosome \n"))+
  ggtitle("Prob of beneficial Chr Fixed, rd = 0.5, h = 0.5") 
SexFixComp

SexFixCompDel <- ggplot(SexFix, aes(x = OSR))+
  facet_grid(h~common)+
  geom_line(aes(y=(M_FdelX), colour =s,))+ #common female (allele 2) X - fixX 
  geom_point(aes(y=(M_FdelX), shape = "Common Female", colour = s))+
  geom_point(aes(y=(F_MdelY), shape = "Common Male", colour = s))+ #common male, allele 1 - fixY
  geom_line(aes(y=(F_MdelY), colour = s))+
  scale_y_continuous(limits = c(-.1,1.0), breaks = my_ticks, labels=my_labels,expand = c(0,0))+
  theme(axis.line.y = element_blank(), plot.margin = unit(c(1,1,1,1), "lines"))+
  ylab(expression("Frequeny of deleterious Chromosome \n"))+
  ggtitle("Prob of beneficial Chr deleted, rd = 0.5, h = 0.5") 
SexFixCompDel
