
setwd("~/GitHub/sex.bias/results")
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


# test <-simple.mdat[simple.mdat$rd == 0.5,]
# test <-test[test$females == 1000,]
# test <-test[test$OSR == 1.0,,]
# test <-test[test$h == 0.5,]
# 
# 
# test2 <-simple.mdat[simple.mdat$rd == 0.2,]
# test2 <-test2[test2$females == 1000,]
# test2 <-test2[test2$OSR == 1.0,,]
# test2 <-test2[test2$h == 0.5,]



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


test4 <-SexFit_h[SexFit_h$rd == 0.2,]
test4 <-test4[test4$common == 1000,]
test4 <-test4[test4$OSR == 1.0,,]
test4 <-test4[test4$h == 0.5,]

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


test3 <-AutFit_h[AutFit_h$rd == 0.5,]
test3 <-test3[test3$common == 1000,]
test3 <-test3[test3$OSR == 1.0,,]
test3 <-test3[test3$h == 0.5,]



