# Priscilla Glenn
# OSR analysis script
# Updated 3-14-2023



#Load libraries
library(gridExtra)
library(grid)
library(ggpubr)
library(ggplot2)
library(viridis)
library(vioplot)

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




#separate by sex-linked or autosome



#### Rare Male Autosome ####

simple.maut<-simple.mdat[simple.mdat$rd == 0.5,] #Autosome

simple.maut$females<-as.factor(simple.maut$females)
simple.maut$h <- as.factor(simple.maut$h)
simple.maut$s <- as.factor(simple.maut$s)

str(simple.maut)

#Looking at rare male so Y-axis needs to be from the female/common sex"

#Autosome set for allele 1, ben male. 1-A gives freq of allele 2, ben female (common sex)
rare.male.aut.freq <- ggplot(simple.maut, aes(y=(1-A), x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~females, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Female - A ben Frequency")+labs("dominance")+ggtitle("'A' Frequency: Rare Male, RD = 0.5")
rare.male.aut.freq


# Fixed at 1 was for allele 1, beneficial for male. Need allele beneficial for female. If allele 1 was deleted, then allele 2 fixed. 
# So looking at prob allele 1 was del, allele 2 fixed"
rare.male.aut.fix <- ggplot(simple.maut, aes(y=(delA), x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~females, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Female - Prob. ben A Fixed")+labs("Number of Common Sex (Female)")+ggtitle("Prob A Fixed: Rare Male, RD = 0.5")
rare.male.aut.fix



#### Violin Plot


maut<-mdat[mdat$rd == 0.5,] #Autosome
#maut<-maut[maut$h != 99,]

maut$females<-as.factor(maut$females)
maut$h <- as.factor(maut$h)
maut$s <- as.factor(maut$s)
maut$OSR<-as.factor(maut$OSR)


#freq of common sex (female) beneficial allele, allele 2. 
#A set up as allele 1, ben male so 1-A - allele 2
# rare.male.aut.freq.v.all <- ggplot(maut, aes(y=(1-A), x=OSR))+ 
#   scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
#   facet_grid(h~females, labeller=label_both)+theme_bw()+ 
#   theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
#   xlab("Operational Sex Ratio")+ylab("Female - A ben Frequency")+labs("dominance")+ggtitle("'A' Frequency: Rare Male, RD = 0.5")+
#   geom_violin(aes(fill = factor(s)))
# rare.male.aut.freq.v.all

rare.male.aut.freq.v.s <- ggplot(maut, aes(y=(1-A), x=OSR))+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~females, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Female - A ben Frequency")+labs("dominance")+ggtitle("'A' Frequency: Rare Male, RD = 0.5")+
  geom_violin(aes(fill=factor(OSR)))
rare.male.aut.freq.v.s



#### Fitness

mfit.aut <- simple.maut[simple.maut$h != 99,]

mfit.aut$females<-as.factor(mfit.aut$females)
mfit.aut$s<-as.factor(mfit.aut$s)
mfit.aut$h<-as.factor(mfit.aut$h)

str(mfit.aut)

#checking that the degree of fitness change for 

sort.mfit <- mfit.aut[order(mfit.aut$malW, mfit.aut$femW),]
#ort.mfit <- mfit.aut[order(mfit.aut$femW, mfit.aut$malW),]
sp <- seq(from=.01, to=.99, length.out=252)
plot(sort.mfit$femW~sp, col = "red", ylab = "Fitness Impact", ylim =c(0:1)) 
points(sort.mfit$malW/1~sp, col="blue")


#in function, male - female, so take the neg of it to have female - male (common - rare)
mfit.aut$OSR<-as.character(mfit.aut$OSR)
mfit.aut$OSR<-as.numeric(mfit.aut$OSR)
rare.male.aut.fit <- ggplot(mfit.aut, aes(y=-(Wdiff), x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~females, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("A Fitness Difference (common - rare)")+labs("Number of Common Sex (Female)")+ggtitle("A Fitness diff: Rare Male, RD = 0.5")
rare.male.aut.fit


mfit.aut$OSR<-as.factor(mfit.aut$OSR)
rare.male.aut.fit.c <- ggplot(mfit.aut, aes(y= femW, x=malW))+ylim(0.5,1)+xlim(0.5,1)+
  geom_line(aes(colour=s), linewidth=1)+ 
  geom_point(aes(shape=OSR, fill=OSR),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(17,18,20:25))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~females, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Male Fitness (rare sex)")+ylab("Female Fitness (common sex)")+ggtitle("Auto Fitness comparison: Rare Male, RD = 0.5")
rare.male.aut.fit.c



#### Rare Male Sex Chr ####

simple.msex<-simple.mdat[simple.mdat$rd == 0.2,] #Sex-linked

simple.msex$females<-as.factor(simple.msex$females)
simple.msex$s<-as.factor(simple.msex$s)
simple.msex$h<-as.factor(simple.msex$h)

str(simple.msex)


#Looking at X chromosome 

# X was based off allele 2, beneficial for female, so no need to adjust
rare.male.sex.freq.X <- ggplot(simple.msex, aes(y=X, x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~females, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Female - X ben Frequency")+labs("Number of Common Sex (Female)")+ggtitle("'X' Frequency: Rare Male, RD = 0.2")
rare.male.sex.freq.X

rare.male.sex.fix.X <- ggplot(simple.msex, aes(y=fixX, x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~females, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Female - Prob. ben X Fixed")+labs("Number of Common Sex (Female)")+ggtitle("Prob X Fixed: Rare Male, RD = 0.2")
rare.male.sex.fix.X


#Looking at Y chromosome 
# Y was based off allele 1, beneficial for male, need to adjust for female
#Females don't have a Y chr so ben allele for females on the Y should have no positive selection
#ben female Y gene only detrimental to male

" Work on these graphs, the Y-axis isn't making sense "
rare.male.sex.freq.Y <- ggplot(simple.msex, aes(y=(1-Y), x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~females, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Female - Y ben Frequency")+labs("Number of Common Sex (Female)")+ggtitle("'Y' Frequency: Rare Male, RD = 0.2")
rare.male.sex.freq.Y

#Y set up for fixation for males ben allele 1, female is common sex so need prob that male Y doesn't fix
rare.male.sex.fix.Y <- ggplot(simple.msex, aes(y=delY, x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~females, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Female - Prob. ben Y Fixed")+labs("Number of Common Sex (Female)")+ggtitle("Prob Y Fixed: Rare Male, RD = 0.2")
rare.male.sex.fix.Y



### Violin Plots


msex<-mdat[mdat$rd == 0.2,] #Sex-linked
#msex<-msex[msex$h != 99,]

msex$females<-as.factor(msex$females)
msex$s<-as.factor(msex$s)
msex$h<-as.factor(msex$h)
msex$OSR<-as.factor(msex$OSR)

# X was based off allele 2, beneficial for female, so don't need to adjust
rare.male.X.freq.v.s <- ggplot(msex, aes(y=(X), x=OSR))+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~females, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Female - X ben Frequency")+labs("dominance")+ggtitle("'X' Frequency: Rare Male, RD = 0.2")+
  geom_violin(aes(fill=factor(OSR)))
rare.male.X.freq.v.s


# Y was based off allele 1, beneficial for male, need to adjust for female
#Females don't have a Y chr so ben allele for females on the Y should have no positive selection
#ben female Y gene only detrimental to male
rare.male.Y.freq.v.s <- ggplot(msex, aes(y=(1-Y), x=OSR))+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~females, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Female - Y ben Frequency")+labs("dominance")+ggtitle("'Y' Frequency: Rare Male, RD = 0.2")+
  geom_violin(aes(fill=factor(OSR)))
rare.male.Y.freq.v.s




#### Fitness

#Exclude h = 99
mfit.sex <- simple.msex[simple.msex$h != 99,]

mfit.sex$females<-as.factor(mfit.sex$females)
mfit.sex$s<-as.factor(mfit.sex$s)
mfit.sex$h<-as.factor(mfit.sex$h)
str(mfit.sex)

#checking that the degree of fitness change is appropriate and follows a similar line for both species
#fitness is very skewed this time!
sort.sex.mfit <- mfit.sex[order(mfit.sex$malW, mfit.sex$femW),]
#sort.sex.mfit <- mfit.sex[order(mfit.sex$femW, mfit.sex$malW),]
sp <- seq(from=.01, to=.99, length.out=252)
plot(sort.sex.mfit$femW~sp, col = "red", ylab = "Fitness Impact", ylim = c(0,1))
points(sort.sex.mfit$malW/1~sp, col="blue")


#Fitness diff o.g. as male - female. To get common sex (fem) - rare sex (mal), mult by -1
mfit.sex$OSR<-as.character(mfit.sex$OSR)
mfit.sex$OSR<-as.numeric(mfit.sex$OSR)
rare.male.sex.fit <- ggplot(mfit.sex, aes(y= -(Wdiff), x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~females, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Sex chr Fitness Difference (common - rare)")+labs("Number of Common Sex (Female)")+ggtitle("'Sex' Fitness: Rare Male, RD = 0.2")
rare.male.sex.fit

mfit.sex$OSR  <-as.factor(mfit.sex$OSR)
rare.male.sex.fit.c <- ggplot(mfit.sex, aes(y= femW, x=malW))+ylim(0.5,1)+xlim(0.5,1)+ 
  geom_line(aes(colour=s), linewidth=1)+ 
  geom_point(aes(shape=OSR, fill=OSR),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(17,18,20:25))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~females, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Male Fitness (rare sex)")+ylab("Female Fitness (common sex)")+ggtitle("Fitness comparison: Rare Male, RD = 0.2")
rare.male.sex.fit.c




#### Rare Female Analysis, Common Male, allele 1 is beneficial for common sex (male) ####


"Rare Female analysis, Common Male, allele 1 is beneficial for common sex"

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


#### Rare Female Autosome ####

#separate by sex-linked or autosome
simple.faut<-simple.fdat[simple.fdat$rd == 0.5,] #Autosome


simple.faut$males<-as.factor(simple.faut$males)
simple.faut$s<-as.factor(simple.faut$s)
simple.faut$h<-as.factor(simple.faut$h)
str(simple.faut)

#Common male so Y-axis uses male beneficial allele

#autosome used allele 1, male ben, so no change necessary
rare.female.aut.freq <- ggplot(simple.faut, aes(y=A, x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~males, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Male ben A freq")+labs("Number of Common Sex (Male)")+ggtitle("'A' Frequency: Rare Female, RD = 0.5")
rare.female.aut.freq

#autosome used allele 1, male ben, so no change necessary
rare.female.aut.fix <- ggplot(simple.faut, aes(y=fixA, x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~males, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Prob. Male ben A Fixed")+labs("Number of Common Sex (Male)")+ggtitle("Prob A Fixed: Rare Female, RD = 0.5")
rare.female.aut.fix


#### Violin

fdat<-read.csv("rare.female.250.iter.csv", as.is=T, header = TRUE)


faut<-fdat[fdat$rd == 0.5,] #Autosome
#faut<-faut[faut$h != 99,]

faut$males<-as.factor(faut$males)
faut$h <- as.factor(faut$h)
faut$s <- as.factor(faut$s)
faut$OSR<-as.factor(faut$OSR)


#freq of common sex (male) beneficial allele, allele 1. 
#A set up as allele 1
# rare.female.aut.freq.v.all <- ggplot(faut, aes(y=(A), x=OSR))+ 
#   scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
#   facet_grid(h~males, labeller=label_both)+theme_bw()+ 
#   theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
#   xlab("Operational Sex Ratio")+ylab("Male - A ben Frequency")+labs("dominance")+ggtitle("'A' Frequency: Rare Female, RD = 0.5")+
#   geom_violin(aes(fill = factor(s)))
# rare.female.aut.freq.v.all

rare.female.aut.freq.v.s <- ggplot(faut, aes(y=(A), x=OSR))+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~males, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Male - A ben Frequency")+labs("dominance")+ggtitle("'A' Frequency: Rare Female, RD = 0.5")+
  geom_violin(aes(fill=factor(OSR)))
rare.female.aut.freq.v.s



#### Fitness

#Exclude h = 99

ffit.aut <- simple.faut[simple.faut$h != 99,]

sort.ffit <- ffit.aut[order(ffit.aut$malW, ffit.aut$femW),]
#sort.ffit <- ffit.aut[order(ffit.aut$femW, ffit.aut$malW),]
sp <- seq(from=.01, to=.99, length.out=252)
plot(sort.ffit$femW~sp, col = "red", ylab = "Fitness Impact", ylim = c(0:1))
points(sort.ffit$malW/1~sp, col="blue")

ffit.aut$males<-as.factor(ffit.aut$males)
ffit.aut$s<-as.factor(ffit.aut$s)
ffit.aut$h<-as.factor(ffit.aut$h)

str(ffit.aut)

#Fit diff set up as male - female so common - rare, no difference needed
ffit.aut$OSR<-as.character(ffit.aut$OSR)
ffit.aut$OSR<-as.numeric(ffit.aut$OSR)
rare.female.aut.fit <- ggplot(ffit.aut, aes(y=Wdiff, x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~males, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("A Fitness Difference (common - rare)")+ggtitle("'A' Fitness: Rare Female, RD = 0.5")
rare.female.aut.fit

ffit.aut$OSR<-as.factor(ffit.aut$OSR)
rare.female.aut.fit.c <- ggplot(ffit.aut, aes(y= malW, x=femW))+ylim(0.5,1)+xlim(0.5,1)+ 
  geom_line(aes(colour=s), linewidth=1)+ 
  geom_point(aes(shape=OSR, fill=OSR),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(17,18,20:25))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~males, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Female Fitness (rare sex)")+ylab("Male Fitness (common sex)")+ggtitle("Fitness comparison: Rare Female, RD = 0.5")
rare.female.aut.fit.c



#### Rare Female Sex Chr ####

simple.fsex<-simple.fdat[simple.fdat$rd == 0.2,] #Sex-linked

simple.fsex$males<-as.factor(simple.fsex$males)
simple.fsex$s <- as.factor(simple.fsex$s)
simple.fsex$h <- as.factor(simple.fsex$h)
str(simple.fsex)

#X chr freq set using allele 2, ben for females, adjust to look at allele 1 freq, ben males
#Looking at X chromosome 
rare.female.sex.freq.X <- ggplot(simple.fsex, aes(y=(1-X), x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~males, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Male ben X Mean Frequency")+labs("Number of Common Sex (Male)")+ggtitle("'X' Frequency: Rare Female, RD = 0.2")
rare.female.sex.freq.X

#Prob that ben X allele for males fixed. Code gave prob of female ben X fixing, prob of female X being deleted is prob of male ben X fixing
rare.female.sex.fix.X <- ggplot(simple.fsex, aes(y=delX, x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~males, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Prob. male ben X Fixed")+labs("Number of Common Sex (Male)")+ggtitle("Prob X Fixed: Rare Female, RD = 0.2")
rare.female.sex.fix.X


#Looking at Y chromosome
#Y chr set up to look at allele 1, ben males, so no change needed
rare.female.sex.freq.Y <- ggplot(simple.fsex, aes(y=Y, x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~males, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("ben male Y Mean Frequency ")+labs("Number of Common Sex (Male)")+ggtitle("'Y' Frequency: Rare Female, RD = 0.2")
rare.female.sex.freq.Y

rare.female.sex.fix.Y <- ggplot(simple.fsex, aes(y=fixY, x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~males, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Prob. ben male Y Fixed")+labs("Number of Common Sex (Male)")+ggtitle("Prob Y Fixed: Rare Female, RD = 0.2")
rare.female.sex.fix.Y




### Violin 

fsex<-fdat[fdat$rd == 0.2,] #Sex-linked
#msex<-msex[msex$h != 99,]

fsex$males<-as.factor(fsex$males)
fsex$h <- as.factor(fsex$h)
fsex$s <- as.factor(fsex$s)
fsex$OSR<-as.factor(fsex$OSR)


# X was based off allele 2, beneficial for female, adjust for male ben allele, allele 1
rare.female.X.freq.v.s <- ggplot(fsex, aes(y=(1-X), x=OSR))+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~males, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Male - X ben Frequency")+labs("dominance")+ggtitle("'X' Frequency: Rare Female, RD = 0.2")+
  geom_violin(aes(fill=factor(OSR)))
rare.female.X.freq.v.s

# Y was based off allele 1, beneficial for male, don't need to adjust
rare.female.Y.freq.v.s <- ggplot(fsex, aes(y=(Y), x=OSR))+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~males, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Male - Y ben Frequency")+labs("dominance")+ggtitle("'Y' Frequency: Rare Female, RD = 0.2")+
  geom_violin(aes(fill=factor(OSR)))
rare.female.Y.freq.v.s




#### Fitness difference
#Exclude h = 99

ffit.sex <- simple.fsex[simple.fsex$h != 99,]

#checking that the degree of fitness change is appropriate and follows a similar line for both species
#fitness is very skewed this time!
sort.f.sex.fit <- ffit.sex[order(ffit.sex$malW, ffit.sex$femW),]
#sort.f.sex.fit <- ffit.sex[order(ffit.sex$femW, ffit.sex$malW),]
sp <- seq(from=.01, to=.99, length.out=252)
plot(sort.f.sex.fit$femW~sp, col = "red", ylab = "Fitness Impact", ylim = c(0:1))
points(sort.f.sex.fit$malW~sp, col="blue")

ffit.sex$males<-as.factor(ffit.sex$males)
ffit.sex$s<-as.factor(ffit.sex$s)
ffit.sex$h<-as.factor(ffit.sex$h)
str(ffit.sex)


#fitness diff calculated o.g. as male - female so no change needed
ffit.sex$OSR<-as.character(ffit.sex$OSR)
ffit.sex$OSR<-as.numeric(ffit.sex$OSR)
rare.female.sex.fit <- ggplot(ffit.sex, aes(y=Wdiff, x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~males, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Fitness Difference (common-rare)")+labs("Number of Common Sex (Female)")+ggtitle("'Sex' Fitness: Rare Female, RD = 0.2")
rare.female.sex.fit

ffit.sex$OSR<-as.factor(ffit.sex$OSR)
rare.female.sex.fit.c <- ggplot(ffit.sex, aes(y= malW, x=femW))+ylim(0.5,1)+xlim(0.5,1)+ 
  geom_line(aes(colour=s), linewidth=1)+ 
  geom_point(aes(shape=OSR, fill=OSR),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.8, option = "viridis")+scale_shape_manual(values=c(25,23,22,21,22,23,24))+scale_fill_viridis_d(begin =0,end=1, option = "viridis")+
  facet_grid(h~males, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Female Fitness (rare sex)")+ylab("Male Fitness (common sex)")+ggtitle("Fitness comparison: Rare Female, RD = 0.2")
rare.female.sex.fit.c




#### Graph Summary ####

#Autosome
#Frequency of beneficial allele graphs are mostly the same, just flipped on h=0 or 1
rare.male.aut.freq
rare.female.aut.freq

#prob of fixation are the same in the graphs, just flipped on the h = 0 and h = 1
rare.male.aut.fix
rare.female.aut.fix




#Sex Chromosome freq and fix 
rare.male.sex.freq.X
rare.male.sex.fix.X

rare.male.sex.freq.Y
rare.male.sex.fix.Y

rare.female.sex.freq.X
rare.female.sex.fix.X

rare.female.sex.freq.Y
rare.female.sex.fix.Y


#Violin
#Autosomes are flipped but very similar
rare.male.aut.freq.v.s
rare.female.aut.freq.v.s

rare.male.X.freq.v.s
rare.female.X.freq.v.s

rare.male.Y.freq.v.s
rare.female.Y.freq.v.s



#Fitness difference
rare.female.sex.fit
rare.male.sex.fit
rare.female.sex.fit.c
rare.male.sex.fit.c


#Autosome graphs are basically equivalent, just flipped on h 0 and 1
rare.female.aut.fit
rare.male.aut.fit
rare.female.aut.fit.c
rare.male.aut.fit.c

#Fitness check
sp <- seq(from=.01, to=.99, length.out=252)

plot(sort.mfit$femW~sp, col = "red", ylab = "Fitness Impact", ylim = c(0:1), main="Rare Male - Auto fitness") 
points(sort.mfit$malW~sp, col="blue")

plot(sort.ffit$femW~sp, col = "red", ylab = "Fitness Impact", ylim = c(0:1), main="Rare Female - Auto fitness")
points(sort.ffit$malW~sp, col="blue")

plot(sort.sex.mfit$femW~sp, col = "red", ylab = "Fitness Impact", ylim = c(0:1), main="Rare Male - Sex fitness")
points(sort.sex.mfit$malW~sp, col="blue")

plot(sort.f.sex.fit$femW~sp, col = "red", ylab = "Fitness Impact", ylim = c(0:1), main="Rare Female - Sex fitness")
points(sort.f.sex.fit$malW~sp, col="blue")


