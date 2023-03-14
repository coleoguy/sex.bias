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

maut<-simple.mdat[simple.mdat$rd == 0.5,] #Autosome

#Exclude h = 99


maut$females<-as.factor(maut$females)
maut$h <- as.factor(maut$h)
maut$s <- as.factor(maut$s)

str(maut)

#Looking at rare male so Y-axis needs to be from the female/common sex"

#Autosome set for allele 1, ben male. 1-A gives freq of allele 2, ben female (common sex)
rare.male.aut.freq.s <- ggplot(maut, aes(y=(1-A), x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~females, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Female - A ben Frequency")+labs("dominance")+ggtitle("'A' Frequency: Rare Male, RD = 0.5")
rare.male.aut.freq.s


# Fixed at 1 was for allele 1, beneficial for male. Need allele beneficial for female. If allele 1 was deleted, then allele 2 fixed. 
# So looking at prob allele 1 was del, allele 2 fixed"
rare.male.aut.ben.s <- ggplot(maut, aes(y=(delA), x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~females, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Female - Prob. ben A Fixed")+labs("Number of Common Sex (Female)")+ggtitle("Prob A Fixed: Rare Male, RD = 0.5")
rare.male.aut.ben.s

#Fitness

mfit.aut <- maut[maut$h != 99,]

mfit.aut$females<-as.factor(mfit.aut$females)
mfit.aut$s<-as.factor(mfit.aut$s)
mfit.aut$h<-as.factor(mfit.aut$h)
mfit.aut$OSR<-as.factor(mfit.aut$OSR)
str(mfit.aut)

#checking that the degree of fitness change for 

sort.mfit <- mfit.aut[order(mfit.aut$malW, mfit.aut$femW),]
#ort.mfit <- mfit.aut[order(mfit.aut$femW, mfit.aut$malW),]
sp <- seq(from=.01, to=.99, length.out=252)
plot(sort.mfit$femW~sp, col = "red", ylab = "Fitness Impact", ylim =c(0:1)) 
points(sort.mfit$malW/1~sp, col="blue")


#in function, male - female, so take the neg of it to have female - male (common - rare)
rare.male.aut.fit.s <- ggplot(mfit.aut, aes(y=-(Wdiff), x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~females, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("A Fitness Difference (common - rare)")+labs("Number of Common Sex (Female)")+ggtitle("A Fitness diff: Rare Male, RD = 0.5")
rare.male.aut.fit.s

rare.male.aut.fit.c <- ggplot(mfit.aut, aes(y= femW, x=malW))+ylim(0.5,1)+xlim(0.5,1)+
  geom_line(aes(colour=s), linewidth=1)+ 
  geom_point(aes(shape=OSR, fill=OSR),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(17,18,20:25))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~females, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Male Fitness (rare sex)")+ylab("Female Fitness (common sex)")+ggtitle("Auto Fitness comparison: Rare Male, RD = 0.5")
rare.male.aut.fit.c



#### Rare Male Sex Chr ####

msex<-simple.mdat[simple.mdat$rd == 0.2,] #Sex-linked

msex$females<-as.factor(msex$females)
msex$s<-as.factor(msex$s)
msex$h<-as.factor(msex$h)

str(msex)


#Looking at X chromosome 

# X was based off allele 2, beneficial for female, so no need to adjust
rare.male.sex.freq.X.s <- ggplot(msex, aes(y=X, x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~females, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Female - X ben Frequency")+labs("Number of Common Sex (Female)")+ggtitle("'X' Frequency: Rare Male, RD = 0.2")
rare.male.sex.freq.X.s

rare.male.sex.ben.X.s <- ggplot(msex, aes(y=fixX, x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~females, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Female - Prob. ben X Fixed")+labs("Number of Common Sex (Female)")+ggtitle("Prob X Fixed: Rare Male, RD = 0.2")
rare.male.sex.ben.X.s


#Looking at Y chromosome 
### given common sex is female and females do not have a Y - it would make no sense to have a female beneficial allele on the Y



#Fitness difference
#Exclude h = 99

mfit.sex <- msex[msex$h != 99,]

mfit.sex$females<-as.factor(mfit.sex$females)
mfit.sex$s<-as.factor(mfit.sex$s)
mfit.sex$h<-as.factor(mfit.sex$h)
mfit.sex$OSR  <-as.factor(mfit.sex$OSR)
str(mfit.sex)

#checking that the degree of fitness change is appropriate and follows a similar line for both species
#fitness is very skewed this time!
sort.sex.mfit <- mfit.sex[order(mfit.sex$malW, mfit.sex$femW),]
#sort.sex.mfit <- mfit.sex[order(mfit.sex$femW, mfit.sex$malW),]
sp <- seq(from=.01, to=.99, length.out=252)
plot(sort.sex.mfit$femW~sp, col = "red", ylab = "Fitness Impact", ylim = c(0,1))
points(sort.sex.mfit$malW/1~sp, col="blue")


#Fitness diff o.g. as male - female. To get common sex (fem) - rare sex (mal), mult by -1
rare.male.sex.fit.s <- ggplot(mfit.sex, aes(y= -(Wdiff), x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~females, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Sex chr Fitness Difference (common - rare)")+labs("Number of Common Sex (Female)")+ggtitle("'Sex' Fitness: Rare Male, RD = 0.2")
rare.male.sex.fit.s

rare.male.sex.fit.c <- ggplot(mfit.sex, aes(y= femW, x=malW))+ 
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
faut<-simple.fdat[simple.fdat$rd == 0.5,] #Autosome


faut$males<-as.factor(faut$males)
faut$s<-as.factor(faut$s)
faut$h<-as.factor(faut$h)
str(faut)

#Common male so Y-axis uses male beneficial allele

#autosome used allele 1, male ben, so no change necessary
rare.female.aut.freq.s <- ggplot(faut, aes(y=A, x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~males, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Male ben A freq")+labs("Number of Common Sex (Male)")+ggtitle("'A' Frequency: Rare Female, RD = 0.5")
rare.female.aut.freq.s

#autosome used allele 1, male ben, so no change necessary
rare.female.aut.ben.s <- ggplot(faut, aes(y=fixA, x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~males, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Prob. Male ben A Fixed")+labs("Number of Common Sex (Male)")+ggtitle("Prob A Fixed: Rare Female, RD = 0.5")
rare.female.aut.ben.s

#Fitness

#Exclude h = 99

ffit.aut <- faut[faut$h != 99,]

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
rare.female.aut.fit.s <- ggplot(ffit.aut, aes(y=Wdiff, x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~males, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("A Fitness Difference (common - rare")+ggtitle("'A' Fitness: Rare Female, RD = 0.5")
rare.female.aut.fit.s



#### Rare Female Sex Chr ####

fsex<-simple.fdat[simple.fdat$rd == 0.2,] #Sex-linked

fsex$males<-as.factor(fsex$males)
fsex$s <- as.factor(fsex$s)
fsex$h <- as.factor(fsex$h)
str(fsex)

#X chr freq set using allele 2, ben for females, adjust to look at allele 1 freq, ben males
#Looking at X chromosome 
rare.female.sex.freq.X.s <- ggplot(fsex, aes(y=(1-X), x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~males, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Male ben X Mean Frequency")+labs("Number of Common Sex (Male)")+ggtitle("'X' Frequency: Rare Female, RD = 0.2")
rare.female.sex.freq.X.s

#Prob that ben X allele for males fixed. Code gave prob of female ben X fixing, prob of female X being deleted is prob of male ben X fixing
rare.female.sex.ben.X.s <- ggplot(fsex, aes(y=delX, x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~males, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Prob. male ben X Fixed")+labs("Number of Common Sex (Male)")+ggtitle("Prob X Fixed: Rare Female, RD = 0.2")
rare.female.sex.ben.X.s


#Looking at Y chromosome
#Y chr set up to look at allele 1, ben males, so no change needed
rare.female.sex.freq.Y.s <- ggplot(fsex, aes(y=Y, x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~males, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("ben male Y Mean Frequency ")+labs("Number of Common Sex (Male)")+ggtitle("'Y' Frequency: Rare Female, RD = 0.2")
rare.female.sex.freq.Y.s

rare.female.sex.ben.Y.s <- ggplot(fsex, aes(y=fixY, x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~males, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Prob. ben male Y Fixed")+labs("Number of Common Sex (Male)")+ggtitle("Prob Y Fixed: Rare Female, RD = 0.2")
rare.female.sex.ben.Y.s


#Fitness difference
#Exclude h = 99

ffit.sex <- fsex[fsex$h != 99,]

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
rare.female.sex.fit.s <- ggplot(ffit.sex, aes(y=Wdiff, x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~males, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Fitness Difference (common-rare)")+labs("Number of Common Sex (Female)")+ggtitle("'Sex' Fitness: Rare Female, RD = 0.2")
rare.female.sex.fit.s



#### Graph Summary ####

#Autosome
#Frequency of beneficial allele graphs are mostly the same, just flipped on h=0 or 1
rare.male.aut.freq.s
rare.female.aut.freq.s

#prob of fixation are the same in the graphs, just flipped on the h = 0 and h = 1
rare.male.aut.ben.s
rare.female.aut.ben.s




#Sex Chromosome - I don't understand what I am seeing...
rare.male.sex.freq.X.s
rare.male.sex.ben.X.s
#rare.male.sex.freq.Y.s
#rare.male.sex.ben.Y.s

rare.female.sex.freq.X.s
rare.female.sex.ben.X.s

rare.female.sex.freq.Y.s
rare.female.sex.ben.Y.s


#Fitness difference
rare.female.sex.fit.s
rare.male.sex.fit.s

#Autosome graphs are equivalent, just flipped on h 0 and 1
rare.female.aut.fit.s
rare.male.aut.fit.s


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


