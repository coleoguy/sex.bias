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



#### Rare Female analysis, Common Female, allele 1 is beneficial for common sex (male) ####

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



#### Graphs ####
#Autosomes are flipped but very similar
rare.male.aut.freq.v.s
rare.female.aut.freq.v.s

rare.male.X.freq.v.s
rare.female.X.freq.v.s

rare.female.Y.freq.v.s
rare.male.Y.freq.v.s