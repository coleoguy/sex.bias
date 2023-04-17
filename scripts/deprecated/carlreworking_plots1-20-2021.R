###Used to generate viridis plots for OSR paper
##need to load data from HB.autosomes
library(gridExtra)
library(grid)
library(ggpubr)
library(ggplot2)
library(viridis)
####rare male data and plots####
dat<-read.csv("rare.male.250.iter.csv", as.is=T, header = TRUE)
str(dat)

#dat$females<-as.factor(dat$females)
#dat$h<-as.factor(dat$h)
#dat$rd<-as.factor(dat$rd)
#dat$s<-as.factor(dat$s)
#dat$OSR<-as.factor(dat$OSR)
#levels(dat$females)
#levels(dat$h)
#levels(dat$rd)
#levels(dat$s)
#levels(dat$OSR)
#want to find average X, Y, A and proportion with fixation for each set of data
#we care about h for levels/facets, X/Y/A for vertical axis, and OSR for horizontal axis

#females has 50, 100, 500, 1000
#h has 0, 0.5, 1, and 99
#rd has 0.2 and 0.5
#s has 0.1, 0.5, and 0.9
#osr has 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1
#we want to separate based on rds
pop<-c(50, 100, 500, 1000)
rd<-c(0.2, 0.5)
s<-c(0.1, 0.5, 0.9)
h<-c(0, 0.5, 1, 99)
osr <- c(1,.8,.6,.4,.2,.1,.05)
# setup mean results
male.aut.result <- as.data.frame(matrix(NA,0,12))
colnames(male.aut.result) <- c("common.num", "OSR","h","s","gens","mean.freq", "prob.ben","prob.del", "X", "X.fit", "Y", "Y.fit")
#com.nums <- c(50,100,500,1000)

#On autosome or linked to sex chr
aut<-dat[dat$rd == 0.5,]
sex<-dat[dat$rd == 0.2,]
row.num<-1
str(sex)


#First 0.5 (autosome)
for(i in 1:4){#pop
  for(j in 1:7){#OSR
    for(k in 1:4){#h
      for(l in 1:3){#s
        temp.cur<-aut[aut$females ==pop[i],] #population number of common sex
        temp.cur<-temp.cur[temp.cur$OSR ==osr[j],] #OSR
        temp.cur<-temp.cur[temp.cur$h == h[k],] #dominance factor
        temp.cur<-temp.cur[temp.cur$s == s[l],] #selection pressure
        male.aut.result[row.num, 1]<-pop[i] #Population of common sex into col 1
        male.aut.result[row.num, 2]<-osr[j] #OSR into col 2
        male.aut.result[row.num, 3]<-h[k] #Dominance factor into col 3
        male.aut.result[row.num, 4]<-s[l] #Selection pressure into col 4
        male.aut.result[row.num, 5]<-mean(temp.cur$gens) #Average number of generations over the simulations
        male.aut.result[row.num, 6]<-mean(temp.cur$A) #Average freq over simulations
        male.aut.result[row.num, 7]<-sum(temp.cur$A == 1)/length(temp.cur$A) #prob beneficial ?
        male.aut.result[row.num, 8]<-sum(temp.cur$A == 0)/length(temp.cur$A) #prob deleted?
        
        male.aut.result[row.num, 9]<-mean(temp.cur$X) #Average X freq over simulations
        x <- mean(temp.cur$X)
        male.aut.result[row.num,10] <- (x^2 * 1) + (2*x*(1-x) * (1+h[k]*s[l])) + ((1-x)^2 * (1+s[l]))
        
        male.aut.result[row.num, 11]<-mean(temp.cur$Y) #Average Y freq over simulations
        y <- mean(temp.cur$Y)
        male.aut.result[row.num, 12]<- (x * y * (1 + s[l])) + (x * (1-y) * (1+h[k]*s[l])) + ((1-x) * y * (1+h[k]*s[l])) + ((1-y) * (1-x) * (1))
        
        row.num<-row.num + 1
      }
    }
  }
}

male.aut.result$common.num<-as.factor(male.aut.result$common.num)
str(male.aut.result)


fitness.male.aut <- male.aut.result[male.aut.result$h != 99,]


#h<-c(0, 0.5, 1)

fit.male.y <-ggplot(fitness.male.aut, aes(y=Y.fit, x=OSR))+
  geom_line(aes(colour=common.num), size=1)+
  geom_point(aes(shape=common.num, fill=common.num), stat="identity", position="identity", size=3)+
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller=label_both)+theme_bw()+ theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Mean Frequency of X")+labs("Number of Common Sex (Female)")+ggtitle("'X' Frequency: Rare Male, RD 0.2")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
fit.male.y



#####rare male autosome frequency plot####
rare.male.aut.freq <- ggplot(male.aut.result, aes(y=mean.freq, x=OSR))+ 
  geom_line(aes(colour=common.num), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=common.num, fill=common.num),
             stat="identity", position="identity", size=3)+ scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Mean Frequency of A")+labs("Number of Common Sex (Female)")+ggtitle("'A' Frequency: Rare Male, RD = 0.5")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.male.aut.freq

####rare male autosome prob plot####
rare.male.aut.prob<-ggplot(male.aut.result, aes(y=prob.ben, x=OSR))+geom_line(aes(colour=common.num), linewidth=1)+ylim(c(0,1))+
  geom_point(aes(shape=common.num, fill=common.num), stat="identity", position="identity", size=3)+
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller=label_both)+theme_bw()+ theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Proportion with fixed 1 for A")+labs("Number of Common Sex (Female)")+ggtitle("Proportion Fixed A: Rare Male, RD = 0.5")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.male.aut.prob

####now to look with rd 0.2 (sex chr) ####


male.x.result <- as.data.frame(matrix(NA,0,12))
colnames(male.x.result) <- c("common.num", "OSR","h","s","gens","mean.freq", "prob.ben","prob.del", "X", "X.fit", "Y", "Y.fit")
sex<-dat[dat$rd == 0.2,]
row.num<-1
str(sex)
for(i in 1:4){#pop
  for(j in 1:7){#OSR
    for(k in 1:4){#h
      for(l in 1:3){#s
        temp.cur<-sex[sex$females ==pop[i],]
        temp.cur<-temp.cur[temp.cur$OSR ==osr[j],]
        temp.cur<-temp.cur[temp.cur$h == h[k],]
        temp.cur<-temp.cur[temp.cur$s == s[l],]
        male.x.result[row.num, 1]<-pop[i]
        male.x.result[row.num, 2]<-osr[j]
        male.x.result[row.num, 3]<-h[k]
        male.x.result[row.num, 4]<-s[l]
        male.x.result[row.num, 5]<-mean(temp.cur$gens)
        male.x.result[row.num, 6]<-mean(temp.cur$X)
        male.x.result[row.num, 7]<-sum(temp.cur$X == 1)/length(temp.cur$X)
        male.x.result[row.num, 8]<-sum(temp.cur$X == 0)/length(temp.cur$X)
        male.x.result[row.num, 9]<-mean(temp.cur$X) #Average X freq over simulations
        x <- mean(temp.cur$X)
        male.x.result[row.num,10] <- (x^2 * 1) + (2*x*(1-x) * (1+h[k]*s[l])) + ((1-x)^2 * (1+s[l]))
        
        male.x.result[row.num, 11]<-mean(temp.cur$Y) #Average Y freq over simulations
        y <- mean(temp.cur$Y)
        male.x.result[row.num, 12]<- (x * y * (1 + s[l])) + (x * (1-y) * (1+h[k]*s[l])) + ((1-x) * y * (1+h[k]*s[l])) + ((1-y) * (1-x) * (1))
        
        row.num<-row.num + 1
      }
    }
  }
}

fitness.male.x <- male.x.result[male.x.result$h != 99,]
fit.male <-ggplot(male.x.result, aes(y=Y.fit, x=OSR))+
  geom_line(aes(colour=common.num), size=1)+
  geom_point(aes(shape=common.num, fill=common.num), stat="identity", position="identity", size=3)+
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller=label_both)+theme_bw()+ theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Mean Frequency of X")+labs("Number of Common Sex (Female)")+ggtitle("'X' Frequency: Rare Male, RD 0.2")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.male.sex.x.freq


#ggplot(male.x.result, aes(x=OSR, y = mean.freq)) +geom_point(aes(shape=common.num, fill=common.num) + facet_grid(h~s, labeller=label_both)


str(male.x.result)
male.x.result$common.num<-as.factor(male.x.result$common.num)
####rare male sex x frequency plot####
rare.male.sex.x.freq<-ggplot(male.x.result, aes(y=mean.freq, x=OSR))+
  geom_line(aes(colour=common.num), size=1)+ylim(c(0,1))+
  geom_point(aes(shape=common.num, fill=common.num), stat="identity", position="identity", size=3)+
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller=label_both)+theme_bw()+ theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Mean Frequency of X")+labs("Number of Common Sex (Female)")+ggtitle("'X' Frequency: Rare Male, RD 0.2")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.male.sex.x.freq

####rare male sex x prob plot####
rare.male.sex.x.prob<-ggplot(male.x.result, aes(y=prob.ben, x=OSR))+geom_line(aes(colour=common.num), size=1)+ylim(c(0,1))+
  geom_point(aes(shape=common.num, fill=common.num), stat="identity", position="identity", size=3)+
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller=label_both)+theme_bw()+ theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Proportion with fixed 1 for X")+labs("Number of Common Sex (Female)")+ggtitle("Proportion Fixed X: Rare Male, RD = 0.2")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.male.sex.x.prob


####rare female data and plots####

datfem<-read.csv("rare.female.250.iter.csv", as.is=T)
str(datfem)

#want to find average X, Y, A and proportion with fixation for each set of data
#we care about h for levels/facets, X/Y/A for vertical axis, and OSR for horizontal axis

#males has 50, 100, 500, 1000
#h has 0, 0.5, 1, and 99
#rd has 0.2 and 0.5
#s has 0.1, 0.5, and 0.9
#osr has 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1
#we want to separate based on rds
pop<-c(50, 100, 500, 1000)
rds<-c(0.2, 0.5)
ss<-c(0.1, 0.5, 0.9)
hs<-c(0, 0.5, 1, 99)
osrs <- c(1,.8,.6,.4,.2,.1,.05)
# setup mean results
mean.results.rare.female.aut <- as.data.frame(matrix(NA,0,8))
colnames(mean.results.rare.female.aut) <- c("common.num", "OSR","h","s","gens","mean.freq", "prob.ben","prob.del")
#com.nums <- c(50,100,500,1000)
femaut<-datfem[datfem$rd == 0.5,]
femsex<-datfem[datfem$rd == 0.2,]
#aut.9<-femdat[(femdat$rd == 0.5) & (femdat$s == 0.9),]
row.num<-1
#str(sex)
i<-j<-k<-l<-1
for(i in 1:4){#pop
  for(j in 1:7){#OSRs
    for(k in 1:4){#hs
      for(l in 1:3){#ss
        temp.cur<-femaut[femaut$females ==pop[i],]
        temp.cur<-temp.cur[temp.cur$OSR ==osrs[j],]
        temp.cur<-temp.cur[temp.cur$h == hs[k],]
        temp.cur<-temp.cur[temp.cur$s == ss[l],]
        mean.results.rare.female.aut[row.num, 1]<-pop[i]
        mean.results.rare.female.aut[row.num, 2]<-osrs[j]
        mean.results.rare.female.aut[row.num, 3]<-hs[k]
        mean.results.rare.female.aut[row.num, 4]<-ss[l]
        mean.results.rare.female.aut[row.num, 5]<-mean(temp.cur$gens)
        mean.results.rare.female.aut[row.num, 6]<-mean(temp.cur$A)
        mean.results.rare.female.aut[row.num, 7]<-sum(temp.cur$A == 1)/length(temp.cur$A)
        mean.results.rare.female.aut[row.num, 8]<-sum(temp.cur$A == 0)/length(temp.cur$A)
        row.num<-row.num + 1
      }
    }
  }
}

str(mean.results.rare.female.aut)
mean.results.rare.female.aut$common.num<-as.factor(mean.results.rare.female.aut$common.num)

#####rare female plots aut freq####
rare.female.aut.freq<-ggplot(mean.results.rare.female.aut, aes(y=mean.freq, x=OSR))+geom_line(aes(colour=common.num), linewidth=1)+ylim(c(0,1))+
  geom_point(aes(shape=common.num, fill=common.num), stat="identity", position="identity", size=3)+
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller = label_both)+theme_bw()+ theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Mean Frequency of A")+labs("Number of Common Sex (Male)")+ggtitle("'A' Frequency: Rare Female, RD = 0.5")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.female.aut.freq

####rare female plots aut prob####
rare.female.aut.prob<-ggplot(mean.results.rare.female.aut, aes(y=prob.ben, x=OSR))+geom_line(aes(colour=common.num), linewidth=1)+ylim(c(0,1))+
  geom_point(aes(shape=common.num, fill=common.num), stat="identity", position="identity", size=3)+
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller=label_both)+theme_bw()+ theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Proportion with fixed 1 for A")+labs("Number of Common Sex (Male)")+ggtitle("Proportion Fixed A: Rare Female, RD = 0.5")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.female.aut.prob

####now to look with rd 0.2 for rare females with x####
mean.results.rare.female.sex.x <- as.data.frame(matrix(NA,0,8))
colnames(mean.results.rare.female.sex.x) <- c("common.num", "OSR","h","s","gens","mean.freq", "prob.ben","prob.del")
femsex<-datfem[datfem$rd == 0.2,]
row.num<-1
i<-j<-k<-l<-1
str(femsex)
for(i in 1:4){#pop
  for(j in 1:7){#OSRs
    for(k in 1:4){#hs
      for(l in 1:3){#ss
        temp.cur<-femsex[femsex$females ==pop[i],]
        temp.cur<-temp.cur[temp.cur$OSR ==osrs[j],]
        temp.cur<-temp.cur[temp.cur$h == hs[k],]
        temp.cur<-temp.cur[temp.cur$s == ss[l],]
        mean.results.rare.female.sex.x[row.num, 1]<-pop[i]
        mean.results.rare.female.sex.x[row.num, 2]<-osrs[j]
        mean.results.rare.female.sex.x[row.num, 3]<-hs[k]
        mean.results.rare.female.sex.x[row.num, 4]<-ss[l]
        mean.results.rare.female.sex.x[row.num, 5]<-mean(temp.cur$gens)
        mean.results.rare.female.sex.x[row.num, 6]<-mean(temp.cur$X)
        mean.results.rare.female.sex.x[row.num, 7]<-sum(temp.cur$X == 1)/length(temp.cur$X)
        mean.results.rare.female.sex.x[row.num, 8]<-sum(temp.cur$X == 0)/length(temp.cur$X)
        row.num<-row.num + 1
      }
    }
  }
}

str(mean.results.rare.female.sex.x)
mean.results.rare.female.sex.x$common.num<-as.factor(mean.results.rare.female.sex.x$common.num)


####rare female sex x frequency plot####
rare.female.sex.x.freq<-ggplot(mean.results.rare.female.sex.x, aes(y=mean.freq, x=OSR))+geom_line(aes(colour=common.num), size=1)+ylim(c(0,1))+
  geom_point(aes(shape=common.num, fill=common.num), stat="identity", position="identity", size=3)+
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller = label_both)+theme_bw()+ theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Mean Frequency of X")+labs("Number of Common Sex (Male)")+ggtitle("'X' Frequency: Rare Female, RD = 0.2")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.female.sex.x.freq

####rare female sex x prob plot####
rare.female.sex.x.prob<-ggplot(mean.results.rare.female.sex.x, aes(y=prob.ben, x=OSR))+geom_line(aes(colour=common.num), size=1)+ylim(c(0,1))+
  geom_point(aes(shape=common.num, fill=common.num), stat="identity", position="identity", size=3)+
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller = label_both)+theme_bw()+ theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Proportion with fixed 1 for X")+labs("Number of Common Sex (Male)")+ggtitle("Proportion Fixed X: Rare Female, RD = 0.2")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.female.sex.x.prob


####now to look with rd 0.2 for rare females with y####
mean.results.rare.female.sex.y <- as.data.frame(matrix(NA,0,8))
colnames(mean.results.rare.female.sex.y) <- c("common.num", "OSR","h","s","gens","mean.freq", "prob.ben","prob.del")
femsex<-datfem[datfem$rd == 0.2,]
row.num<-1
i<-j<-k<-l<-1
str(femsex)
for(i in 1:4){#pop
  for(j in 1:7){#OSRs
    for(k in 1:4){#hs
      for(l in 1:3){#ss
        temp.cur<-sex[femsex$males ==pop[i],]
        temp.cur<-temp.cur[temp.cur$OSR ==osrs[j],]
        temp.cur<-temp.cur[temp.cur$h == hs[k],]
        temp.cur<-temp.cur[temp.cur$s == ss[l],]
        mean.results.rare.female.sex.y[row.num, 1]<-pop[i]
        mean.results.rare.female.sex.y[row.num, 2]<-osrs[j]
        mean.results.rare.female.sex.y[row.num, 3]<-hs[k]
        mean.results.rare.female.sex.y[row.num, 4]<-ss[l]
        mean.results.rare.female.sex.y[row.num, 5]<-mean(temp.cur$gens)
        mean.results.rare.female.sex.y[row.num, 6]<-mean(temp.cur$Y)
        mean.results.rare.female.sex.y[row.num, 7]<-sum(temp.cur$Y == 1)/length(temp.cur$Y)
        mean.results.rare.female.sex.y[row.num, 8]<-sum(temp.cur$Y == 0)/length(temp.cur$Y)
        row.num<-row.num + 1
      }
    }
  }
}

str(mean.results.rare.female.sex.y)
mean.results.rare.female.sex.y$common.num<-as.factor(mean.results.rare.female.sex.y$common.num)


####rare female sex y frequency plot####
rare.female.sex.y.freq<-ggplot(mean.results.rare.female.sex.y, aes(y=mean.freq, x=OSR))+geom_line(aes(colour=common.num), size=1)+ylim(c(0,1))+
  geom_point(aes(shape=common.num, fill=common.num), stat="identity", position="identity", size=3)+
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller = label_both)+theme_bw()+ theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Mean Frequency of Y")+labs("Number of Common Sex (Male)")+ggtitle("'Y' frequency: Rare Female, RD = 0.2")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.female.sex.y.freq

####rare female sex y prob plot####
rare.female.sex.y.prob<-ggplot(mean.results.rare.female.sex.y, aes(y=prob.ben, x=OSR))+geom_line(aes(colour=common.num), size=1)+ylim(c(0,1))+
  geom_point(aes(shape=common.num, fill=common.num), stat="identity", position="identity", size=3)+
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller=label_both)+theme_bw()+ theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Proportion with fixed 1 for Y")+labs("Number of Common Sex (Male)")+ggtitle("Proportion Fixed Y: Rare Female, RD = 0.2")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.female.sex.y.prob






####making new plots for paper####

#### I want to use rare females with s = 0.9 and no h =99 and rd = 0.5 (autosome)
datfem<-read.csv("rare.female.250.iter.csv", as.is=T)
str(datfem)


str(datfem)
datfem$h<-as.numeric(datfem$h)
#males has 50, 100, 500, 1000
#h has 0, 0.5, 1, and 99
#rd has 0.2 and 0.5
#s has 0.1, 0.5, and 0.9
#osr has 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1
#we want to separate based on rds
pop<-c(50, 100, 500, 1000)
rds<-c(0.2, 0.5)
ss<-c(0.1, 0.5, 0.9)
hs<-c(0, 0.5, 1)
osrs <- c(1,.8,.6,.4,.2,.1,.05)

# setup mean results
female.aut.dat <- as.data.frame(matrix(NA,0,7))
colnames(female.aut.dat) <- c("common.num", "OSR","h","gens","mean.freq", "prob.ben","prob.del")
#com.nums <- c(50,100,500,1000)

femaut.9<-datfem[(datfem$rd == 0.5) & (datfem$s == 0.9),]
str(femaut.9)
femaut.9.h<-femaut.9[(femaut.9$h != 99),]
#including h options 0,0.5, and 1 still

str(femaut.9.h)
row.num<-1
i<-j<-k<-l<-1
for(i in 1:4){#pop
  for(j in 1:7){#OSRs
    for(k in 1:3){#hs
        temp.cur<-femaut.9.h[femaut.9.h$females ==pop[i],]
        temp.cur<-temp.cur[temp.cur$OSR ==osrs[j],]
        temp.cur<-temp.cur[temp.cur$h == hs[k],]
        female.aut.dat[row.num, 1]<-pop[i]
        female.aut.dat[row.num, 2]<-osrs[j]
        female.aut.dat[row.num, 3]<-hs[k]
        female.aut.dat[row.num, 4]<-mean(temp.cur$gens)
        female.aut.dat[row.num, 5]<-mean(temp.cur$A)
        female.aut.dat[row.num, 6]<-sum(temp.cur$A == 1)/length(temp.cur$A)
        female.aut.dat[row.num, 7]<-sum(temp.cur$A == 0)/length(temp.cur$A)
        row.num<-row.num + 1
      }
    }
  }
str(female.aut.dat)

female.aut.dat$common.num<-as.factor(female.aut.dat$common.num)



g1<-ggplot(female.aut.dat, aes(y=prob.ben, x=OSR))+geom_line(aes(colour=common.num), linewidth=1)+
  geom_point(aes(shape=common.num, fill=common.num), stat="identity", position="identity", size=3)+
  scale_color_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~ .)+theme_bw()+ theme(text=element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+theme(axis.title.y = element_blank())+
  scale_size(range=c(1,4))+xlab("Operational Sex Ratio")+ylab("Proporion of simulations with fixation")+labs("Number of common sex")+ggtitle("Proportion with fixation")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
g1

g2<-ggplot(female.aut.dat, aes(y=mean.freq, x=OSR))+geom_line(aes(colour=common.num), linewidth=1)+
  geom_point(aes(shape=common.num, fill=common.num), stat="identity", position="identity", size=3)+
  scale_color_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~ .)+theme_bw()+ theme(text=element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+theme(axis.title.y = element_blank())+
  scale_size(range=c(1,4))+xlab("Operational Sex Ratio")+ggtitle("Frequency of allele")+ylim(c(0.00,1.00))+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
g2

final<-ggarrange(g1, g2, ncol=2, legend="bottom", common.legend=TRUE)
final
