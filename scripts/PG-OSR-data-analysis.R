"I have fixed the simulation code

Now I have to analyze it. Instead of working with the previous scripts
I'm going to write a pseudo code and analyze everything step by step!


Old comments
#want to find average X, Y, A and proportion with fixation for each set of data
#we care about h for levels/facets, X/Y/A for vertical axis, and OSR for horizontal axis


"



#Load libraries
library(gridExtra)
library(grid)
library(ggpubr)
library(ggplot2)
library(viridis)


#set values
pop<-c(50, 100, 500, 1000) #population size
rd<-c(0.2, 0.5) #sex-linked or autosome
s<-c(0.1, 0.5, 0.9) #selection pressure
h<-c(0, 0.5, 1, 99) #dominance
osr <- c(1,.8,.6,.4,.2,.1,.05) #operation sex ratio




"Rare Male analysis, Common Female"

dat<-read.csv("rare.male.250.iter.csv", as.is=T, header = TRUE)
#separate by sex-linked or autosome

aut<-dat[dat$rd == 0.5,] #Autosome
sex<-dat[dat$rd == 0.2,] #Sex-linked



#### Male Autosome ####

"Autosome first"

# setup mean results
male.aut.result <- as.data.frame(matrix(NA,0,13))
colnames(male.aut.result) <- c("common.num", "OSR","h","s","gens","A","X","Y", "prob.ben","prob.del", "Male.Fit", "Female.Fit", "Fit.Diff")

#Because it is autosomal - only using A frequency to calculate male and female fitness

#Can I convert this into a autosome function
#And create a sex chr function

row.num <- 1
for(i in 1:4){#pop
  for(j in 1:7){#OSR
    for(k in 1:4){#h
      for(l in 1:3){#s
        
        #Grab a subset of the data in which we want to average
        temp.cur<-aut[aut$females ==pop[i],] #population number of common sex
        temp.cur<-temp.cur[temp.cur$OSR ==osr[j],] #OSR
        temp.cur<-temp.cur[temp.cur$h == h[k],] #dominance factor
        temp.cur<-temp.cur[temp.cur$s == s[l],] #selection pressure
        
        #Subset information
        male.aut.result[row.num, 1]<-pop[i] #Population of common sex into col 1
        male.aut.result[row.num, 2]<-osr[j] #OSR into col 2
        male.aut.result[row.num, 3]<-h[k] #Dominance factor into col 3
        male.aut.result[row.num, 4]<-s[l] #Selection pressure into col 4
        
        #Averages
        male.aut.result[row.num, 5]<-mean(temp.cur$gens) #Average number of generations over the simulations
        male.aut.result[row.num, 6]<- a <-mean(temp.cur$A) #Average freq over simulations
        male.aut.result[row.num, 7]<- x <-mean(temp.cur$X) #Average freq over simulations
        male.aut.result[row.num, 8]<- y <-mean(temp.cur$Y) #Average freq over simulations
        
        male.aut.result[row.num, 9]<-sum(temp.cur$A == 1)/length(temp.cur$A) #prob beneficial/fixed
        male.aut.result[row.num, 10]<-sum(temp.cur$A == 0)/length(temp.cur$A) #prob deleterious/removed
        

        #things are going to be janky for anything for h = 99, will exclude in later steps when looking at fitness
        #Since rd = 0.5, looking only at fitness using the autosome frequency
        male.aut.result[row.num,11] <- ((a^2) * 1) + ((2*a*(1-a)) * (1+h[k]*s[l])) + ((1-a)^2 * (1+s[l])) #Male fitness, 1.9 is max
        male.aut.result[row.num,12] <- ((a^2) * 1) + ((2*a*(1-a)) * (1/(1+h[k]*s[l]))) + ((1-a)^2 * (1/(1+s[l]))) #Female fitness, 1 is max
        male.aut.result[row.num,13] <- ((male.aut.result[row.num,11]/(1+male.aut.result[row.num,4])) - (male.aut.result[row.num,12]))
        
        row.num<-row.num + 1
      }
    }
  }
}


#Exclude h = 99

fitness <- male.aut.result[male.aut.result$h != 99,]

#checking that the degree of fitness change is appropriate and follows a similar line for both species
sort.mfit <- fitness[order(fitness$Male.Fit, fitness$Female.Fit),]
sp <- seq(from=.01, to=.99, length.out=252)
plot(1/sort.mfit$Female.Fit~sp, col = "red", ylab = "Fitness Impact")
points(sort.mfit$Male.Fit/1~sp, col="blue")



male.aut.result$common.num<-as.factor(male.aut.result$common.num)
fitness$common.num<-as.factor(fitness$common.num)

str(male.aut.result)
str(fitness)

male.aut.result$h <- as.factor(male.aut.result$h)
male.aut.result$s <- as.factor(male.aut.result$s)


#Only common.num as factor
rare.male.aut.freq <- ggplot(male.aut.result, aes(y=A, x=OSR))+ 
  geom_line(aes(colour=s), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=s, fill=s),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~common.num, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Mean Frequency of A")+labs("dominance")+ggtitle("'A' Frequency: Rare Male, RD = 0.5")
rare.male.aut.freq

rare.male.aut.ben <- ggplot(male.aut.result, aes(y=prob.ben, x=OSR))+ 
  geom_line(aes(colour=common.num), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=common.num, fill=common.num),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Prob. Fixed")+labs("Number of Common Sex (Female)")+ggtitle("Prob A Fixed: Rare Male, RD = 0.5")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.male.aut.ben


rare.male.aut.fit <- ggplot(fitness, aes(y=Fit.Diff, x=OSR))+ 
  geom_line(aes(colour=common.num), linewidth=1)+ 
  geom_point(aes(shape=common.num, fill=common.num),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Fitness Difference")+labs("Number of Common Sex (Female)")+ggtitle("'A' Fitness: Rare Male, RD = 0.5")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.male.aut.fit




#### Male Sex Chr ####

" Sex chr"

# setup mean results
male.sex.result <- as.data.frame(matrix(NA,0,15))
colnames(male.sex.result) <- c("common.num", "OSR","h","s","gens","A","X","Y", "prob.X.ben","prob.X.del","prob.Y.ben","prob.Y.del", "Male.Fit", "Female.Fit", "Fit.Diff")

#Because it is autosomal - only using A frequency to calculate male and female fitness

#Can I convert this into a autosome function
#And create a sex chr function

i<- j <- k <- l <- 1

row.num <- 1
for(i in 1:4){#pop
  for(j in 1:7){#OSR
    for(k in 1:4){#h
      for(l in 1:3){#s
        
        #Grab a subset of the data in which we want to average
        temp.cur<-sex[sex$females ==pop[i],] #population number of common sex
        temp.cur<-temp.cur[temp.cur$OSR ==osr[j],] #OSR
        temp.cur<-temp.cur[temp.cur$h == h[k],] #dominance factor
        temp.cur<-temp.cur[temp.cur$s == s[l],] #selection pressure
        
        #Subset information
        male.sex.result[row.num, 1]<-pop[i] #Population of common sex into col 1
        male.sex.result[row.num, 2]<-osr[j] #OSR into col 2
        male.sex.result[row.num, 3]<-h[k] #Dominance factor into col 3
        male.sex.result[row.num, 4]<-s[l] #Selection pressure into col 4
        
        #Averages
        male.sex.result[row.num, 5]<-mean(temp.cur$gens) #Average number of generations over the simulations
        male.sex.result[row.num, 6]<- a <-mean(temp.cur$A) #Average freq over simulations
        male.sex.result[row.num, 7]<- x <-mean(temp.cur$X) #Average freq over simulations
        male.sex.result[row.num, 8]<- y <-mean(temp.cur$Y) #Average freq over simulations
        
        male.sex.result[row.num, 9]<-sum(temp.cur$X == 1)/length(temp.cur$X) #prob X beneficial/fixed
        male.sex.result[row.num, 10]<-sum(temp.cur$X == 0)/length(temp.cur$X) #prob X deleterious/removed
        
        male.sex.result[row.num, 11]<-sum(temp.cur$Y == 1)/length(temp.cur$Y) #prob Y beneficial/fixed
        male.sex.result[row.num, 12]<-sum(temp.cur$Y == 0)/length(temp.cur$Y) #prob Y deleterious/removed
        
        
        #things are going to be janky for anything for h = 99, will exclude in later steps when looking at fitness
        #Since rd = 0.5, looking only at fitness using the autosome frequency
        male.sex.result[row.num,13] <- ((x*y) * 1) + ((x*(1-y)) * (1+h[k]*s[l]))+ (((1-x)*y) * (1+h[k]*s[l])) + ((1-y)*(1-x) * (1+s[l])) #Male fitness, 1.9 is max
        male.sex.result[row.num,14] <- ((x^2) * 1) + ((2*x*(1-x)) * (1/(1+h[k]*s[l]))) + ((1-x)^2 * (1/(1+s[l]))) #Female fitness, 1 is max
        male.sex.result[row.num,15] <- ((male.sex.result[row.num,13]/(1+male.sex.result[row.num,4])) - (male.sex.result[row.num,14]))
        
        row.num<-row.num + 1
      }
    }
  }
}


#Exclude h = 99

fitness.male.Sex <- male.sex.result[male.sex.result$h != 99,]

#checking that the degree of fitness change is appropriate and follows a similar line for both species
#fitness is very skewed this time!
sort.mfit <- fitness.male.Sex[order(fitness.male.Sex$Male.Fit, fitness.male.Sex$Female.Fit),]
sort.mfit <- fitness.male.Sex[order(fitness.male.Sex$Female.Fit, fitness.male.Sex$Male.Fit),]
sp <- seq(from=.01, to=.99, length.out=252)
plot(1/sort.mfit$Female.Fit~sp, col = "red", ylab = "Fitness Impact")
points(sort.mfit$Male.Fit/1~sp, col="blue")



male.sex.result$common.num<-as.factor(male.sex.result$common.num)
fitness.male.Sex$common.num<-as.factor(fitness.male.Sex$common.num)

str(fitness.male.Sex)


#Looking at X chromosome first
rare.male.sex.freq.X <- ggplot(male.sex.result, aes(y=X, x=OSR))+ 
  geom_line(aes(colour=common.num), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=common.num, fill=common.num),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("X Mean Frequency")+labs("Number of Common Sex (Female)")+ggtitle("'X' Frequency: Rare Male, RD = 0.2")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.male.sex.freq.X

rare.male.sex.ben.X <- ggplot(male.sex.result, aes(y=prob.X.ben, x=OSR))+ 
  geom_line(aes(colour=common.num), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=common.num, fill=common.num),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Prob. Fixed")+labs("Number of Common Sex (Female)")+ggtitle("Prob X Fixed: Rare Male, RD = 0.2")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.male.sex.ben.X




#Looking at Y chromosome first
rare.male.sex.freq.Y <- ggplot(male.sex.result, aes(y=Y, x=OSR))+ 
  geom_line(aes(colour=common.num), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=common.num, fill=common.num),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Y Mean Frequency ")+labs("Number of Common Sex (Female)")+ggtitle("'Y' Frequency: Rare Male, RD = 0.2")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.male.sex.freq.Y

rare.male.sex.ben.Y <- ggplot(male.sex.result, aes(y=prob.Y.ben, x=OSR))+ 
  geom_line(aes(colour=common.num), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=common.num, fill=common.num),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Prob. Fixed")+labs("Number of Common Sex (Female)")+ggtitle("Prob Y Fixed: Rare Male, RD = 0.2")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.male.sex.ben.Y




#Fitness difference
rare.male.sex.fit <- ggplot(fitness.male.Sex, aes(y=Fit.Diff, x=OSR))+ 
  geom_line(aes(colour=common.num), linewidth=1)+ 
  geom_point(aes(shape=common.num, fill=common.num),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Fitness Difference")+labs("Number of Common Sex (Female)")+ggtitle("'Sex' Fitness: Rare Male, RD = 0.2")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.male.sex.fit








#### Female Autosome ####


#set values
pop<-c(50, 100, 500, 1000) #population size
rd<-c(0.2, 0.5) #sex-linked or autosome
s<-c(0.1, 0.5, 0.9) #selection pressure
h<-c(0, 0.5, 1, 99) #dominance
osr <- c(1,.8,.6,.4,.2,.1,.05) #operation sex ratio



"Rare Female analysis, Common Male"

fdat<-read.csv("rare.female.250.iter.csv", as.is=T, header = TRUE)
#separate by sex-linked or autosome

faut<-fdat[fdat$rd == 0.5,] #Autosome
fsex<-fdat[fdat$rd == 0.2,] #Sex-linked




"Autosome first"

# setup mean results
fm.aut.result <- as.data.frame(matrix(NA,0,13))
colnames(fm.aut.result) <- c("common.num", "OSR","h","s","gens","A","X","Y", "prob.ben","prob.del", "Male.Fit", "Female.Fit", "Fit.Diff")

#Because it is autosomal - only using A frequency to calculate male and female fitness

#Can I convert this into a autosome function
#And create a sex chr function

row.num <- 1
for(i in 1:4){#pop
  for(j in 1:7){#OSR
    for(k in 1:4){#h
      for(l in 1:3){#s
        
        #Grab a subset of the data in which we want to average
        temp.cur<-faut[faut$females ==pop[i],] #population number of common sex 
        ####GO BACK AND fix females -> males
        temp.cur<-temp.cur[temp.cur$OSR ==osr[j],] #OSR
        temp.cur<-temp.cur[temp.cur$h == h[k],] #dominance factor
        temp.cur<-temp.cur[temp.cur$s == s[l],] #selection pressure
        
        #Subset information
        fm.aut.result[row.num, 1]<-pop[i] #Population of common sex into col 1
        fm.aut.result[row.num, 2]<-osr[j] #OSR into col 2
        fm.aut.result[row.num, 3]<-h[k] #Dominance factor into col 3
        fm.aut.result[row.num, 4]<-s[l] #Selection pressure into col 4
        
        #Averages
        fm.aut.result[row.num, 5]<-mean(temp.cur$gens) #Average number of generations over the simulations
        fm.aut.result[row.num, 6]<- a <-mean(temp.cur$A) #Average freq over simulations
        fm.aut.result[row.num, 7]<- x <-mean(temp.cur$X) #Average freq over simulations
        fm.aut.result[row.num, 8]<- y <-mean(temp.cur$Y) #Average freq over simulations
        
        fm.aut.result[row.num, 9]<-sum(temp.cur$A == 1)/length(temp.cur$A) #prob beneficial/fixed
        fm.aut.result[row.num, 10]<-sum(temp.cur$A == 0)/length(temp.cur$A) #prob deleterious/removed
        
        
        #things are going to be janky for anything for h = 99, will exclude in later steps when looking at fitness
        #Since rd = 0.5, looking only at fitness using the autosome frequency
        fm.aut.result[row.num,11] <- ((a^2) * 1) + ((2*a*(1-a)) * (1+h[k]*s[l])) + ((1-a)^2 * (1+s[l])) #Male fitness, 1.9 is max
        fm.aut.result[row.num,12] <- ((a^2) * 1) + ((2*a*(1-a)) * (1/(1+h[k]*s[l]))) + ((1-a)^2 * (1/(1+s[l]))) #Female fitness, 1 is max
        fm.aut.result[row.num,13] <- ((fm.aut.result[row.num,12]- (fm.aut.result[row.num,11]/(1+fm.aut.result[row.num,4]))))
        
        row.num<-row.num + 1
      }
    }
  }
}


#Exclude h = 99

FAfitness <- fm.aut.result[fm.aut.result$h != 99,]

#checking that the degree of fitness change is appropriate and follows a similar line for both species
sort.mfit <- FAfitness[order(FAfitness$Male.Fit, FAfitness$Female.Fit),]
sp <- seq(from=.01, to=.99, length.out=252)
plot(1/sort.mfit$Female.Fit~sp, col = "red", ylab = "Fitness Impact")
points(sort.mfit$Male.Fit/1~sp, col="blue")



fm.aut.result$common.num<-as.factor(fm.aut.result$common.num)
FAfitness$common.num<-as.factor(FAfitness$common.num)

str(FAfitness)


#Only common.num as factor
rare.female.aut.freq <- ggplot(fm.aut.result, aes(y=A, x=OSR))+ 
  geom_line(aes(colour=common.num), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=common.num, fill=common.num),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Mean Frequency of A")+labs("Number of Common Sex (Male)")+ggtitle("'A' Frequency: Rare Female, RD = 0.5")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.female.aut.freq

rare.female.aut.ben <- ggplot(fm.aut.result, aes(y=prob.ben, x=OSR))+ 
  geom_line(aes(colour=common.num), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=common.num, fill=common.num),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Prob. Fixed")+labs("Number of Common Sex (Male)")+ggtitle("Prob A Fixed: Rare Female, RD = 0.5")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.female.aut.ben


rare.female.aut.fit <- ggplot(FAfitness, aes(y=Fit.Diff, x=OSR))+ 
  geom_line(aes(colour=common.num), linewidth=1)+ 
  geom_point(aes(shape=common.num, fill=common.num),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Fitness Difference")+labs("Number of Common Sex (Male)")+ggtitle("'A' Fitness: Rare Female, RD = 0.5")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.female.aut.fit




#### Feale Sex Chr ####

" Sex chr"

# setup mean results
fm.sex.result <- as.data.frame(matrix(NA,0,15))
colnames(fm.sex.result) <- c("common.num", "OSR","h","s","gens","A","X","Y", "prob.X.ben","prob.X.del","prob.Y.ben","prob.Y.del", "Male.Fit", "Female.Fit", "Fit.Diff")

#Because it is autosomal - only using A frequency to calculate male and female fitness

#Can I convert this into a autosome function
#And create a sex chr function

i<- j <- k <- l <- 1

row.num <- 1
for(i in 1:4){#pop
  for(j in 1:7){#OSR
    for(k in 1:4){#h
      for(l in 1:3){#s
        
        #Grab a subset of the data in which we want to average
        temp.cur<-sex[sex$females ==pop[i],] #population number of common sex
        temp.cur<-temp.cur[temp.cur$OSR ==osr[j],] #OSR
        temp.cur<-temp.cur[temp.cur$h == h[k],] #dominance factor
        temp.cur<-temp.cur[temp.cur$s == s[l],] #selection pressure
        
        #Subset information
        fm.sex.result[row.num, 1]<-pop[i] #Population of common sex into col 1
        fm.sex.result[row.num, 2]<-osr[j] #OSR into col 2
        fm.sex.result[row.num, 3]<-h[k] #Dominance factor into col 3
        fm.sex.result[row.num, 4]<-s[l] #Selection pressure into col 4
        
        #Averages
        fm.sex.result[row.num, 5]<-mean(temp.cur$gens) #Average number of generations over the simulations
        fm.sex.result[row.num, 6]<- a <-mean(temp.cur$A) #Average freq over simulations
        fm.sex.result[row.num, 7]<- x <-mean(temp.cur$X) #Average freq over simulations
        fm.sex.result[row.num, 8]<- y <-mean(temp.cur$Y) #Average freq over simulations
        
        fm.sex.result[row.num, 9]<-sum(temp.cur$X == 1)/length(temp.cur$X) #prob X beneficial/fixed
        fm.sex.result[row.num, 10]<-sum(temp.cur$X == 0)/length(temp.cur$X) #prob X deleterious/removed
        
        fm.sex.result[row.num, 11]<-sum(temp.cur$Y == 1)/length(temp.cur$Y) #prob Y beneficial/fixed
        fm.sex.result[row.num, 12]<-sum(temp.cur$Y == 0)/length(temp.cur$Y) #prob Y deleterious/removed
        
        
        #things are going to be janky for anything for h = 99, will exclude in later steps when looking at fitness
        #Since rd = 0.5, looking only at fitness using the autosome frequency
        fm.sex.result[row.num,13] <- ((x*y) * 1) + ((x*(1-y)) * (1+h[k]*s[l]))+ (((1-x)*y) * (1+h[k]*s[l])) + ((1-y)*(1-x) * (1+s[l])) #Male fitness, 1.9 is max
        fm.sex.result[row.num,14] <- ((x^2) * 1) + ((2*x*(1-x)) * (1/(1+h[k]*s[l]))) + ((1-x)^2 * (1/(1+s[l]))) #Female fitness, 1 is max
        fm.sex.result[row.num,15] <- ((fm.sex.result[row.num,14])-(fm.sex.result[row.num,13]/(1+fm.sex.result[row.num,4])))
        
        row.num<-row.num + 1
      }
    }
  }
}


#Exclude h = 99

fitness.female.Sex <- fm.sex.result[fm.sex.result$h != 99,]

#checking that the degree of fitness change is appropriate and follows a similar line for both species
#fitness is very skewed this time!
sort.mfit <- fitness.female.Sex[order(fitness.female.Sex$Male.Fit, fitness.female.Sex$Female.Fit),]
sort.mfit <- fitness.female.Sex[order(fitness.female.Sex$Female.Fit, fitness.female.Sex$Male.Fit),]
sp <- seq(from=.01, to=.99, length.out=252)
plot(1/sort.mfit$Female.Fit~sp, col = "red", ylab = "Fitness Impact")
points(sort.mfit$Male.Fit/1~sp, col="blue")



fm.sex.result$common.num<-as.factor(fm.sex.result$common.num)
fitness.female.Sex$common.num<-as.factor(fitness.female.Sex$common.num)

str(fitness.female.Sex)


#Looking at X chromosome first
rare.female.sex.freq.X <- ggplot(fm.sex.result, aes(y=X, x=OSR))+ 
  geom_line(aes(colour=common.num), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=common.num, fill=common.num),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("X Mean Frequency")+labs("Number of Common Sex (Male)")+ggtitle("'X' Frequency: Rare Female, RD = 0.2")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.female.sex.freq.X

rare.female.sex.ben.X <- ggplot(fm.sex.result, aes(y=prob.X.ben, x=OSR))+ 
  geom_line(aes(colour=common.num), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=common.num, fill=common.num),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Prob. Fixed")+labs("Number of Common Sex (Male)")+ggtitle("Prob X Fixed: Rare Female, RD = 0.2")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.female.sex.ben.X




#Looking at Y chromosome first
rare.female.sex.freq.Y <- ggplot(fm.sex.result, aes(y=Y, x=OSR))+ 
  geom_line(aes(colour=common.num), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=common.num, fill=common.num),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Y Mean Frequency ")+labs("Number of Common Sex (Male)")+ggtitle("'Y' Frequency: Rare Female, RD = 0.2")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.female.sex.freq.Y

rare.female.sex.ben.Y <- ggplot(fm.sex.result, aes(y=prob.Y.ben, x=OSR))+ 
  geom_line(aes(colour=common.num), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=common.num, fill=common.num),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Prob. Fixed")+labs("Number of Common Sex (Male)")+ggtitle("Prob Y Fixed: Rare Female, RD = 0.2")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.female.sex.ben.Y




#Fitness difference
rare.female.sex.fit <- ggplot(fitness.female.Sex, aes(y=Fit.Diff, x=OSR))+ 
  geom_line(aes(colour=common.num), linewidth=1)+ 
  geom_point(aes(shape=common.num, fill=common.num),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Fitness Difference")+labs("Number of Common Sex (Female)")+ggtitle("'Sex' Fitness: Rare Female, RD = 0.2")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.female.sex.fit





