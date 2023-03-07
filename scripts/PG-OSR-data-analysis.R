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



"Autosome first"

# setup mean results
male.aut.result <- as.data.frame(matrix(NA,0,12))
colnames(male.aut.result) <- c("common.num", "OSR","h","s","gens","A","X","Y", "prob.ben","prob.del", "Male.Fit", "Female.Fit")

#Because it is autosomal - only using A frequency to calculate male and female fitness

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
        male.aut.result[row.num, 7]<- a <-mean(temp.cur$X) #Average freq over simulations
        male.aut.result[row.num, 8]<- a <-mean(temp.cur$Y) #Average freq over simulations
        
        male.aut.result[row.num, 9]<-sum(temp.cur$A == 1)/length(temp.cur$A) #prob beneficial/fixed
        male.aut.result[row.num, 10]<-sum(temp.cur$A == 0)/length(temp.cur$A) #prob deleterious/removed
        

        #things are going to be janky for anything for h = 99, will exclude in later steps when looking at fitness
        #Since rd = 0.5, looking only at fitness using the autosome frequency
        male.aut.result[row.num,11] <- ((a^2) * 1) + ((2*a*(1-a)) * (1+h[k]*s[l])) + ((1-a)^2 * (1+s[l])) #Male fitness, 1.9 is max
        male.aut.result[row.num,12] <- ((a^2) * 1) + ((2*a*(1-a)) * (1/(1+h[k]*s[l]))) + ((1-a)^2 * (1/(1+s[l]))) #Female fitness, 1 is max

        row.num<-row.num + 1
      }
    }
  }
}


#Exclude h = 99

fitness <- male.aut.result[male.aut.result$h != 99,]

#checking that the degree of fitness change is appropriate and follows a similar line for both species
sort.mfit <- fitness[order(fitness$Male.Fit, fitness$Female.Fit),]
s <- seq(from=.01, to=.99, length.out=252)
plot(1/sort.mfit$Female.Fit~s, col = "red", ylab = "Fitness Impact")
points(sort.mfit$Male.Fit/1~s, col="blue")



male.aut.result$common.num<-as.factor(male.aut.result$common.num)
fitness$common.num<-as.factor(fitness$common.num)

str(fitness)


#Only common.num as factor
rare.male.aut.freq <- ggplot(male.aut.result, aes(y=A, x=OSR))+ 
  geom_line(aes(colour=common.num), linewidth=1)+ylim(c(0,1)) + 
  geom_point(aes(shape=common.num, fill=common.num),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Operational Sex Ratio")+ylab("Mean Frequency of A")+labs("Number of Common Sex (Female)")+ggtitle("'A' Frequency: Rare Male, RD = 0.5")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.male.aut.freq


rare.m.fitness.total <- ggplot(fitness, aes(y=Male.Fit, x=Female.Fit))+ 
  geom_line(aes(colour=common.num),linewidth=1) + 
  geom_point(aes(shape=common.num, fill=common.num),stat="identity", position="identity", size=3)+ 
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~s, labeller=label_both)+theme_bw()+ 
  theme(text= element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+
  xlab("Female Fitness")+ylab("Male Fitness")+labs("Number of Common Sex (Female)")+ggtitle("Autosome - Rare Male - Male vs Female Fitness" )+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
rare.m.fitness.total


#Look at subset
commonN50 <- fitness[fitness$common.num ==1000,]
commonN50$OSR<-as.factor(commonN50$OSR)
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


ggplot(commonN50, aes(y=Male.Fit, x = Female.Fit)) + 
  geom_point(aes(shape=OSR, color=OSR), size = 3) +
  facet_grid(h~s, labeller=label_both) +
  scale_color_manual(values = cbPalette)+
  scale_shape_manual(values = c(15,16,17,18,19,16,15)) +
  scale_fill_manual(values = cbPalette)+
  ggtitle("1000 Common Sex, fitness comparison")
  
  
#further subset
s.9 <- fitness[fitness$s == 0.9,]
s.9$OSR <- as.factor(s.9$OSR)
ggplot(s.9, aes(y=Male.Fit, x = Female.Fit)) +
  geom_point(aes(shape=OSR, color=OSR), size = 3) +
  facet_grid(h~common.num, labeller=label_both) +
  scale_color_manual(values = cbPalette)+
  scale_shape_manual(values = c(15,16,17,18,19,16,15)) +
  scale_fill_manual(values = cbPalette)

#Fitness difference = male - 1/female ?


fitness$Fit.Diff <- fitness$Male.Fit - (1/fitness$Female.Fit)



 fitness[1,11] * fitness[1,12]





male.aut.result[row.num, 9]<- x <-mean(temp.cur$X) #Average X freq over simulations
#x <- mean(temp.cur$X) 
male.aut.result[row.num, 10]<- y <-mean(temp.cur$Y) #Average Y freq over simulations
#y <- mean(temp.cur$Y)
male.aut.result[row.num, 12]<- (x * y * (1 + s[l])) + (x * (1-y) * (1+h[k]*s[l])) + ((1-x) * y * (1+h[k]*s[l])) + ((1-y) * (1-x) * (1))





#### Other old graph scripts
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

