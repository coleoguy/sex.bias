# fate of sex limited deleterious mutations
# common sex: 50, 100, 500, 1000
# OSR: 0.2, 0.1
# s: 0.01, 0.05, 0.1
library(doMC)
library(evobiR)
source("../functions/functions.sex.lim.faster.R")
iter <- 10000
gens <- 2000
cores <- 12
s <- .3
# s .5 males 500 females 5 produces fixation in both but much higher in males
# s=.5 equal males
# s=.3 females
getNe(males=100, females=20)
getNe(males=33, females=33)
com.male <- runSim(males = 100, females = 20,
                   s = s, iter = iter,
                   gens = gens, cores = cores)

com.fema <- runSim(males = 20, females = 100,
                   s = s, iter = iter,
                   gens = gens, cores = cores)

com.equa <- runSim(males = 33, females = 33,
                   s = s, iter = iter,
                   gens = gens, cores = cores)
mean(com.male[,1000:2000])
mean(com.fema[,1000:2000])
mean(com.equa[,1000:2000])

colFixations <- function(x){
  fixed <- function(z){sum(z==1)}
  apply(x, MARGIN = 2, FUN=fixed)
}
#run1 <- list(com.male, com.fema, com.equa)
plot(colMeans(com.male), col="blue", type="l",
     main=paste(s,"mal=500.fem=5"),ylab="freq",xlab="gen")
lines(colMeans(com.fema), col="red")
lines(colMeans(com.equa), type="l", col="green")
fix.mal <- colFixations(com.male)
fix.fem <- colFixations(com.fema)
fix.equ <- colFixations(com.equa)

fix.2 <- rbind(fix.mal,fix.fem,fix.equ)
write.csv(fix.2, file="fix.2.csv")
plot(fix.mal, col="blue", type="l")
lines(fix.equ,col="green")
lines(fix.fem,col="red")
x<- 1:20000
lm(fix.mal~x)
lm(fix.fem~x)
lm(fix.equ~x)

plot(com.male[1,],type="l",ylim=c(0,1))
plot(com.fema[1,],type="l",ylim=c(0,1))
plot(com.equa[1,],type="l",ylim=c(0,1))

hist(com.male[,1000], xlim=c(0,.1), breaks=100)
hist(com.fema[,1000], xlim=c(0,.1), breaks=100)
lines(density(com.equa[,1000]), xlim=c(0,.1),col="red")

ptm <- proc.time()
source("../functions/functions.sex.lim.R")

# this produces a quicker fitness collapse in one sex
males <- c(500) # i
s <- .5        # j
osr <- .01      # k
iter <- 5000
gens<-1000

# this produces a stable quilibria
males <- c(100) # i
s <- .3        # j
osr <- .2      # k
iter <- 50000
gens<-1000

library(doMC)
registerDoMC(14)
females <- males * osr
males.num <- males
x <- foreach(n=1:iter, .combine = "rbind") %dopar%{
  resultA <- c()
  pop <- makeGenomes(females = females,
                     males = males.num,
                     freqs = c(c(females, 0, 0, 0),
                               c(males.num, 0, 0, 0)))
  for(k in 1:gens){
    resultA[k] <- GetFreq(pop, chrom = "A", allele = 2,
                          females = females, males = males.num)
    pop <- Generation(pop, females = females, males = males.num,
                      rd = .5, h = .5, s = s, impacted.sex = "rare")
  }
  # this puts the final result of our sim together
  resultA
}
xmalcom2 <- x





males.num <- males * osr
females <- males
x <- foreach(n=1:iter, .combine = "rbind") %dopar%{
  resultA <- c()
  pop <- makeGenomes(females = females,
                     males = males.num,
                     freqs = c(c(females, 0, 0, 0),
                               c(males.num, 0, 0, 0)))
  for(k in 1:gens){
    resultA[k] <- GetFreq(pop, chrom = "A", allele = 2,
                          females = females, males = males.num)
    pop <- Generation(pop, females = females, males = males.num,
                      rd = .5, h = .5, s = s, impacted.sex = "rare")
  }
  # this puts the final result of our sim together
  resultA
}
x.fem.com2 <- x

vers2 <- proc.time() - ptm








plot(colMeans(x.fem.com2),col="red",type="l")
lines(colMeans(xmalcom2), col="lightblue")


#



















plot(x[1,], ylim=c(0,1),type="l",lwd=.5, col=rgb(1,0,0,.5))
for(i in 2:100) lines(x[i,], type="l",lwd=.5, col=rgb(1,0,0,.5))
for(i in 1:100) lines(xmalcom[i,], type="l",lwd=.5, , col=rgb(0,0,1,.5))

# mean(xmalcom[,3000])
# 0.05926238, 0.07470297, 0.07838416, 0.0847901
# 0.07321485, 0.06726832, 0.06026436, 0.07322277

# mean(x.fem.com[,3000])
# 0.02214554, 0.02835842, 0.02915149, 0.02779307
# 0.01996931, 0.02558317, 0.01716139, 0.02223069
