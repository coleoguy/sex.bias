
load("~/Desktop/Dropbox/gitrepos/sex.bias/results/sexant.RData")
results <- list(autosomes.1, autosomes.5, autosomes.9,
                sexchroms.1, sexchroms.5, sexchroms.9)
# add sex fitness
# function to calculate mean fitness in each sex
fits <- function(x){
  h <- x$h
  s <- x$s
  if(x$rd == .5){
    p <- x$A
    fgenos <- mgenos <- c(p^2, 2*p*(1-p), (1-p)^2)
  }
  if(x$rd < .5){
    px <- x$X
    py <- x$Y
    fgenos <- c(px^2, 2*px*(1-px), (1-px)^2)
    mgenos <- c(px*py, 
                px*(1-py) + (1-px)*py, 
                (1-px)*(1-py))
  }
  if(h != 99){
    malh <- h
    femh <- 1-h
    males <- c(1+s, 1+h*s, 1)
    females <- c(1, 1+femh*s, 1+s)
  }
  if(h == 99){
    males <- c(1+s, 1+s, 1)
    females <- c(1, 1+s, 1+s)
  }
  females <- females/(1+s)
  males <- males/(1+s)
  return(c(sum(mgenos*males), sum(fgenos*females)))
}

# use fits function to calculate absolute mean fitness of each sex
for(j in 1:6){
  male.fit <- fem.fit <- vector(length=nrow(results[[j]]))
  for(i in 1:nrow(results[[j]])){
    print(i)
    curfit <- fits(results[[j]][i,])
    male.fit[i] <- curfit[1]
    fem.fit[i] <- curfit[2]
  }
  # combine results
  results[[j]] <- data.frame(results[[j]],male.fit,fem.fit)
}
# clean up the space
rm(curfit, male.fit, fem.fit, i, fits, j, foo)

# setup mean results
mean.results <- as.data.frame(matrix(NA,0,11))
colnames(mean.results) <- c("common.num", "OSR","h","s","gens","mal.fit","fem.fit","rd","freq","prob.ben","prob.del")
com.nums <- c(50,100,500,1000)
osrs <- c(1,.8,.6,.4,.2,.1,.05)
hs <- c(0,.5,1,"SSD")
row.counter <- 1
ss <- c(.1,.5,.9,.1,.5,.9)
rds <- c(.5,.5,.5,.2,.2,.2)
for(ii in 1:6){
  cur <- results[[ii]]
  for(i in 1:4){ #com.nums
    for(j in 1:7){ # OSRs
      for(k in 1:3){ # h
          temp.cur <- cur[cur$males == com.nums[i], ]
          temp.cur <- temp.cur[temp.cur$OSR == osrs[j], ]
          temp.cur <- temp.cur[temp.cur$h == hs[k], ]
          temp.cur <- temp.cur[temp.cur$s == ss[ii], ]
          mean.results[row.counter, 1] <- com.nums[i]
          mean.results[row.counter, 2] <- osrs[j]
          mean.results[row.counter, 3] <- hs[k]
          mean.results[row.counter, 4] <- ss[ii]
          mean.results[row.counter, 5] <- mean(temp.cur$gens)
          mean.results[row.counter, 6] <- mean(temp.cur$male.fit)
          mean.results[row.counter, 7] <- mean(temp.cur$fem.fit)
          mean.results[row.counter, 9] <- mean(temp.cur$A)
          mean.results[row.counter, 10] <- sum(temp.cur$A == 1)/length(temp.cur$A)
          mean.results[row.counter, 11] <- sum(temp.cur$A == 0)/length(temp.cur$A)
          mean.results[row.counter, 8] <- rds[ii]
          row.counter <- row.counter + 1
          print(row.counter)
      }
    }
  }
}

fem.means <- data.frame(mean.results[,-6], rep("rare", 504))
mal.means <- data.frame(mean.results[,-7], rep("common", 504))
colnames(fem.means)[c(6,11)] <- colnames(mal.means)[c(6,11)] <- c("fitness", "sex")
means.all <- rbind(fem.means, mal.means)
means.all.auto<- means.all[means.all$rd==.5,]
means.all.sex<- means.all[means.all$rd==.2,]
means.all.auto.9 <- means.all.auto[means.all.auto$s==.9,]
means.all.auto.5 <- means.all.auto[means.all.auto$s==.5,]
means.all.auto.1 <- means.all.auto[means.all.auto$s==.1,]

GetNe <- function(x){
  Nc <- x$common.num
  Nr <- x$common.num * x$OSR
  Ne <- (4*Nc*Nr) / (Nc + Nr)
  return(Ne)
}

Ne.1 <- Ne.5 <- Ne.9 <- c()


for(i in 1:nrow(means.all.auto.1)){
  Ne.1[i] <- GetNe(means.all.auto.1[i,]) 
  Ne.5[i] <- GetNe(means.all.auto.5[i,]) 
  Ne.9[i] <- GetNe(means.all.auto.9[i,]) 
}
means.all.auto.1 <- data.frame(means.all.auto.1, Ne.1)
means.all.auto.5 <- data.frame(means.all.auto.5, Ne.5)
means.all.auto.9 <- data.frame(means.all.auto.9, Ne.9)
colnames(means.all.auto.1)[12] <- 
  colnames(means.all.auto.5)[12] <- 
  colnames(means.all.auto.9)[12] <- "Ne"

foo <- means.all.auto.9
foo$Ne <- round(foo$Ne)
ggplot(foo, aes(y=fitness, x=OSR)) + 
  geom_point(aes(shape=as.factor(sex), colour=Ne), 
             stat="identity", position="identity", alpha=0.7, size=2) + 
  scale_color_viridis_c(trans="log", breaks=c(10,50,250,1000),name=expression(paste(N[e])))+
  facet_grid(h ~ common.num) + theme_bw() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", 
                          size=13, hjust=0.5, vjust=0.5)) + 
  scale_size(range=c(1, 4)) + guides(shape=guide_legend(title="sex")) + 
  xlab("Operational Sex Ratio") + ylab("Absolute Fitness")
# export 4"x8"


foo.freq <- data.frame(foo[,-c(9,10)],"frequency")
foo.prob.ben <- data.frame(foo[,-c(8,10)],"prob.ben")
foo.prob.del <- data.frame(foo[,-c(8,9)],"prob.del")
colnames(foo.freq)[c(8,11)] <- colnames(foo.prob.ben)[c(8,11)] <- 
  colnames(foo.prob.del)[c(8,11)] <- c("value","type")
foo2 <- rbind(foo.freq, foo.prob.ben, foo.prob.del)
ggplot(foo2, aes(y=value, x=OSR)) + 
  geom_point(aes(colour=type), 
             stat="identity", position="identity", alpha=0.35, size=2) + 
  scale_color_viridis_d()+
  facet_grid(h ~ common.num) + theme_bw() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", 
                          size=13, hjust=0.5, vjust=0.5)) + 
  scale_size(range=c(1, 4)) + guides(shape=guide_legend(title="sex")) + 
  xlab("Operational Sex Ratio") + ylab("Male benefit allele")




foo3 <- foo2[foo2$type=="frequency", ]
fit.means <- foo3$fitness[foo3$sex=="common"] - foo3$fitness[foo3$sex=="rare"]
Nes <- foo3$Ne[foo3$sex=="common"]
foo3 <- data.frame(fit.means, Nes, 
                   foo3$common.num[foo3$sex=="common"],
                   foo3$h[foo3$sex=="common"],
                   foo3$OSR[foo3$sex=="common"])
colnames(foo3) <- c("Fitness", "Ne","commmon", "h","OSR")
gen.avgd <- as.data.frame(matrix(NA,84,4))
counter <-1
for(i in seq(from=1, by=3, length.out = 84)){
  gen.avgd[counter, 1] <- mean(foo3[i:(i+2), 1])
  gen.avgd[counter, 2:4] <- foo3[i,c(2:3,5)]
  counter<- counter+1
}
colnames(gen.avgd) <- c("fitness", "Ne","common","OSR")
gen.avgd<- gen.avgd[1:28,]
gen.avgd$common <- factor(gen.avgd$common, levels=c("50","100","500","1000"))
#gen.avgd$OSR <- factor(gen.avgd$OSR, levels=c("0,05","0.10", "0.20","0.40","0.60","0.80","1.00"))


ggplot(gen.avgd, aes(y=fitness, x=OSR)) + 
  geom_point(aes(colour=common), stat="identity", position="identity", alpha=1, size=3) + 
  geom_line(aes(colour=common), stat="identity", position="identity", alpha=1) + theme_bw() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  scale_size(range=c(1, 3)) + xlab("Operational Sex Ratio") + ylab("fitness difference (common - rare)")
