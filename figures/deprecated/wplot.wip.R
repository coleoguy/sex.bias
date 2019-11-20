# 
# category <- c("a", "b", "c", "d", "e", "f", "g")
# wbar.m <- c(0.8, 0.83, 0.835, 0.85, 0.867, 0.872, 0.85)
# wbar.f <- c(0.78, 0.817, 0.82, 0.8, 0.79, 0.76, 0.75)
# max.w <- rep((1), 7)
# category2 <- c("A", "B", "C", "D", "E", "F", "G")
# rm(wbar.plot)
# wbar.plot <- cbind(category, wbar.m, category2, wbar.f, max.w)
# wbar.plot <- as.data.frame(wbar.plot)
# 
# ggplot(wbar.plot, aes(x=category, y=max.w)) +
#   geom_segment(aes(x=category, xend=category, y=wbar.m, yend=max.w), 
#                size=1, data=wbar.plot, colour="black", linetype="solid", alpha= 0.5) +
#   geom_segment(aes(x=category2, xend=category2, y=wbar.f, yend=max.w), 
#                size=1, data=wbar.plot, colour="blue", linetype="solid", alpha= 0.5) +
#   geom_point(aes(x= category, y=wbar.m), data= wbar.plot, size=4, colour= "black", alpha = 0.5) +
#   geom_point(aes(x= category, y=max.w), data= wbar.plot, size=4, colour= "black", alpha = 0.5) +
#   geom_point(aes(x= category2, y=wbar.f), data= wbar.plot, size=4, colour= "blue", alpha = 0.5) +
#   geom_point(aes(x= category2, y=max.w), data= wbar.plot, size=4, colour= "blue", alpha = 0.5) 
# 
# category <- c("OSR1", "OSR.8", "OSR.6", "OSR.4", "OSR.2", "OSR.1", "OSR.05")
# wbar.m <- c(0.8, 0.83, 0.835, 0.85, 0.867, 0.872, 0.85)
# wbar.f <- c(0.78, 0.817, 0.82, 0.8, 0.79, 0.76, 0.75)
# max.w <- rep((1), 7)
# category2 <- c("osr1", "osr.8", "osr.6", "osr.4", "osr.2", "osr.1", "osr.05")
# rm(wbar.plot)
# wbar.plot <- cbind(category, wbar.m, category2, wbar.f, max.w)
# wbar.plot <- as.data.frame(wbar.plot)


library(viridis)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)

# First let's parse the results from the XY model simulations
# We want to calculate the mean fitness for each combination of parameters
# of interest (e.g. dominance factor, selection coefficient, rd) and fill
# each cell for a given Common.sex:OSR condition in separate tables for males
# and females.


source("../../sex.bias/rscripts/functions/analysis.functions.R")
load("../../sex.bias/results/base.comp.model.RData")
res.base <- results
load("../../sex.bias/results/rare.male.model.RData")
rare.male <- results
load("../../sex.bias/results/rare.female.model.RData")
rare.female <- results
rm(results)

# create a data structure to store results
results.fit.mal <- matrix(, 7, 4)

# insert column and row names
colnames(results.fit.mal) <- c("50", "100", 
                               "500", "1000")
rownames(results.fit.mal) <- c("OSR1","OSR.8","OSR.6",
                               "OSR.4","OSR.2","OSR.1",
                               "OSR.05")

# create a data structure to store results
results.fit.fem <- matrix(, 7, 4)

# insert column and row names
colnames(results.fit.fem) <- c("50", "100", 
                               "500", "1000")
rownames(results.fit.fem) <- c("OSR1","OSR.8","OSR.6",
                               "OSR.4","OSR.2","OSR.1",
                               "OSR.05")

h <- 0.5
s <- 0.5
loc <- "auto"

for(i in 1:4){
  # cycle through OSR values
  for(j in 1:7){
    if(j == 1){
      # looks at base model
      fit <- getFitness(data = res.base[[i]]$rd0.5$h0.5$s0.5, loc, h, s)
      results.fit.mal[j, i] <- mean(fit$wbar.m)
    }
    if(j > 1){
      # looks at bias model
      fit <-  getFitness(data = rare.female[[i]][[j-1]]$rd0.5$h0.5$s0.5, loc, h, s)
      results.fit.mal[j, i] <- mean(fit$wbar.m)
    }
  }
}

for(i in 1:4){
  # cycle through OSR values
  for(j in 1:7){
    if(j == 1){
      # looks at base model
      fit <- getFitness(data = res.base[[i]]$rd0.5$h0.5$s0.5, loc, h, s)
      results.fit.fem[j, i] <- mean(fit$wbar.f)
    }
    if(j > 1){
      # looks at bias model
      fit <-  getFitness(data = rare.female[[i]][[j-1]]$rd0.5$h0.5$s0.5, loc, h, s)
      results.fit.fem[j, i] <- mean(fit$wbar.f)
    }
  }
}

# Now, let's reshape the data to make it long and add a column with maximum fitness
# values for each sex based on the fitness functions we used in our model.

wbar.fem.h05 <- melt(results.fit.fem[1:7,])
wbar.mal.h05 <- melt(results.fit.mal[1:7,])
max.w <- rep(c(1), 28)
category <- c("a", "b", "c", "d", "e", "f", "g")
category2 <- c("A", "B", "C", "D", "E", "F", "G")

wbar.h05 <- cbind(wbar.fem.h05, wbar.mal.h05, max.w)
colnames(wbar.h05) <- c("OSR1", "Common.sex", "wbar.f", "OSR2", "Common.sex", "wbar.m", "max.w")

wbar.plot <- wbar.h05[(wbar.h05$Common.sex == 100),]

wbar.plot <- wbar.plot[,-c(1:2)]
wbar.plot <- wbar.plot[,-c(2:3)]

# category <- c("OSR1", "OSR.8", "OSR.6", "OSR.4", "OSR.2", "OSR.1", "OSR.05")
# category2 <- c("osr1", "osr.8", "osr.6", "osr.4", "osr.2", "osr.1", "osr.05")

wbar <- cbind(wbar.plot, category, category2)
# wbar.plot$category <- factor(wbar.plot$category, 
#                              levels = c("OSR1", "OSR.8", "OSR.6", "OSR.4", "OSR.2", "OSR.1", "OSR.05"))
# wbar.plot$category2 <- factor(wbar.plot$category2, 
#                               levels = c("osr1", "osr.8", "osr.6", "osr.4", "osr.2", "osr.1", "osr.05"))


# order <- c(1, 2, 3, 4, 5, 6, 7)
# wbar.plot <- cbind(wbar.plot, order)
# wbar.plot$category <- factor(wbar.plot$category, levels = wbar.plot$category[order(wbar.plot$order)])
# wbar.plot$category2 <- factor(wbar.plot$category2, levels = wbar.plot$category2[order(wbar.plot$order)])

male.c <- ggplot(wbar, aes(x=category, y=max.w)) +
  geom_segment(aes(x=category, xend=category, y=wbar.m, yend=max.w), 
               size=1, data=wbar, colour="#29788E", linetype="solid", alpha= 0.35) +
  geom_segment(aes(x=category2, xend=category2, y=wbar.f, yend=max.w), 
               size=1, data=wbar, colour="#440154", linetype="solid", alpha= 0.35) +
  geom_point(aes(x= category, y=wbar.m), data= wbar, size=4, colour= "#29788E", alpha = 0.8) +
  geom_point(aes(x= category, y=max.w), data= wbar, size=4, colour= "#29788E", alpha = 0.8) +
  geom_point(aes(x= category2, y=wbar.f), data= wbar, size=4, colour= "#440154", alpha = 0.8) +
  geom_point(aes(x= category2, y=max.w), data= wbar, size=4, colour= "#440154", alpha = 0.8) +
  theme_light() + theme(axis.ticks.x=element_blank()) + ylim(0.7, 1.05) +
  xlab("OSR") + ylab("Mean fitness") + ggtitle("Mean fitness (Common sex = males)")



  
# wbar.plot <- wbar.plot[,-c(2, 3)]
# wbar.plot <- cbind(wbar.plot, category, category2)
# wbar.plot$category <- factor(wbar.plot$category, levels = wbar.plot$category[order(wbar.plot$wbar.m)])
# wbar.plot$category2 <- factor(wbar.plot$category2, levels = wbar.plot$category2[order(wbar.plot$wbar.f)])

# Now Sex chromosomes. Females = common

# create a data structure to store results
results.fit.mal <- matrix(, 7, 4)

# insert column and row names
colnames(results.fit.mal) <- c("50", "100", 
                               "500", "1000")
rownames(results.fit.mal) <- c("OSR1","OSR.8","OSR.6",
                               "OSR.4","OSR.2","OSR.1",
                               "OSR.05")

# create a data structure to store results
results.fit.fem <- matrix(, 7, 4)

# insert column and row names
colnames(results.fit.fem) <- c("50", "100", 
                               "500", "1000")
rownames(results.fit.fem) <- c("OSR1","OSR.8","OSR.6",
                               "OSR.4","OSR.2","OSR.1",
                               "OSR.05")

loc <- "sex"

for(i in 1:4){
  # cycle through OSR values
  for(j in 1:7){
    if(j == 1){
      # looks at base model
      fit <- getFitness(data = res.base[[i]]$rd0.1$h0.5$s0.5, loc, h, s)
      results.fit.mal[j, i] <- mean(fit$wbar.m)
    }
    if(j > 1){
      # looks at bias model
      fit <-  getFitness(data = rare.male[[i]][[j-1]]$rd0.1$h0.5$s0.5, loc, h, s)
      results.fit.mal[j, i] <- mean(fit$wbar.m)
    }
  }
}

for(i in 1:4){
  # cycle through OSR values
  for(j in 1:7){
    if(j == 1){
      # looks at base model
      fit <- getFitness(data = res.base[[i]]$rd0.1$h0.5$s0.5, loc, h, s)
      results.fit.fem[j, i] <- mean(fit$wbar.f)
    }
    if(j > 1){
      # looks at bias model
      fit <-  getFitness(data = rare.male[[i]][[j-1]]$rd0.1$h0.5$s0.5, loc, h, s)
      results.fit.fem[j, i] <- mean(fit$wbar.f)
    }
  }
}

wbar.fem.h05 <- melt(results.fit.fem[1:7,])
wbar.mal.h05 <- melt(results.fit.mal[1:7,])
max.w <- rep(c(1), 28)
category <- c("a", "b", "c", "d", "e", "f", "g")
category2 <- c("A", "B", "C", "D", "E", "F", "G")

wbar.h05 <- cbind(wbar.fem.h05, wbar.mal.h05, max.w)
colnames(wbar.h05) <- c("OSR1", "Common.sex", "wbar.f", "OSR2", "Common.sex", "wbar.m", "max.w")

wbar.plot <- wbar.h05[(wbar.h05$Common.sex == 100),]

wbar.plot <- wbar.plot[,-c(1:2)]
wbar.plot <- wbar.plot[,-c(2:3)]

wbar <- cbind(wbar.plot, category, category2)

fem.sexchrom <- ggplot(wbar, aes(x=category, y=max.w)) +
  geom_segment(aes(x=category, xend=category, y=wbar.m, yend=max.w), 
               size=1, data=wbar, colour="#29788E", linetype="solid", alpha= 0.35) +
  geom_segment(aes(x=category2, xend=category2, y=wbar.f, yend=max.w), 
               size=1, data=wbar, colour="#440154", linetype="solid", alpha= 0.35) +
  geom_point(aes(x= category, y=wbar.m), data= wbar, size=4, colour= "#29788E", alpha = 0.8) +
  geom_point(aes(x= category, y=max.w), data= wbar, size=4, colour= "#29788E", alpha = 0.8) +
  geom_point(aes(x= category2, y=wbar.f), data= wbar, size=4, colour= "#440154", alpha = 0.8) +
  geom_point(aes(x= category2, y=max.w), data= wbar, size=4, colour= "#440154", alpha = 0.8) +
  theme_light() + theme(axis.ticks.x=element_blank()) + ylim(0.7, 1.05) +
  xlab("OSR") + ylab("Mean fitness") + ggtitle("S.A.L on sex chromosome. (Common sex = females)")


# Now Sex chromosomes. Males = common

# create a data structure to store results
results.fit.mal <- matrix(, 7, 4)

# insert column and row names
colnames(results.fit.mal) <- c("50", "100", 
                               "500", "1000")
rownames(results.fit.mal) <- c("OSR1","OSR.8","OSR.6",
                               "OSR.4","OSR.2","OSR.1",
                               "OSR.05")

# create a data structure to store results
results.fit.fem <- matrix(, 7, 4)

# insert column and row names
colnames(results.fit.fem) <- c("50", "100", 
                               "500", "1000")
rownames(results.fit.fem) <- c("OSR1","OSR.8","OSR.6",
                               "OSR.4","OSR.2","OSR.1",
                               "OSR.05")

loc <- "sex"

for(i in 1:4){
  # cycle through OSR values
  for(j in 1:7){
    if(j == 1){
      # looks at base model
      fit <- getFitness(data = res.base[[i]]$rd0.1$h0.5$s0.5, loc, h, s)
      results.fit.mal[j, i] <- mean(fit$wbar.m)
    }
    if(j > 1){
      # looks at bias model
      fit <-  getFitness(data = rare.female[[i]][[j-1]]$rd0.1$h0.5$s0.5, loc, h, s)
      results.fit.mal[j, i] <- mean(fit$wbar.m)
    }
  }
}

for(i in 1:4){
  # cycle through OSR values
  for(j in 1:7){
    if(j == 1){
      # looks at base model
      fit <- getFitness(data = res.base[[i]]$rd0.1$h0.5$s0.5, loc, h, s)
      results.fit.fem[j, i] <- mean(fit$wbar.f)
    }
    if(j > 1){
      # looks at bias model
      fit <-  getFitness(data = rare.female[[i]][[j-1]]$rd0.1$h0.5$s0.5, loc, h, s)
      results.fit.fem[j, i] <- mean(fit$wbar.f)
    }
  }
}

wbar.fem.h05 <- melt(results.fit.fem[1:7,])
wbar.mal.h05 <- melt(results.fit.mal[1:7,])
max.w <- rep(c(1), 28)
category <- c("a", "b", "c", "d", "e", "f", "g")
category2 <- c("A", "B", "C", "D", "E", "F", "G")

wbar.h05 <- cbind(wbar.fem.h05, wbar.mal.h05, max.w)
colnames(wbar.h05) <- c("OSR1", "Common.sex", "wbar.f", "OSR2", "Common.sex", "wbar.m", "max.w")

wbar.plot <- wbar.h05[(wbar.h05$Common.sex == 100),]

wbar.plot <- wbar.plot[,-c(1:2)]
wbar.plot <- wbar.plot[,-c(2:3)]

wbar <- cbind(wbar.plot, category, category2)
osr <- c("OSR1", "OSR.8", "OSR.6", "OSR.4", "OSR.2", "OSR.1", "OSR.05")
wbar <- cbind(wbar, osr)

mal.sexchrom <- ggplot(wbar, aes(x=category, y=max.w)) +
  geom_segment(aes(x=category, xend=category, y=wbar.m, yend=max.w), 
               size=1, data=wbar, colour="#29788E", linetype="solid", alpha= 0.35) +
  geom_segment(aes(x=category2, xend=category2, y=wbar.f, yend=max.w), 
               size=1, data=wbar, colour="#440154", linetype="solid", alpha= 0.35) +
  geom_point(aes(x= category, y=wbar.m), 
             data= wbar, size=4, colour= "#29788E", alpha = 0.8) +  
  geom_point(aes(x= category, y=max.w), 
             data= wbar, size=4, colour= "#29788E", alpha = 0.8, shape = wbar.m) +
  geom_point(aes(x= category2, y=wbar.f), 
             data= wbar, size=4, colour= "#440154", alpha = 0.8) +
  geom_point(aes(x= category2, y=max.w), 
             data= wbar, size=4, colour= "#440154", alpha = 0.8, shape=wbar.f) +
  theme_light() + theme(axis.ticks.x=element_blank()) + ylim(0.7, 1.05) +
  xlab("OSR") + ylab("Mean fitness") + ggtitle("S.A.L on sex chromosome. (Common sex = males)") +
  theme(legend.position="bottom") 

library(cowplot)
g1 <- ggplotGrob(male.c)
g2 <- ggplotGrob(fem.sexchrom)
g3 <- ggplotGrob(mal.sexchrom)


plot_grid(
  g1, g2, g3, ncol=3,
  labels = c('A)', 'B)', 'C)'),
  align="hv"
)

library(gridExtra)
grid.arrange(g1, g2, g3, ncol = 3)


