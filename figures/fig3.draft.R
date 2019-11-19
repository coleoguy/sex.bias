
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

wbar <- cbind(wbar.plot, category, category2)
order <- c(1:7)
wbar <- cbind(wbar, order)

col.vir <- viridis(7)

wbar <- cbind(wbar, col.vir)

auto.mal <- ggplot(wbar, aes(x=as.factor(category), y=max.w)) +
  geom_segment(aes(x= category, xend=category, y=wbar.m, yend=max.w), 
               size=1, data=wbar, colour="red", linetype="solid", alpha= 0.35) +
  geom_segment(aes(x=category2, xend=category2, y=wbar.f, yend=max.w), 
               size=1, data=wbar, colour="black", linetype="solid", alpha= 0.35) +
  geom_point(aes(fill=as.factor(order), x= category, y=wbar.m), data= wbar, size=2, 
             alpha = 1, , colour = col.vir) +
  geom_point(aes(fill=as.factor(order),x= category, y=max.w), data= wbar, size=2, 
             alpha = 1, colour = col.vir) +
  geom_point(aes(fill=as.factor(order), x= category2, y=wbar.f), data= wbar, size=2,
             alpha = 1, colour = col.vir) +
  geom_point(aes(fill=as.factor(order), x= category2, y=max.w), data= wbar, size=2,
             alpha = 1, colour = col.vir) +
  theme_light() + theme(axis.ticks.x=element_blank(), legend.position = 'none') + ylim(0.7, 1.05) +
  xlab("OSR") + ylab("Mean fitness") + ggtitle("Autosomes", subtitle = 'Common sex: Males')

#########
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
col.vir <- viridis(7)

wbar <- cbind(wbar, col.vir)

sexchrom.fem <- ggplot(wbar, aes(x=as.factor(category), y=max.w)) +
  geom_segment(aes(x= category, xend=category, y=wbar.m, yend=max.w), 
               size=1, data=wbar, colour="red", linetype="solid", alpha= 0.35) +
  geom_segment(aes(x=category2, xend=category2, y=wbar.f, yend=max.w), 
               size=1, data=wbar, colour="black", linetype="solid", alpha= 0.35) +
  geom_point(aes(fill=as.factor(order), x= category, y=wbar.m), data= wbar, size=2, 
             alpha = 1, , colour = col.vir) +
  geom_point(aes(fill=as.factor(order),x= category, y=max.w), data= wbar, size=2, 
             alpha = 1, colour = col.vir) +
  geom_point(aes(fill=as.factor(order), x= category2, y=wbar.f), data= wbar, size=2,
             alpha = 1, colour = col.vir) +
  geom_point(aes(fill=as.factor(order), x= category2, y=max.w), data= wbar, size=2,
             alpha = 1, colour = col.vir) +
  theme_light() + theme(axis.ticks.x=element_blank(), legend.position = 'none') + ylim(0.7, 1.05) +
  xlab("OSR") + ylab("Mean fitness") + ggtitle("Sex chromosomes", subtitle="Common sex: Females")


#######
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
col.vir <- viridis(7)

wbar <- cbind(wbar, col.vir)

sexchrom.mal <- ggplot(wbar, aes(x=as.factor(category), y=max.w)) +
  geom_segment(aes(x= category, xend=category, y=wbar.m, yend=max.w), 
               size=1, data=wbar, colour="red", linetype="solid", alpha= 0.35) +
  geom_segment(aes(x=category2, xend=category2, y=wbar.f, yend=max.w), 
               size=1, data=wbar, colour="black", linetype="solid", alpha= 0.35) +
  geom_point(aes(fill=as.factor(order), x= category, y=wbar.m), data= wbar, size=2, 
             alpha = 1, , colour = col.vir) +
  geom_point(aes(fill=as.factor(order),x= category, y=max.w), data= wbar, size=2, 
             alpha = 1, colour = col.vir) +
  geom_point(aes(fill=as.factor(order), x= category2, y=wbar.f), data= wbar, size=2,
             alpha = 1, colour = col.vir) +
  geom_point(aes(fill=as.factor(order), x= category2, y=max.w), data= wbar, size=2,
             alpha = 1, colour = col.vir) +
  theme_light() + theme(axis.ticks.x=element_blank(), legend.position = 'none') + ylim(0.7, 1.05) +
  xlab("OSR") + ylab("Mean fitness") + ggtitle("Sex chromosomes", subtitle="Common sex: Males")

library(cowplot)
g1 <- ggplotGrob(auto.mal)
g2 <- ggplotGrob(sexchrom.fem)
g3 <- ggplotGrob(sexchrom.mal)


plot_grid(
  g1, g2, g3, ncol=3,
  labels = c('A)', 'B)', 'C)'),
  align="hv"
)

library(gridExtra)
grid.arrange(g1, g2, g3, ncol = 3)


