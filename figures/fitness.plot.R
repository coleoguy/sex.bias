# Fitness plots


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
max.fit <- rep(c(1), 28)
wbar.h05 <- cbind(wbar.fem.h05, wbar.mal.h05$value, max.fit)

colnames(wbar.h05) <- c("OSR", "Common.sex", "wbar.f", "wbar.m", "maxw")


wbar.plot <- wbar.h05[(wbar.h05$Common.sex == 1000),]


ggplot(wbar.plot, aes(x=OSR, y=maxw)) + 
  geom_segment(aes(x=OSR, xend=OSR, y=wbar.m, yend=maxw), 
               size=1, data=wbar.plot, colour="black", linetype="solid", alpha = 0.5) + 
  geom_segment(aes(x=OSR, xend=OSR, y=wbar.f, yend=maxw), 
               size=1, data=wbar.plot, colour="blue", linetype="dashed", alpha = 0.5) + 
  geom_point(aes(x= OSR, y=wbar.m), data= wbar.plot, size=5, colour= "black", alpha = 0.5) +
  geom_point(aes(x= OSR, y=maxw), data= wbar.plot, size=5, colour= "black", alpha = 0.5) +
  geom_point(aes(x= OSR, y=wbar.f), data= wbar.plot, size=5, colour= "blue", alpha = 0.5) +
  geom_point(aes(x= OSR, y=maxw), data= wbar.plot, size=5, colour= "blue", alpha = 0.5) +
  theme_light() +
  guides(fill=guide_legend(title="Common.sex")) + xlab("OSR") + ylab("Mean fitness")

ggraptR(wbar.plot)

ggplot(wbar.fem, aes(x=OSR, y=wbar)) + 
  geom_segment(aes(x=OSR, xend=OSR, y=wbar, yend=maxw), 
               size=1, data=wbar.fem, colour="black", linetype="solid") + 
  geom_point(aes(x= OSR, y=wbar), data= wbar.fem, size=5, colour= "black") +
  geom_point(aes(x= OSR, y=maxw), data= wbar.fem, size=5, colour= "black") +
  guides(fill=guide_legend(title="Common.sex")) + xlab("OSR") + ylab("Mean fitness")

library(viridis)
ggplot(wbar.fem.h05, aes(x = as.factor(Common.sex), y = wbar)) +
  geom_point(aes(x = OSR, y = wbar, colour = factor(Common.sex)),
    data = wbar.fem.h05, size = 6.5, alpha = .5) +
  geom_point(aes(x = OSR, y = maxw), data = wbar.fem.h05, size = 6.5,
    alpha = .2, colour = "black") +
  guides(fill = guide_legend(title = "Common.sex")) + xlab("OSR") + ylab("Mean fitness") +
  theme_light() +
  ggtitle("Mean vs Maximal Fitness - Common sex males") +
  scale_color_viridis(discrete=TRUE, option = "D")


ggplot(wbar.fem.h05, aes(x = as.factor(Common.sex), y = wbar)) +
  geom_point(aes(x = OSR, y = wbar, colour = factor(Common.sex)),
             data = wbar.fem.h05, size = 6.5, alpha = .5) +
  geom_point(aes(x = OSR, y = maxw), data = wbar.fem.h05, size = 6.5,
             alpha = .2, colour = "black") +
  
  guides(fill = guide_legend(title = "Common.sex")) + xlab("OSR") + ylab("Mean fitness") +
  theme_light() +
  ggtitle("Mean vs Maximal Fitness - Common sex males") +
  scale_color_viridis(discrete=TRUE, option = "D")


ggraptR(wbar.fem.h05)
  

# load data to be plotted
h0.auto <- read.csv("XY.model/csv.rare.mal/h0.autosomes.csv", 
                    header = TRUE, sep = ",", as.is = T, check.names = F,
                    row.names = 1)

# make data long
melted_result <- melt(h0.auto[1:7,])

osr <- rep(c("OSR1", "OSR.8", "OSR.6", "OSR.4", "OSR.2", "OSR.1", "OSR.05"), 4)
OSR <- rep(c(1, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05), 4)
tidy.h0 <- cbind(melted_result, OSR, osr)
p1000 <- mean(res.base[["pop1000"]][["rd0.1"]][["h0.5"]][["s0.5"]][,3])


wm.a11 <- 1 + s
wm.a12 <- 1 + h * s
wm.a22 <- 1
wf.a11 <- 1 / (1 + s)
wf.a12 <- 1 / (1 + h * s)
wf.a22 <- 1
                  