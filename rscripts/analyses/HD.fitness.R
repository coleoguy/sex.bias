# Julio Rincones Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com
# P.I. Dr. Heath Blackmon
# coleoguy@gmail.com

# Script to calculate fitness differences between males and females under our haplodiploidy model.
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggraptR)
library(viridis)

# load functions
source("../functions/fitnessHD.R")
load("../results/hd.rare.females.RData")
res.rare.fem <- results
load("../results/hd.rare.males.RData")
res.rare.mal <- results
rm(results)


# We will test first when males are the common sex.
# first we test when the dominance factor = 0 (i.e. recessive)

h <- 0
s <- .5
dom <- "h=0"

# create a data structure to store results
results.hd.0 <- matrix(, 7, 4)

colnames(results.hd.0) <- c("1000", "500", "100", "50")

rownames(results.hd.0) <- c("1", "0.8", "0.6",
                            "0.4", "0.2", "0.1",
                            "0.5")

# results structure res.rare.fem[["males1000"]][["osr1"]][["s0.1"]]

for (i in 1:4) {
  for (j in 1:7) {
    for (k in 1:4) {
      if (k == 3) {
        results <- getFitness(data = res.rare.fem[[i]][[j]][[k]], dom, h, s)
        results.hd.0[j, i] <- getDifference.mal(results)
      }
    }
  }
}

results.hd.0 <- as.data.frame(results.hd.0)
hd0 <-
  gather(results.hd.0,
         "1000", "500", "100", "50",
         key = "Common sex number",
         value = "Fitness difference")

osr <- rep(c(1, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05), 4)
h <- rep(0, 28)
s <- rep(0.5, 28)
hd0 <- cbind(hd0, osr, h, s)


# Then when the dominance factor = 0.5 (i.e. additive)
h <- 0.5
s <- .5
dom <- "h=0.5"

results.hd.5 <- matrix(, 7, 4)

colnames(results.hd.5) <- c("1000", "500", "100", "50")

rownames(results.hd.5) <- c("1", "0.8", "0.6",
                            "0.4", "0.2", "0.1",
                            "0.5")

for (i in 1:4) {
  for (j in 1:7) {
    for (k in 1:4) {
      if (k == 3) {
        results <- getFitness(data = res.rare.fem[[i]][[j]][[k]], dom, h, s)
        results.hd.5[j, i] <- getDifference.mal(results)
      }
    }
  }
}

results.hd.5 <- as.data.frame(results.hd.5)
hd.5 <-
  gather(results.hd.5,
         "1000", "500", "100", "50",
         key = "Common sex number",
         value = "Fitness difference")

osr <- rep(c(1, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05), 4)
h <- rep(0.5, 28)
s <- rep(0.5, 28)
hd.5 <- cbind(hd.5, osr, h, s)

# Then when the dominance factor = 1 (i.e. dominant)

h <-  1
s <- .5
dom <- "h=1"

results.hd.1 <- matrix(, 7, 4)

colnames(results.hd.1) <- c("1000", "500", "100", "50")

rownames(results.hd.1) <- c("1", "0.8", "0.6",
                            "0.4", "0.2", "0.1",
                            "0.5")

for (i in 1:4) {
  for (j in 1:7) {
    for (k in 1:4) {
      if (k == 3) {
        results <- getFitness(data = res.rare.fem[[i]][[j]][[k]], dom, h, s)
        results.hd.1[j, i] <- getDifference.mal(results)
      }
    }
  }
}

results.hd.1 <- as.data.frame(results.hd.1)
hd.1 <-
  gather(results.hd.1,
         "1000", "500", "100", "50",
         key = "Common sex number",
         value = "Fitness difference")

osr <- rep(c(1, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05), 4)
h <- rep(1, 28)
s <- rep(0.5, 28)
hd.1 <- cbind(hd.1, osr, h, s)

# Bind the results, transform to tidy format and save.
hd.male.comm <- rbind(hd0, hd.5, hd.1)
hd.male.comm <- as.data.frame(hd.male.comm)


resultsplot <- hd.male.comm[(hd.male.comm$`Common sex number` == 1000),]
resultsplot <- hd.male.comm[(hd.male.comm$`Common sex number` == 500),]
resultsplot <- hd.male.comm[(hd.male.comm$`Common sex number` == 100),]
resultsplot <- hd.male.comm[(hd.male.comm$`Common sex number` == 50),]


ggplot(resultsplot, aes(y = `Fitness difference`, x = osr)) + geom_point(
  aes(colour = as.factor(h)),
  stat = "identity",
  position = "identity",
  alpha = 0.75,
  size = 3
) + theme_light() + theme(text = element_text(
  family = "sans",
  face = "plain",
  color = "#000000",
  size = 15,
  hjust = 0.5,
  vjust = 0.5
)) + scale_size(range = c(1, 3)) + xlab("osr") + ylab("Fitness.difference") + scale_x_reverse() +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) 

##
# Now let's test when females are the common sex

# first we test when the dominance factor = 0 (i.e. recessive)

h <- 0
s <- .5
dom <- "h=0"

# create a data structure to store results
results.hd.0 <- matrix(, 7, 4)

colnames(results.hd.0) <- c("1000", "500", "100", "50")

rownames(results.hd.0) <- c("1", "0.8", "0.6",
                            "0.4", "0.2", "0.1",
                            "0.5")


for (i in 1:4) {
  for (j in 1:7) {
    for (k in 1:4) {
      if (k == 3) {
        results <- getFitness(data = res.rare.mal[[i]][[j]][[k]], dom, h, s)
        results.hd.0[j, i] <- getDifference.fem(results)
      }
    }
  }
}

results.hd.0 <- as.data.frame(results.hd.0)
hd0 <-
  gather(results.hd.0,
         "1000", "500", "100", "50",
         key = "Common sex number",
         value = "Fitness difference")

osr <- rep(c(1, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05), 4)
h <- rep(0, 28)
s <- rep(0.5, 28)
hd0 <- cbind(hd0, osr, h, s)


# Then when the dominance factor = 0.5 (i.e. additive)
h <- 0.5
s <- .5
dom <- "h=0.5"

results.hd.5 <- matrix(, 7, 4)

colnames(results.hd.5) <- c("1000", "500", "100", "50")

rownames(results.hd.5) <- c("1", "0.8", "0.6",
                            "0.4", "0.2", "0.1",
                            "0.5")

for (i in 1:4) {
  for (j in 1:7) {
    for (k in 1:4) {
      if (k == 3) {
        results <- getFitness(data = res.rare.mal[[i]][[j]][[k]], dom, h, s)
        results.hd.5[j, i] <- getDifference.fem(results)
      }
    }
  }
}

results.hd.5 <- as.data.frame(results.hd.5)
hd.5 <-
  gather(results.hd.5,
         "1000", "500", "100", "50",
         key = "Common sex number",
         value = "Fitness difference")

osr <- rep(c(1, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05), 4)
h <- rep(0.5, 28)
s <- rep(0.5, 28)
hd.5 <- cbind(hd.5, osr, h, s)

# Then when the dominance factor = 1 (i.e. dominant)

h <-  1
s <- .5
dom <- "h=1"

results.hd.1 <- matrix(, 7, 4)

colnames(results.hd.1) <- c("1000", "500", "100", "50")

rownames(results.hd.1) <- c("1", "0.8", "0.6",
                            "0.4", "0.2", "0.1",
                            "0.5")

for (i in 1:4) {
  for (j in 1:7) {
    for (k in 1:4) {
      if (k == 3) {
        results <- getFitness(data = res.rare.mal[[i]][[j]][[k]], dom, h, s)
        results.hd.1[j, i] <- getDifference.fem(results)
      }
    }
  }
}

results.hd.1 <- as.data.frame(results.hd.1)
hd.1 <-
  gather(results.hd.1,
         "1000", "500", "100", "50",
         key = "Common sex number",
         value = "Fitness difference")

osr <- rep(c(1, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05), 4)
h <- rep(1, 28)
s <- rep(0.5, 28)
hd.1 <- cbind(hd.1, osr, h, s)

# Bind the results, transform to tidy format and save.
hd.female.comm <- rbind(hd0, hd.5, hd.1)
hd.female.comm <- as.data.frame(hd.female.comm)

resultsplot <- hd.female.comm[(hd.female.comm$`Common sex number` == 1000),]
resultsplot <- hd.female.comm[(hd.female.comm$`Common sex number` == 500),]
resultsplot <- hd.female.comm[(hd.female.comm$`Common sex number` == 100),]
resultsplot <- hd.female.comm[(hd.female.comm$`Common sex number` == 50),]


ggplot(resultsplot, aes(y = `Fitness difference`, x = osr)) + geom_point(
  aes(colour = as.factor(h)),
  stat = "identity",
  position = "identity",
  alpha = 0.5,
  size = 3
) + theme_light() + theme(text = element_text(
  family = "sans",
  face = "plain",
  color = "#000000",
  size = 15,
  hjust = 0.5,
  vjust = 0.5
)) + scale_size(range = c(1, 3)) + xlab("osr") + ylab("Fitness.difference") + scale_x_reverse() + 
  scale_fill_viridis(discrete = TRUE) + scale_color_viridis(discrete = TRUE)
