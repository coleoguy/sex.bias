
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
results.fit.mal <- matrix(NA, 7, 4)

# insert column and row names
colnames(results.fit.mal) <- c("50", "100", 
                               "500", "1000")
rownames(results.fit.mal) <- c("OSR1","OSR.8","OSR.6",
                               "OSR.4","OSR.2","OSR.1",
                               "OSR.05")

# create a data structure to store results
results.fit.fem <- matrix(NA, 7, 4)

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
      fit <-  getFitness(data = rare.female[[i]][[j-1]]$rd0.5$h0.5$s0.5, 
                         loc, h, s)
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
      fit <-  getFitness(data = rare.female[[i]][[j-1]]$rd0.5$h0.5$s0.5, 
                         loc, h, s)
      results.fit.fem[j, i] <- mean(fit$wbar.f)
    }
  }
}

# Now, let's reshape the data to make it long 

max.w <- rep(c(1), 28)
wbar <- cbind(melt(results.fit.fem[1:7,]), 
              melt(results.fit.mal[1:7,]), max.w)
colnames(wbar) <- c("OSR1", "Common.sex", "wbar.f", 
                        "OSR2", "Common.sex", "wbar.m", "max.w")
wbar <- wbar[(wbar$Common.sex == 100),]

# Don't need the OSR1, OSR2, and Common.sex columns.
wbar <- wbar[,-c(1:2, 4, 5)]

# Retain this matrix so that you can calculate rowMeans and then
# use it to introduce the breaks in the plot.
osr <- cbind(seq(from=.1,by=.5,length.out=7), 
             seq(from=.2,by=.5,length.out=7))

# create a separate OSR df to bind with the wbar df being created below.
OSR <- melt(osr)
wbar <- melt(wbar, id.vars = c("max.w"), 
             measure.vars = c("wbar.f", "wbar.m"))
wbar <- cbind(wbar, OSR$value)
colnames(wbar) <- c("max.w","sex","wbar", "OSR")

# store point size and alpha variables rather than hard code them.
point.size <- 2
alpha <- .7


auto.mal <- ggplot(wbar, aes(y=wbar, x=OSR)) + 
  geom_segment(aes(x=OSR, xend=OSR, y=wbar, yend=max.w),
               size=0.5, data=wbar, linetype="solid", 
               alpha= 0.5) + 
  geom_point(aes(colour=sex), stat="identity", position="identity", 
             alpha=alpha, size=point.size) + 
  geom_point(aes(x= OSR, y=max.w), data= wbar, size=point.size, 
           alpha = 1) +
  theme_light() + theme(legend.position='none',
                        legend.background = 
                          element_rect(fill="white",
                                       size=0.5, linetype="solid", 
                                       colour ="black")) + 
  ylim(0.7, 1) +
  scale_x_continuous(breaks= rowMeans(osr[,1:2]), 
                     labels= c("1",".8",".6",".4",".2",".1",".05")) +
  xlab("OSR") + 
  ylab("Mean fitness") + 
  ggtitle("Autosomes", subtitle = 'Common sex: Males')

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


# Now, let's reshape the data to make it long 

max.w <- rep(c(1), 28)
wbar <- cbind(melt(results.fit.fem[1:7,]), 
              melt(results.fit.mal[1:7,]), max.w)
colnames(wbar) <- c("OSR1", "Common.sex", "wbar.f", 
                    "OSR2", "Common.sex", "wbar.m", "max.w")
wbar <- wbar[(wbar$Common.sex == 100),]

# Don't need the OSR1, OSR2, and Common.sex columns.
wbar <- wbar[,-c(1:2, 4, 5)]

# Retain this matrix so that you can calculate rowMeans and then
# use it to introduce the breaks in the plot.
osr <- cbind(seq(from=.1,by=.5,length.out=7), 
             seq(from=.2,by=.5,length.out=7))

# create a separate OSR df to bind with the wbar df being created below.
OSR <- melt(osr)
wbar <- melt(wbar, id.vars = c("max.w"), 
             measure.vars = c("wbar.f", "wbar.m"))
wbar <- cbind(wbar, OSR$value)
colnames(wbar) <- c("max.w","sex","wbar", "OSR")

# store point size and alpha variables rather than hard code them.
point.size <- 2
alpha <- .7


sexchrom.fem <- ggplot(wbar, aes(y=wbar, x=OSR)) + 
  geom_segment(aes(x=OSR, xend=OSR, y=wbar, yend=max.w),
               size=0.5, data=wbar, linetype="solid", 
               alpha= 0.5) + 
  geom_point(aes(colour=sex), stat="identity", position="identity", 
             alpha=alpha, size=point.size) + 
  geom_point(aes(x= OSR, y=max.w), data= wbar, size=point.size, 
             alpha = 1) +
  theme_light() + theme(legend.position='none',
                        legend.background = 
                          element_rect(fill="white",
                                       size=0.5, linetype="solid", 
                                       colour ="black")) + 
  ylim(0.7, 1) +
  scale_x_continuous(breaks= rowMeans(osr[,1:2]), 
                     labels= c("1",".8",".6",".4",".2",".1",".05")) +
  xlab("OSR") + ylab("Mean fitness") + 
  ggtitle("Sex chromosomes", subtitle="Common sex: Females")


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

# Now, let's reshape the data to make it long 

max.w <- rep(c(1), 28)
wbar <- cbind(melt(results.fit.fem[1:7,]), 
              melt(results.fit.mal[1:7,]), max.w)
colnames(wbar) <- c("OSR1", "Common.sex", "wbar.f", 
                    "OSR2", "Common.sex", "wbar.m", "max.w")
wbar <- wbar[(wbar$Common.sex == 100),]

# Don't need the OSR1, OSR2, and Common.sex columns.
wbar <- wbar[,-c(1:2, 4, 5)]

# Retain this matrix so that you can calculate rowMeans and then
# use it to introduce the breaks in the plot.
osr <- cbind(seq(from=.1,by=.5,length.out=7), 
             seq(from=.2,by=.5,length.out=7))

# create a separate OSR df to bind with the wbar df being created below.
OSR <- melt(osr)
wbar <- melt(wbar, id.vars = c("max.w"), 
             measure.vars = c("wbar.f", "wbar.m"))
wbar <- cbind(wbar, OSR$value)
colnames(wbar) <- c("max.w","sex","wbar", "OSR")

# store point size and alpha variables rather than hard code them.
point.size <- 2
alpha <- .7


sexchrom.mal <- ggplot(wbar, aes(y=wbar, x=OSR)) + 
  geom_segment(aes(x=OSR, xend=OSR, y=wbar, yend=max.w),
               size=0.5, data=wbar, linetype="solid", 
               alpha= 0.5) + 
  geom_point(aes(colour=sex), stat="identity", position="identity", 
             alpha=alpha, size=point.size) + 
  geom_point(aes(x= OSR, y=max.w), data= wbar, size=point.size, 
             alpha = 1) +
  theme_light() + theme(legend.position= 'none',
                        legend.background = 
                          element_rect(fill="white",
                                       size=0.5, linetype="solid", 
                                       colour ="black")) + 
  ylim(0.7, 1) +
  scale_x_continuous(breaks= rowMeans(osr[,1:2]), 
                     labels= c("1",".8",".6",".4",".2",".1",".05")) +
  xlab("OSR") + ylab("Mean fitness") + 
  ggtitle("Sex chromosomes", subtitle="Common sex: Males")


library(cowplot)
g1 <- ggplotGrob(auto.mal)
g2 <- ggplotGrob(sexchrom.fem)
g3 <- ggplotGrob(sexchrom.mal)


plot_grid(
  g1, g2, g3, ncol=3,
  labels = c('A)', 'B)', 'C)'),
  align="hv"
)

# grab the legend from the dummy plot below 
# so it can be sized and positioned in adobe acrobat.

dummy <- ggplot(wbar, aes(y=wbar, x=OSR)) + 
  geom_segment(aes(x=OSR, xend=OSR, y=wbar, yend=max.w),
               size=1, data=wbar, linetype="solid", 
               alpha= alpha) + 
  geom_point(aes(colour=sex), stat="identity", position="identity", 
             alpha=alpha, size=point.size) + 
  geom_point(aes(x= OSR, y=max.w), data= wbar, size=point.size, 
             alpha = 1) +
  theme_light() + theme(legend.position= 'right',
                        legend.background = 
                          element_rect(fill="white",
                                       size=0.5, linetype="solid", 
                                       colour ="black")) + 
  ylim(0.7, 1) +
  scale_x_continuous(breaks= rowMeans(osr[,1:2]), 
                     labels= c("1",".8",".6",".4",".2",".1",".05")) +
  xlab("OSR") + ylab("Mean fitness") + 
  ggtitle("Sex chromosomes", subtitle="Common sex: Males")

legend <- get_legend(dummy)

library(gridExtra)
# plot and save the legend
grid.arrange(legend)

########
