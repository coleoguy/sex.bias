# Script for pulling the results from HD model and plotting difference 
# between mean male and female fitness
# J. Rincones-Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com

# Might be deprecated


# load functions
source("../functions/analysis.functions.HD.R")
load("../results/hd.rare.females.RData")
load("../results/hd.rare.males.RData")
res.rare.fem <- results
rm(results)

# set these values at first to determine the genetic architecture 
# and selection coefficient you wish
# to test.
s <- c(.1, .2, .5, .9)
h <- c(0, 0.5, 1)
osr <- c("OSR1","OSR.8","OSR.6",
         "OSR.4","OSR.2","OSR.1",
         "OSR.05")
comm.sex.num <- c("1000", "500", 
                  "100", "50")

# create a data structure to store results
results.fit <- as.data.frame(matrix(, 0, 8))
colnames(results.fit) <- c("mal.fit", "fem.fit", "mal-fem", "fem-mal", "comm.num",
                           "osr", "h", "s")


for(i in 1:length(comm.sex.num)){
  print(i)
  for(j in 1:length(osr)){
    for(k in 1:length(s)){
      curr.sims <- res.rare.fem[[i]][[j]][[k]]
      cur.results <- getFitness(data = curr.sims, h, s[k], comm.sex.num[i], osr[j])
      results.fit <- rbind(results.fit, cur.results)
    }
  }
}
    
  
foo1 <- data.frame(results.fit[,c(1,5:8)],
                rep("male", nrow(results.fit)))
foo2 <- data.frame(results.fit[,c(2,5:8)],
                   rep("female", nrow(results.fit)))
colnames(foo1)[c(1,6)] <- colnames(foo2)[c(1,6)] <- c("fit","sex")
xfoo <- rbind(foo1,foo2)
ggraptR(xfoo)
xfoos.9 <- xfoo[xfoo$s==0.9,]
ggraptR(xfoos.9)


library(viridis)
library(ggplot2)
library(reshape2)

melted_result <- melt(results.fit[1:7,])

melted_result$Var2 <- as.factor(melted_result$Var2)

p <- ggplot(data = melted_result, 
            aes(x=Var1, y=Var2, fill=value)) +
  xlab("Operational sex ratio") +
  ylab("Common sex number (Males)") +
  geom_tile(color="white") +
  scale_fill_viridis(name = "Fitness difference", 
                     limit = c(-.25,.25), guide = "colourbar") +
  theme(legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7)) 
ggsave(filename = "Mal.HD_h_.5_s_0.5.pdf",
       plot=p, width=5, height=4, units="in")

write.csv(results.fit, file= "Mal.HD.h.5.s0.5.csv")
