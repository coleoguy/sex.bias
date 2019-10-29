# Script for pulling the results from HD model and plotting difference 
# between mean male and female fitness
# J. Rincones-Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com

# Might be deprecated


# load functions
source("../functions/analysis.functions.HD.R")
load("../results/hd.rare.males.RData")
res.rare.mal <- results
rm(results)

#this formal also has to be changed manually. The redundancy is due to the fact that this value of
# h feeds into the fitness functions used to calculate genotype-specific fitnesses. While h2 is used to 
# pull the relevant allele frequency from the results list from our simulations.
h <- .5 
# We have simulated different selection coefficients, you can choose from 0.1, 0.2, 0.5, and 0.9.
s <- .5 
#this  formal has to be changed manually for the function to work. The genetic architectures
# are "h=0", "h=0.5", and "h=1"
h2 <- "h=0.5" 

# create a data structure to store results
results.fit <- matrix(, 7, 4)


# insert column and row names. Verify the data structure of your source of results as it is otherwise
# easy to be careless and fill this table incorrectly.
colnames(results.fit) <- c("1000", "500", 
                           "100", "50")
rownames(results.fit) <- c("OSR1","OSR.8","OSR.6",
                           "OSR.4","OSR.2","OSR.1",
                           "OSR.05")


for(i in 1:7){
  for(j in 1:4){
    results <- getFitness(data = res.rare.mal[[j]][[i]]$s0.5, h2, h, s)
    results.fit[i, j] <- getDifference.fem(results) 
  }

}


library(viridis)
library(ggplot2)
library(reshape2)

melted_result <- melt(results.fit[1:7,])

melted_result$Var2 <- as.factor(melted_result$Var2)

p <- ggplot(data = melted_result, 
            aes(x=Var1, y=Var2, fill=value)) +
  xlab("Operational sex ratio") +
  ylab("Common sex number (Female)") +
  geom_tile(color="white") +
  scale_fill_viridis(name = "Fitness difference", 
                     limit = c(-.2,.2), guide = "colourbar") +
  theme(legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7)) 
ggsave(filename = "Fem.HD_h_.5_s_0.5.pdf",
       plot=p, width=5, height=4, units="in")

# change the file name accordingly if you are creating a table with different genetic architecture and 
# selection coefficient.
write.csv(results.fit, file= "Fem.HD.h.5.s0.5.csv")


