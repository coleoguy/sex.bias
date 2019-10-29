# this script analyzes results of the OSR-SA project
# utilizing functions in the analysis.functions.R 
# file.
# H. Blackmon
# coleoguy@gmail.com
# 5 July 2019
# Additional code J. Rincones-Gamboa
# jgamboa@bio.tamu.edu

# load functions
source("../sex.bias/rscripts/functions/analysis.functions.R")
load("../sex.bias/rscripts/results/base.comp.model.RData")
res.base <- results
load("../sex.bias/rscripts/results/rare.female.model.RData")
res.rare.fem <- results
load("../sex.bias/rscripts/results/rare.male.model.RData")
res.rare.mal <- results
rm(results)


# create a data structure to store results
results.fit <- matrix(, 7, 4)

# insert column and row names
colnames(results.fit) <- c("50", "100", 
                           "500", "1000")
rownames(results.fit) <- c("OSR1","OSR.8","OSR.6",
                           "OSR.4","OSR.2","OSR.1",
                           "OSR.05")

# structure of results
# res.rare.fem$males50$females0.8$rd0$h0$s0.1
# res.base$pop50$rd0$h0$s0.1
loc <- "auto"
h <- 1
s <- .5


# THIS CODE EVALUATES AUTOSOMAL LOCI
# cycle through population sizes
for(i in 1:4){
  # cycle through OSR values
  for(j in 1:7){
    if(j == 1){
      # looks at base model
      fit <- getFitness(data = res.base[[i]]$rd0.5$h0$s0.5, loc, h, s)
      results.fit[j, i] <- getDifference.mal(fit) 
    }
    if(j > 1){
      # looks at bias model
      fit <- getFitness(data = res.rare.fem[[i]][[j-1]]$rd0.5$h0$s0.5, loc, h, s)
      results.fit[j, i] <- getDifference.mal(fit) 
    }
  }
}
## Heath stopped editing here

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
  scale_fill_viridis(name = "Net Masculinization", 
                     limit = c(-.35,.35), guide = "colourbar") +
  theme(legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7)) 
ggsave(filename = "mal.com.auto_h_1_s_0.5.pdf",
       plot=p, width=5, height=4, units="in")

write.csv(results.fit, file= "mal.com.auto.h1.s05.csv")

# Now set the location to a sex chromosome with respective h and s values
loc <- "sex"
h <- 0.5
s <- .5

# THIS CODE EVALUATES X and Y LOCI
# cycle through population sizes
for(i in 1:4){
  # cycle through OSR values
  for(j in 1:7){
    if(j == 1){
      # looks at base model
      fit <- getFitness(data = res.base[[i]]$rd0.1$h0.5$s0.5, loc, h, s)
      results.fit[j, i] <- getDifference.mal(fit) 
    }
    if(j > 1){
      # looks at bias model
      fit <- getFitness(data = res.rare.fem[[i]][[j-1]]$rd0.1$h0.5$s0.5, loc, h, s)
      results.fit[j, i] <- getDifference.mal(fit) 
    }
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
  ylab("Common sex number (Males)") +
  geom_tile(color="white") +
  scale_fill_viridis(name = "Net Masculinization", 
                     limit = c(-.35,.35), guide = "colourbar") +
  theme(legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7)) 
ggsave(filename = "mal.com.XY_h_0.5_s_0.5.pdf",
       plot=p, width=5, height=4, units="in")

write.csv(results.fit, file= "mal.com.XY.h0.5.s05.csv")
