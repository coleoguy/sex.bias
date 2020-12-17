
load("../results/hd.rare.males.RData")
results$
results$comm.sex1000$osr1
results$comm.sex1000$osr1$s0.1
results$comm.sex1000$osr1$s0.1$`h=0`

# Convert the list of results above as data frame and assign column names
result <- as.data.frame(matrix(,112000,5))
colnames(result) <- c("freq", "comm", "osr", "s", "h")

# What was the purpose of creating this vector and storing a value of one again?
rs <- 1

# nested for loops to iterate through the results

for(i in 1:length(results)){ # iterate through popsize
  print(i)
  for(j in 1:length(results$comm.sex1000)){ # iterate through osr
    for(k in 1:length(results$comm.sex1000$osr1)){ # iterate through s
      for(m in 1:length(results$comm.sex1000$osr1$s0.1)){ # iterate through h
        # Results of these loops saved in vector x 
        x <- results[[i]][[j]][[k]][[m]]
        # Need a little comment here. 
        re <- rs + length(x) -1
        # Store frequencies in column 1 of this data frame
        result[rs:re, 1] <- x
        # Stores the population size to column 2
        result[rs:re, 2] <- rep(names(results)[i], 1000)
        # Stores the number of the common sex to column 3
        result[rs:re, 3] <- rep(names(results$comm.sex1000)[j], 1000)
        # Stores the selection coefficient to column 4
        result[rs:re, 4] <- rep(names(results$comm.sex1000$osr1)[k], 1000)
        # Stores dominance type to column 5
        result[rs:re, 5] <- rep(names(results$comm.sex1000$osr1$s0.1)[m], 1000)
        rs <- rs + length(x)
      }
    }
  }
}

# subset by dominance
h0.5 <- result[result$h == "h=0.5",]
h0 <- result[result$h == "h=0", ]
h1 <- result[result$h == "h=1", ]


write.csv(h0.5, file="../../../../Haplodiploidy/parsed.HD.h05.rare.male.csv")
write.csv(h0, file="../../../../Haplodiploidy/parsed.HD.h0.rare.male.csv")
write.csv(h1, file="../../../../Haplodiploidy/parsed.HD.h1.rare.male.csv")
write.csv(result, file="parsed.HD.rare.male.csv")
library(ggraptR)
library(viridis)
ggraptR(result)

ggplot(result, aes(x= freq)) + 
  geom_histogram(aes(fill= as.factor(s)), stat="bin", position="dodge", alpha=0.5, bins=10) + 
  facet_grid(comm ~ osr) + 
  theme_grey() + theme(text=element_text(family="sans", face="plain", color="#000000", size=15, 
                                              hjust=0.5, vjust=0.5)) + 
  guides(fill=guide_legend(title="s")) + xlab("frequency of the allele benefitting the common sex (females)") + 
  ylab("simulation count") + scale_x_continuous(breaks=c(0,.5,1),labels=c(0,.5,1), name="frequency of allele benefitting the common sex (females)") +
  ggtitle("Fate of sexually antagonistic mutations when males are the rare sex") +
  scale_colour_viridis_d(aesthetics = "fill")


######################

result2 <- result
result2$comm <-  as.numeric(gsub("[a-z//.]",replacement = "" ,x= result2$comm))
result2$osr <- as.numeric(gsub("[a-z]",replacement = "" ,x= result2$osr))

ploth0 <- h0 
ploth0$comm <-  as.numeric(gsub("[a-z//.]",replacement = "" ,x= ploth0$comm))
ploth0$osr <- as.numeric(gsub("[a-z]",replacement = "" ,x= ploth0$osr))

ploth0.5 <- h0.5 
ploth0.5$comm <-  as.numeric(gsub("[a-z//.]",replacement = "" ,x= ploth0.5$comm))
ploth0.5$osr <- as.numeric(gsub("[a-z]",replacement = "" ,x= ploth0.5$osr))

ploth1 <- h1 
ploth1$comm <-  as.numeric(gsub("[a-z//.]",replacement = "" ,x= ploth1$comm))
ploth1$osr <- as.numeric(gsub("[a-z]",replacement = "" ,x= ploth1$osr))

ggplot(result2, aes(x = freq)) + 
  geom_histogram(aes(fill= as.factor(s)), stat="bin", position ="dodge", alpha=0.5, bins=10) + 
  facet_grid(comm ~ osr) + 
  theme_grey() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) +
  guides(fill=guide_legend(title="s")) + xlab("Frequency of the allele benefitting the common sex") + 
  ylab("count") +
  scale_x_continuous(breaks = c(0,.5,1), labels=c(0, .5, 1), name = "frequency of allele benefitting the common sex") +
  scale_colour_viridis_d(aesthetics = "fill")

ggplot(ploth0, aes(x = freq)) + 
  geom_histogram(aes(fill= as.factor(s)), stat="bin", position ="dodge", alpha=0.5, bins=10) + 
  facet_grid(comm ~ osr) + 
  theme_grey() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) +
  guides(fill=guide_legend(title="s")) + xlab("Frequency of the allele benefitting the common sex (females)") + 
  ylab("simulation count") +
  scale_x_continuous(breaks = c(0,.5,1), labels=c(0, .5, 1), name = "frequency of allele benefitting the common sex (females)") +
  ggtitle("Fate of sexually antagonistic mutations in a haplodiploidy scenario: h=0") +
  scale_colour_viridis_d(aesthetics = "fill")
  
ggplot(ploth0.5, aes(x = freq)) + 
    geom_histogram(aes(fill= as.factor(s)), stat="bin", position ="dodge", alpha=0.5, bins=10) + 
    facet_grid(comm ~ osr) + 
    theme_grey() + 
    theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) +
    guides(fill=guide_legend(title="s")) + xlab("Frequency of the allele benefitting the common sex (females)") + 
    ylab("simulation count") +
    scale_x_continuous(breaks = c(0,.5,1), labels=c(0, .5, 1), name = "frequency of allele benefitting the common sex (females)") +
    ggtitle("Fate of sexually antagonistic mutations in a haplodiploidy scenario:h=0.5") +
  scale_colour_viridis_d(aesthetics = "fill")
  
ggplot(ploth1, aes(x = freq)) + 
    geom_histogram(aes(fill= as.factor(s)), stat="bin", position ="dodge", alpha=0.5, bins=10) + 
    facet_grid(comm ~ osr) + 
    theme_grey() + 
    theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) +
    guides(fill=guide_legend(title="s")) + xlab("Frequency of the allele benefitting the common sex (females)") + 
    ylab("simulation count") +
    scale_x_continuous(breaks = c(0,.5,1), labels=c(0, .5, 1), name = "frequency of allele benefitting the common sex (females)") +
    ggtitle("Fate of sexually antagonistic mutations in a haplodiploidy scenario: h=1") +
  scale_colour_viridis_d(aesthetics = "fill")

h0.5 <- result[result$h == "h=0.5",]
h0 <- result[result$h == "h=0", ]
h1 <- result[result$h == "h=1", ]

ggraptR(h1)
ggraptR(h0)
ggraptR(h0.5)

#testing mean frequencies
dominant <- h1[h1$comm == "comm.sex100", ]
domosr <- dominant[dominant$osr == "osr0.05", ]
library(ggraptR)
ggraptR(domosr)

dominant2 <- h1[h1$comm == "comm.sex50", ]
domosr2 <- dominant2[dominant2$osr == "osr0.05", ]
ggraptR(domosr2)


recessive <- h0[h0$comm == "comm.sex100", ]
resosr <- recessive[recessive$osr == "osr0.05", ]
ggraptR(resosr)

recessive2 <- h0[h0$comm == "comm.sex50", ]
resosr2 <- recessive2[recessive2$osr == "osr0.05", ]
ggraptR(resosr2)

additive <- h0.5[h0.5$comm == "comm.sex100", ]
addosr <- additive[additive$osr == "osr0.05", ]
ggraptR(addosr)

additive2 <- h0.5[h0.5$comm == "comm.sex100", ]
addosr2 <- additive2[additive2$osr == "osr0.05", ]
ggraptR(addosr2)

