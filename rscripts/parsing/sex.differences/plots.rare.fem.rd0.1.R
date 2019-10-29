# Julio Rincones-Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com
# P.I. Dr. Heath Blackmon
# coleoguy@gmail.com
library(reshape)
library(viridis)
library(ggplot2)
# This script takes the fitness differences between males and females when females are the rare sex.


# Plot the recessive genetic architecture
rfem.xy.h0 <- read.csv("../../sex.bias/rscripts/parsing/sex.differences/h0.rfem.sexchrom.csv", 
                       row.names = 1, as.is = T, header = T, check.names = F)

# Force as matrix in order to generate melted results with Pop. size and OSR as separate variables
rfem.xy.h0 <- as.matrix(rfem.xy.h0)

melted_result <- melt(rfem.xy.h0[1:7,])

melted_result$Var2 <- as.factor(melted_result$Var2)

p <- ggplot(data = melted_result, 
            aes(x=Var1, y=Var2, fill=value)) +
  xlab("Operational sex ratio") +
  ylab("Common sex number (Males)") +
  geom_tile(color="white") +
  scale_fill_viridis(name = "Fitness difference", 
                     limit = c(-.35,.35), guide = "colourbar") +
  theme(legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7)) 
ggsave(filename = "mal.com.XY_h_0_s_0.5.pdf",
       plot=p, width=5, height=4, units="in")


# Plot the additive genetic architecture
rfem.xy.h0.5 <- read.csv("../sex.bias/rscripts/parsing/sex.differences/h0.5.rfem.sexchrom.csv", 
                       row.names = 1, as.is = T, header = T, check.names = F)

# Force as matrix in order to generate melted results with Pop. size and OSR as separate variables
rfem.xy.h0.5 <- as.matrix(rfem.xy.h0.5)

melted_result <- melt(rfem.xy.h0.5[1:7,])

melted_result$Var2 <- as.factor(melted_result$Var2)

p <- ggplot(data = melted_result, 
            aes(x=Var1, y=Var2, fill=value)) +
  xlab("Operational sex ratio") +
  ylab("Common sex number (Males)") +
  geom_tile(color="white") +
  scale_fill_viridis(name = "Fitness difference", 
                     limit = c(-.35,.35), guide = "colourbar") +
  theme(legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7)) 
ggsave(filename = "mal.com.XY_h_0.5_s_0.5.pdf",
       plot=p, width=5, height=4, units="in")

# Plot the dominant genetic architecture
rfem.xy.h1 <- read.csv("../sex.bias/rscripts/parsing/sex.differences/h1.rfem.sexchrom.csv", 
                         row.names = 1, as.is = T, header = T, check.names = F)

# Force as matrix in order to generate melted results with Pop. size and OSR as separate variables
rfem.xy.h1 <- as.matrix(rfem.xy.h1)

melted_result <- melt(rfem.xy.h1[1:7,])

melted_result$Var2 <- as.factor(melted_result$Var2)

p <- ggplot(data = melted_result, 
            aes(x=Var1, y=Var2, fill=value)) +
  xlab("Operational sex ratio") +
  ylab("Common sex number (Males)") +
  geom_tile(color="white") +
  scale_fill_viridis(name = "Fitness difference", 
                     limit = c(-.35,.35), guide = "colourbar") +
  theme(legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7)) 
ggsave(filename = "mal.com.XY_h_1_s_0.5.pdf",
       plot=p, width=5, height=4, units="in")
