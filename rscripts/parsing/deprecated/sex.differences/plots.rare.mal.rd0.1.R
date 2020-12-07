# Julio Rincones-Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com
# P.I. Dr. Heath Blackmon
# coleoguy@gmail.com

library(viridis)
library(ggplot2)
library(reshape2)

# This script takes the fitness differences between males and females when males are the rare sex.

# Plot the recessive genetic architecture
rmal.xy.h0 <- read.csv("../sex.bias/rscripts/parsing/sex.differences/h0.rmal.sexchrom.csv", 
                       row.names = 1, as.is = T, header = T, check.names = F)

# Force as matrix in order to generate melted results with Pop. size and OSR as separate variables
rmal.xy.h0 <- as.matrix(rmal.xy.h0)

melted_result <- melt(rmal.xy.h0[1:7,])

melted_result$Var2 <- as.factor(melted_result$Var2)

p <- ggplot(data = melted_result, 
            aes(x=Var1, y=Var2, fill=value)) +
  xlab("Operational sex ratio") +
  ylab("Common sex number (Females)") +
  geom_tile(color="white") +
  scale_fill_viridis(name = "Fitness difference", 
                     limit = c(-.35,.35), guide = "colourbar") +
  theme(legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7)) 
ggsave(filename = "fem.com.XY_h_0_s_0.5.pdf",
       plot=p, width=5, height=4, units="in")


# Now plot the additive genetic architecture
rmal.xy.h0.5 <- read.csv("../sex.bias/rscripts/parsing/sex.differences/h0.5.rmal.sexchrom.csv", 
                       row.names = 1, as.is = T, header = T, check.names = F)

# Force as matrix in order to generate melted results with Pop. size and OSR as separate variables
rmal.xy.h0.5 <- as.matrix(rmal.xy.h0.5)

melted_result <- melt(rmal.xy.h0.5[1:7,])

melted_result$Var2 <- as.factor(melted_result$Var2)

p <- ggplot(data = melted_result, 
            aes(x=Var1, y=Var2, fill=value)) +
  xlab("Operational sex ratio") +
  ylab("Common sex number (Females)") +
  geom_tile(color="white") +
  scale_fill_viridis(name = "Fitness difference", 
                     limit = c(-.35,.35), guide = "colourbar") +
  theme(legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7)) 
ggsave(filename = "fem.com.XY_h_.5_s_0.5.pdf",
       plot=p, width=5, height=4, units="in")


# Now plot the dominant genetic architecture
rmal.xy.h1 <- read.csv("../sex.bias/rscripts/parsing/sex.differences/h1.rmal.sexchrom.csv", 
                         row.names = 1, as.is = T, header = T, check.names = F)

# Force as matrix in order to generate melted results with Pop. size and OSR as separate variables
rmal.xy.h1 <- as.matrix(rmal.xy.h1)

melted_result <- melt(rmal.xy.h1[1:7,])

melted_result$Var2 <- as.factor(melted_result$Var2)

p <- ggplot(data = melted_result, 
            aes(x=Var1, y=Var2, fill=value)) +
  xlab("Operational sex ratio") +
  ylab("Common sex number (Females)") +
  geom_tile(color="white") +
  scale_fill_viridis(name = "Fitness difference", 
                     limit = c(-.35,.35), guide = "colourbar") +
  theme(legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7)) 
ggsave(filename = "fem.com.XY_h_1_s_0.5.pdf",
       plot=p, width=5, height=4, units="in")

