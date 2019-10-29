# Julio Rincones-Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com 
# P.I. Dr. Heath Blackmon
# coleoguy@gmail.com

# Autosomal loci plots for rare females

library(viridis)
library(ggplot2)
library(reshape2)

# load data to be plotted


# plot recessive 

melted_result <- melt(h0.autosomes[1:7,])

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
ggsave(filename = "mal.com.auto_h_0_s_0.5.pdf",
       plot=p, width=5, height=4, units="in")


# plot additive 

melted_result <- melt(h.5.autosomes[1:7,])

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
ggsave(filename = "mal.com.auto_h_.5_s_0.5.pdf",
       plot=p, width=5, height=4, units="in")

# plot dominant

melted_result <- melt(h1.autosomes[1:7,])

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
ggsave(filename = "mal.com.auto_h_1_s_0.5.pdf",
       plot=p, width=5, height=4, units="in")





