# Julio Rincones Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com
# P.I. Dr. Heath Blackmon
# coleoguy@gmail.com

# Script to plot figure 2 for the manuscript
library(ggplot2)
library(viridis)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(ggpubr)
library(grid)
library(lattice)

# Set your working directory to 'figures' folder

# We will do the first elements by plotting them for the XY model (dominance and selection on the first 
# two plots)


resultsplot <- read.csv("../../sex.bias/rscripts/parsing/Autosome.rare.male.csv", row.names = 1, 
                        as.is = T, header = T, check.names = F)
resultsplot <- resultsplot[(resultsplot$rd == 0.5),]
resultsplot <- resultsplot[(resultsplot$h != "h99"),]
resultsplot <- resultsplot[(resultsplot$s == "s0.5"),]
resultsplot <- resultsplot[(resultsplot$comm == 500),]

resultsplot2 <- resultsplot
resultsplot2$s <- c(0.5)
colnames(resultsplot2) <- c("frequency", "common.sex", "osr", "rd", "h", "s")


# Facet by dominance to show the effect of different genetic architectures. Comm.sex = 500. OSR = 1
resultsplot3 <- resultsplot2[resultsplot2$osr == 1,]

dominance <- ggplot(resultsplot3, aes(y=frequency, x=as.factor(h))) +
  ylim(0, 1) +
  geom_hline(yintercept = .5, alpha = .5) +
  geom_violin(aes(fill=as.factor(h)), stat="ydensity", position="dodge", 
              alpha=0.8, trim=TRUE, scale="area") +
  theme_light() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none") +
  ylab("Allele frequency") +
  ggtitle("Effect of genetic architecture") +
  scale_fill_viridis(discrete=TRUE, option = "B") +
  scale_color_viridis(discrete=TRUE, option = "B")

# Facet by recombination distance but only plot h=0.5. Comm.sex = 500, OSR = 1

resultsplot <- read.csv("../../sex.bias/rscripts/parsing/Autosome.rare.male.csv", row.names = 1, 
                        as.is = T, header = T, check.names = F)
resultsplot <- resultsplot[(resultsplot$h != "h99"),]
resultsplot <- resultsplot[(resultsplot$rd != 0),]
resultsplot <- resultsplot[(resultsplot$s == "s0.5"),]
resultsplot <- resultsplot[(resultsplot$h == "h0.5"),]
resultsplot <- resultsplot[(resultsplot$comm == 500),]
resultsplot <- resultsplot[(resultsplot$osr == 1),]

resultsplot2 <- resultsplot
resultsplot2$s <- c(0.5)
colnames(resultsplot2) <- c("frequency", "common.sex", "osr", "rd", "h", "s")

rd <- ggplot(resultsplot2, aes(y=1-frequency, x=as.factor(rd))) +
  ylim(0, 1) +
  geom_hline(yintercept = .5, alpha = .5) +
  geom_violin(aes(fill=as.factor(rd)), stat="ydensity", position="dodge", 
              alpha=0.8, trim=TRUE, scale="area") +
  theme_light() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none") +
  ylab("Allele frequency") +
  ggtitle("Effect of recombination") +
  scale_fill_viridis(discrete=TRUE, option = "B") +
  scale_color_viridis(discrete=TRUE, option = "B")


# Facet by selection coefficients but only plot h=0.5. Comm.sex = 100. OSR = 0.2

resultsplot <- read.csv("../../sex.bias/rscripts/parsing/Autosome.rare.male.csv", row.names = 1, 
                        as.is = T, header = T, check.names = F)
resultsplot <- resultsplot[(resultsplot$rd == 0.5),]
resultsplot <- resultsplot[(resultsplot$h == "h0.5"),]
resultsplot <- resultsplot[(resultsplot$comm == 500),]
colnames(resultsplot) <- c("frequency", "common.sex", "osr", "rd", "h", "s")

resultsplot4 <- resultsplot[(resultsplot$osr == 0.1),]


selection <- ggplot(resultsplot4, aes(y=1-frequency, x=as.factor(common.sex))) + 
  ylim(0, 1) +
  geom_hline(yintercept = .5, alpha = .5) +
  geom_violin(aes(fill=as.factor(s)), stat="ydensity", position="dodge", 
              alpha=0.8, trim=TRUE, scale="area") + 
  theme_light() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none") +
  ylab("Allele frequency") +
  ggtitle("Effect of selection") +
  scale_fill_viridis(discrete=TRUE, option = "B") +
  scale_color_viridis(discrete=TRUE, option = "B")
print(selection)



# Autosomal. Common sex (females) = 500, s=0.5, h=0.5, rd=0.5

resultsplot <- read.csv("../../sex.bias/rscripts/parsing/Autosome.rare.male.csv", row.names = 1, 
                        as.is = T, header = T, check.names = F)
resultsplot <- resultsplot[(resultsplot$rd == 0.5),]
resultsplot <- resultsplot[(resultsplot$h == "h0.5"),]
resultsplot <- resultsplot[(resultsplot$s == "s0.5"),]
resultsplot <- resultsplot[(resultsplot$comm == 500),]

resultsplot2 <- resultsplot
resultsplot2$s <- c(0.5)
colnames(resultsplot2) <- c("frequency", "common.sex", "osr", "rd", "h", "s")


autosomal <- ggplot(resultsplot2, aes(y=1-frequency, x=as.factor(common.sex))) + 
  ylim(0, 1) +
  geom_hline(yintercept = .5, alpha = .5) +
  geom_violin(aes(fill=as.factor(osr)), stat="ydensity", position="dodge", 
              alpha=0.8, trim=TRUE, scale="area") + 
  theme_light() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none") + 
  ylab("Allele frequency") +
  ggtitle("Autosomal") +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE)


# X chromosome. Common sex (females) = 500, s=0.5, h=0.5, rd=0.1

resultsplot <- read.csv("../../sex.bias/rscripts/parsing/Xchrom.rare.male.csv", row.names = 1, 
                        as.is = T, header = T, check.names = F)
resultsplot <- resultsplot[(resultsplot$rd == 0.1),]
resultsplot <- resultsplot[(resultsplot$h == "h0.5"),]
resultsplot <- resultsplot[(resultsplot$s == "s0.5"),]
resultsplot <- resultsplot[(resultsplot$comm == 500),]

resultsplot2 <- resultsplot
resultsplot2$s <- c(0.5)
colnames(resultsplot2) <- c("frequency", "common.sex", "osr", "rd", "h", "s")

xchrom <- ggplot(resultsplot2, aes(y=frequency, x=as.factor(common.sex))) +
  ylim(0, 1) +
  geom_hline(yintercept = .5, alpha = .5) +
  geom_violin(aes(fill=as.factor(osr)), stat="ydensity", position="dodge", 
              alpha=0.8, trim=TRUE, scale="area") + 
  theme_light() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none") +
  ylab("Allele frequency") +
  ggtitle("X chromosome") +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE)

# Y chromosome. Common sex (males) = 500, s=0.5, h=0.5, rd=0.1

resultsplot <- read.csv("../../sex.bias/rscripts/parsing/Ychrom.rare.female.csv", row.names = 1, 
                        as.is = T, header = T, check.names = F)
resultsplot <- resultsplot[(resultsplot$rd == 0.1),]
resultsplot <- resultsplot[(resultsplot$h == "h0.5"),]
resultsplot <- resultsplot[(resultsplot$s == "s0.5"),]
resultsplot <- resultsplot[(resultsplot$comm == 500),]

resultsplot2 <- resultsplot
resultsplot2$s <- c(0.5)
colnames(resultsplot2) <- c("frequency", "common.sex", "osr", "rd", "h", "s")

ychrom <- ggplot(resultsplot2, aes(y=frequency, x=as.factor(common.sex))) + 
  ylim(0, 1) +
  geom_hline(yintercept = .5, alpha = .5) +
  geom_violin(aes(fill=as.factor(osr)), stat="ydensity", position="dodge", 
              alpha=0.8, trim=TRUE, scale="area") + 
  theme_light() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none") +
  ylab("Allele frequency") +
  ggtitle("Y chromosome") +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE)

# Haplodiploidy rare males (Common sex = 500, s=0.5, h=0.5)

hd.rare.mal <- read.csv("../../sex.bias/figures/HD/parsed.HD.rare.male.csv", 
                        header = TRUE, sep = ",", as.is = T, check.names = F,
                        row.names = 1)

result.s <- hd.rare.mal[(hd.rare.mal$s == "s0.5"),]
result.h <- result.s[(result.s$h == "h=0.5"),]
result.h$h <- c(0.5)
result.h$s <- c(0.5)
results.500 <- result.h[(result.h$comm == "comm.sex500"),]
results.500$comm <- rep(500, length(results.500$comm))
results.500$osr[results.500$osr == "osr1"] <- 1
results.500$osr[results.500$osr == "osr0.8"] <- 0.8
results.500$osr[results.500$osr == "osr0.6"] <- 0.6
results.500$osr[results.500$osr == "osr0.4"] <- 0.4
results.500$osr[results.500$osr == "osr0.2"] <- 0.2
results.500$osr[results.500$osr == "osr0.1"] <- 0.1
results.500$osr[results.500$osr == "osr0.05"] <- 0.05


hd.fem <- ggplot(results.500, aes(y=freq, x=as.factor(comm))) +
  ylim(0, 1) +
  geom_hline(yintercept = .5, alpha = .5) +
  geom_violin(aes(fill=as.factor(osr)), stat="ydensity", position="dodge", alpha=0.8, trim=TRUE, 
              scale="area") + theme_light() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none") +
  ylab("Allele frequency") +
  ggtitle("HD model - Rare males") +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) 


# Haplodiploidy rare females (Common sex = 500, s=0.5, h=0.5)
hd.rare.fem <- read.csv("../../sex.bias/figures/HD/parsed.HD.rare.female.csv", 
                        header = TRUE, sep = ",", as.is = T, check.names = F,
                        row.names = 1)

result.s <- hd.rare.fem[(hd.rare.fem$s == "s0.5"),]
result.h <- result.s[(result.s$h == "h=0.5"),]
result.h$h <- c(0.5)
result.h$s <- c(0.5)
results.500 <- result.h[(result.h$comm == "comm.sex500"),]
results.500$comm <- rep(500, length(results.500$comm))
results.500$osr[results.500$osr == "osr1"] <- 1
results.500$osr[results.500$osr == "osr0.8"] <- 0.8
results.500$osr[results.500$osr == "osr0.6"] <- 0.6
results.500$osr[results.500$osr == "osr0.4"] <- 0.4
results.500$osr[results.500$osr == "osr0.2"] <- 0.2
results.500$osr[results.500$osr == "osr0.1"] <- 0.1
results.500$osr[results.500$osr == "osr0.05"] <- 0.05


hd.mal <- ggplot(results.500, aes(y=1-freq, x=as.factor(comm))) + 
  ylim(0, 1) +
  geom_hline(yintercept = .5, alpha = .5) +
  geom_violin(aes(fill=as.factor(osr)), stat="ydensity", position="dodge", alpha=0.8, trim=TRUE, 
              scale="area") + theme_light() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none") +
  ylab("Allele frequency") +
  ggtitle("HD model - Rare females") +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) 


# ESD (Common sex (females) = 500, s=0.5, h=0.5)

esd.rare.mal <- read.csv("../results/ESD.raremal.csv", 
                         header = TRUE, sep = ",", as.is = T, check.names = F,
                         row.names = 1)

result.plot <- esd.rare.mal[(esd.rare.mal$s == 0.5),]
result.plot <- result.plot[(result.plot$h == 0.5),]
results.500 <- result.plot[(result.plot$num.com == 500),]

esd.rare.mal <- ggplot(results.500, aes(y=freq0, x=as.factor(num.com))) + 
  ylim(0, 1) +
  geom_hline(yintercept = .5, alpha=.5) +
  geom_violin(aes(fill=as.factor(OSR)), stat="ydensity", position="dodge", alpha=0.8, trim=TRUE, 
              scale="area") + 
  theme_light() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none") +
  ylab("Allele frequency") +
  ggtitle("ESD") +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) 

sum(results.500$freq0[results.500$OSR == 0.05] == 1)
sum(results.500$freq0[results.500$OSR == 0.05] == 0)
print(esd.rare.mal)

library(cowplot)
g1 <- ggplotGrob(dominance)
g2 <- ggplotGrob(rd)
g3 <- ggplotGrob(selection)
g4 <- ggplotGrob(autosomal)
g5 <- ggplotGrob(xchrom)
g6 <- ggplotGrob(ychrom)
g7 <- ggplotGrob(hd.fem)
g8 <- ggplotGrob(hd.mal)
g9 <- ggplotGrob(esd.rare.mal)


plot_grid(
  g1, g2, g3, g4, g5, g6, g7, g8, g9, ncol=3,
  labels = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'),
  align="hv"
)