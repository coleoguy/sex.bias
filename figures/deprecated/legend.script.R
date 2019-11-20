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
dom.leg <- ggplot(resultsplot3, aes(y=frequency, x=as.factor(h))) +
  ylim(0, 1) +
  geom_hline(yintercept = .5, alpha = .5) +
  geom_violin(aes(fill=as.factor(h)), stat="ydensity", position="dodge", 
              alpha=0.8, trim=TRUE, scale="area") +
  theme_light() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom", 
        legend.background = element_rect(fill= NA),
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5,)) +
  guides(fill=guide_legend(title="Dominance factor")) +
  ylab("Allele frequency") +
  ggtitle("Effect of genetic architecture") +
  scale_fill_viridis(discrete=TRUE, option = "B") +
  scale_color_viridis(discrete=TRUE, option = "B")

d.leg <- get_legend(dom.leg)

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

rd.legend <- ggplot(resultsplot2, aes(y=1-frequency, x=as.factor(rd))) +
  ylim(0, 1) +
  geom_hline(yintercept = .5, alpha = .5) +
  geom_violin(aes(fill=as.factor(rd)), stat="ydensity", position="dodge", 
              alpha=0.8, trim=TRUE, scale="area") +
  theme_light() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom",
        legend.background = element_rect(fill= NA),
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5,)) +
  guides(fill=guide_legend(title="Recombination distance")) +
  ylab("Allele frequency") +
  ggtitle("Effect of recombination") +
  scale_fill_viridis(discrete=TRUE, option = "B") +
  scale_color_viridis(discrete=TRUE, option = "B")

rd.leg <- get_legend(rd.legend)


# Facet by selection coefficients but only plot h=0.5. Comm.sex = 100. OSR = 0.2

resultsplot <- read.csv("../../sex.bias/rscripts/parsing/Autosome.rare.male.csv", row.names = 1, 
                        as.is = T, header = T, check.names = F)
resultsplot <- resultsplot[(resultsplot$rd == 0.5),]
resultsplot <- resultsplot[(resultsplot$h == "h0.5"),]
resultsplot <- resultsplot[(resultsplot$comm == 500),]
colnames(resultsplot) <- c("frequency", "common.sex", "osr", "rd", "h", "s")

resultsplot4 <- resultsplot[(resultsplot$osr == 0.1),]

selection.legend <- ggplot(resultsplot4, aes(y=1-frequency, x=as.factor(common.sex))) + 
  ylim(0, 1) +
  geom_hline(yintercept = .5, alpha = .5) +
  geom_violin(aes(fill=as.factor(s)), stat="ydensity", position="dodge", 
              alpha=0.8, trim=TRUE, scale="area") + 
  theme_light() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom",
        legend.background = element_rect(fill= NA),
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) +
  guides(fill=guide_legend(title="selection")) +
  ylab("Allele frequency") +
  ggtitle("Effect of selection") +
  scale_fill_viridis(discrete=TRUE, option = "B") +
  scale_color_viridis(discrete=TRUE, option = "B")

s.leg <- get_legend(selection.legend)


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

auto.legend <- ggplot(resultsplot2, aes(y=1-frequency, x=as.factor(common.sex))) + 
  ylim(0, 1) +
  geom_hline(yintercept = .5, alpha = .5) +
  geom_violin(aes(fill=as.factor(osr)), stat="ydensity", position="dodge", 
              alpha=0.8, trim=TRUE, scale="area") + 
  theme_light() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom",
        legend.background = element_rect(fill= NA),
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        text=element_text(family="sans", face="plain", color="#000000", 
                          size=15, hjust=0.5, vjust=0.5)) +
  guides(fill= guide_legend(title= "OSR")) +
  ylab("Allele frequency") +
  ggtitle("Autosomal") +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE)

au.leg <- get_legend(auto.legend)

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

xchrom.legend <- ggplot(resultsplot2, aes(y=frequency, x=as.factor(common.sex))) +
  ylim(0, 1) +
  geom_hline(yintercept = .5, alpha = .5) +
  geom_violin(aes(fill=as.factor(osr)), stat="ydensity", position="dodge", 
              alpha=0.8, trim=TRUE, scale="area") + 
  theme_light() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "right",
        legend.background = element_rect(fill= NA),
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) +
  guides(fill= guide_legend(title= "OSR")) + xlab("Common sex (females)") +
  ylab("Allele frequency") +
  ggtitle("X chromosome") +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE)
print(xchrom.legend)


x.leg <- get_legend(xchrom.legend)
OSR.leg <- get_legend(xchrom.legend)

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

ychrom.legend <- ggplot(resultsplot2, aes(y=frequency, x=as.factor(common.sex))) + 
  ylim(0, 1) +
  geom_hline(yintercept = .5, alpha = .5) +
  geom_violin(aes(fill=as.factor(osr)), stat="ydensity", position="dodge", 
              alpha=0.8, trim=TRUE, scale="area") + 
  theme_light() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom",
        legend.background = element_rect(fill= NA),
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) +
  guides(fill= guide_legend(title= "OSR")) + xlab("Common sex (Males)") +
  ylab("Allele frequency") +
  ggtitle("Y chromosome") +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE)

y.leg <- get_legend(ychrom.legend)

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

hd.fem.legend <- ggplot(results.500, aes(y=freq, x=as.factor(comm))) +
  ylim(0, 1) +
  geom_hline(yintercept = .5, alpha = .5) +
  geom_violin(aes(fill=as.factor(osr)), stat="ydensity", position="dodge", alpha=0.8, trim=TRUE, 
              scale="area") + theme_light() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom",
        legend.background = element_rect(fill= NA),
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13), 
        text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) +
  guides(fill= guide_legend(title= "OSR")) + xlab("Common sex (females)") +
  ylab("Allele frequency") +
  ggtitle("HD model - Rare males") +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) 


hd.fem.leg <- get_legend(hd.fem.legend)

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


hd.mal.legend <- ggplot(results.500, aes(y=1-freq, x=as.factor(comm))) + 
  ylim(0, 1) +
  geom_hline(yintercept = .5, alpha = .5) +
  geom_violin(aes(fill=as.factor(osr)), stat="ydensity", position="dodge", alpha=0.8, trim=TRUE, 
              scale="area") + theme_light() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom",
        legend.background = element_rect(fill= NA),
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13), 
        text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) +
  guides(fill= guide_legend(title= "OSR")) + xlab("Common sex (Males)") +
  ylab("Allele frequency") +
  ggtitle("HD model - Rare females") +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) 


hd.mal.leg <- get_legend(hd.mal.legend)


# ESD (Common sex (females) = 500, s=0.5, h=0.5)

esd.rare.mal <- read.csv("../results/ESD.raremal.csv", 
                         header = TRUE, sep = ",", as.is = T, check.names = F,
                         row.names = 1)

result.plot <- esd.rare.mal[(esd.rare.mal$s == 0.5),]
result.plot <- result.plot[(result.plot$h == 0.5),]
results.500 <- result.plot[(result.plot$num.com == 500),]

esd.rare.mal.legend <- ggplot(results.500, aes(y=freq0, x=as.factor(num.com))) + 
  ylim(0, 1) +
  geom_hline(yintercept = .5, alpha=.5) +
  geom_violin(aes(fill=as.factor(OSR)), stat="ydensity", position="dodge", alpha=0.8, trim=TRUE, 
              scale="area") + 
  theme_light() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom",
        legend.background = element_rect(fill= NA),
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13), 
        text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) +
  guides(fill= guide_legend(title= "OSR")) +
  ylab("Allele frequency") +
  ggtitle("ESD") +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) 

esd.leg <- get_legend(esd.rare.mal.legend)


grid.arrange(d.leg, rd.leg, s.leg, au.leg, OSR.leg,
             ncol = 2, nrow = 3)
