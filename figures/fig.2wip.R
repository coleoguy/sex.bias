# Julio Rincones Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com
# P.I. Dr. Heath Blackmon
# coleoguy@gmail.com

# Script to plot figure 2 for the manuscript
library(ggplot2)
library(viridis)


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
        legend.position = "bottom",
        legend.background = element_rect(fill= NA),
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5,)) + 
  #guides(fill=guide_legend(title="Dominance factor")) + 
  ylab("Allele frequency") +
  ggtitle("Effect of genetic architecture") +
  scale_fill_viridis(discrete=TRUE, option = "B") +
  scale_color_viridis(discrete=TRUE, option = "B")

print(dominance)

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

print(rd)


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
        legend.position = "bottom",
        legend.background = element_rect(fill= NA),
        legend.key.size = unit(0.1, "cm"),
        legend.key.width = unit(0.25,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill= guide_legend(title= "OSR")) +
  ylab("Allele frequency") +
  ggtitle("Autosomal") +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE)

print(autosomal)

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
        legend.position = c(0.935, 0.25),
        legend.background = element_rect(fill= NA),
        legend.key.size = unit(0.1, "cm"),
        legend.key.width = unit(0.25,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill= guide_legend(title= "OSR")) + xlab("Common sex (females)") + 
  ylab("Allele frequency") +
  ggtitle("X chromosome") +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE)

print(xchrom)

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
        legend.position = c(0.935, 0.25),
        legend.background = element_rect(fill= NA),
        legend.key.size = unit(0.1, "cm"),
        legend.key.width = unit(0.25,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill= guide_legend(title= "OSR")) + xlab("Common sex (Males)") + 
  ylab("Allele frequency") +
  ggtitle("Y chromosome") +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE)

print(ychrom)

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

haplodiploidy.fem <- ggplot(results.500, aes(y=freq, x=as.factor(comm))) +
  ylim(0, 1) +
  geom_hline(yintercept = .5, alpha = .5) +
  geom_violin(aes(fill=as.factor(osr)), stat="ydensity", position="dodge", alpha=0.8, trim=TRUE, 
              scale="area") + theme_light() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = c(0.935, 0.25),
        legend.background = element_rect(fill= NA),
        legend.key.size = unit(0.1, "cm"),
        legend.key.width = unit(0.25,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill= guide_legend(title= "OSR")) + xlab("Common sex (females)") + 
  ylab("Allele frequency") +
  ggtitle("HD model - Rare males") +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) 

print(haplodiploidy.fem)

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


haplodiploidy.mal <- ggplot(results.500, aes(y=1-freq, x=as.factor(comm))) + 
  ylim(0, 1) +
  geom_hline(yintercept = .5, alpha = .5) +
  geom_violin(aes(fill=as.factor(osr)), stat="ydensity", position="dodge", alpha=0.8, trim=TRUE, 
              scale="area") + theme_light() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = c(0.935, 0.25),
        legend.background = element_rect(fill= NA),
        legend.key.size = unit(0.1, "cm"),
        legend.key.width = unit(0.25,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill= guide_legend(title= "OSR")) + xlab("Common sex (Males)") + 
  ylab("Allele frequency") +
  ggtitle("HD model - Rare females") +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) 

print(haplodiploidy.mal)

# ESD (Common sex (females) = 500, s=0.5, h=0.5)

load("../results/esd.RData")
esd.rare.mal <- results

result.plot <- esd.rare.mal[(esd.rare.mal$s == 0.5),]
result.plot <- result.plot[(result.plot$h == 0.5),]
results.500 <- result.plot[(result.plot$num.com == 500),]

esd <- ggplot(results.500, aes(y=freq0, x=as.factor(num.com))) + 
  ylim(0, 1) +
  geom_hline(yintercept = .5, alpha=.5) +
  geom_violin(aes(fill=as.factor(OSR)), stat="ydensity", position="dodge", alpha=0.8, trim=TRUE, 
              scale="area") + 
  theme_light() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = c(0.935, 0.25),
        legend.background = element_rect(fill= NA),
        legend.key.size = unit(0.1, "cm"),
        legend.key.width = unit(0.25,"cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill= guide_legend(title= "OSR")) +
  ylab("Allele frequency") +
  ggtitle("ESD") +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) 

print(esd)

# use grid arrange to combine all of these plots into a single figure. Number of rows = 3, so there  
# will be three plots per row.

# fig.2 <- grid.arrange(dominance, rd, selection, autosomal, xchrom, ychrom, haplodiploidy.fem,
#                       haplodiploidy.mal, esd, nrow = 3, respect = TRUE)
# 
# 
# library(ggpubr)
# fig.2 <- ggarrange(dominance, rd, selection, autosomal, xchrom, ychrom, 
#                    haplodiploidy.fem, haplodiploidy.mal, esd, ncol = 3, nrow = 3,
#                    labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"))
# 
# library(grid)
# lay <- rbind(c(1,2,3),
#              c(4,5,6),
#              c(7,8,9))
# 
# 
# fig.2 <- print(grid.arrange(arrangeGrob(dominance, left = textGrob("a)", x = unit(1, "npc"), 
#                                                        y = unit(1, "npc"))), 
#                    arrangeGrob(rd, left =textGrob("b)", x = unit(1, "npc"), 
#                                                       y = unit(1, "npc"))),
#                    arrangeGrob(selection, left=textGrob("c)", x = unit(1, "npc"), 
#                                                      y = unit(1, "npc"))),
#                    arrangeGrob(autosomal, left=textGrob("d)", x = unit(1, "npc"), 
#                                                         y = unit(1, "npc"))),
#                    arrangeGrob(xchrom, left=textGrob("e)", x = unit(1, "npc"), 
#                                                         y = unit(1, "npc"))),
#                    arrangeGrob(ychrom, left=textGrob("f)", x = unit(1, "npc"), 
#                                                         y = unit(1, "npc"))),
#                    arrangeGrob(haplodiploidy.fem, left=textGrob("g)", x = unit(1, "npc"), 
#                                                         y = unit(1, "npc"))),
#                    arrangeGrob(haplodiploidy.mal, left=textGrob("h)", x = unit(1, "npc"), 
#                                                         y = unit(1, "npc"))),
#                    arrangeGrob(esd, left=textGrob("i)", x = unit(1, "npc"), 
#                                                         y = unit(1, "npc"))),
#                    layout_matrix = lay))
# 
# library(gtable)
# g1 <- ggplotGrob(dominance)
# g2 <- ggplotGrob(rd)
# g3 <- ggplotGrob(selection)
# g4 <- ggplotGrob(autosomal)
# g5 <- ggplotGrob(xchrom)
# g6 <- ggplotGrob(ychrom)
# g7 <- ggplotGrob(haplodiploidy.fem)
# g8 <- ggplotGrob(haplodiploidy.mal)
# g9 <- ggplotGrob(esd)
# 
# g <- rbind(g1, g2, g3, g4, g5, g6, g7, g8, g9, size ="max")
# 
# g$widths <- unit.pmax(g1$widths, g2$widths, g3$widths, g4$widths, g5$widths, g6$widths, g7$widths,
#                         g8$widths, g9$widths)
# grid.newpage()
# grid.draw(g)


library(cowplot)
g1 <- ggplotGrob(dominance)
g2 <- ggplotGrob(rd)
g3 <- ggplotGrob(selection)
g4 <- ggplotGrob(autosomal)
g5 <- ggplotGrob(xchrom)
g6 <- ggplotGrob(ychrom)
g7 <- ggplotGrob(haplodiploidy.fem)
g8 <- ggplotGrob(haplodiploidy.mal)
g9 <- ggplotGrob(esd)


plot_grid(
  g1, g2, g3, g4, g5, g6, g7, g8, g9, ncol=3,
  labels = c('A)', 'B)', 'C)', 'D)', 'E)', 'F)', 'G)', 'H)', 'I)'),
  align="hv"
)



ggsave(filename = "fig2.final.pdf",
       plot= fig.2 , width=11, height=13, units="cm")
