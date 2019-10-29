# Julio Rincones Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com
# P.I. Dr. Heath Blackmon
# coleoguy@gmail.com

# Script to plot allele frequencies for the ESD model.
library(ggraptR)
library(viridis)

# Set your working directory to be sex.bias
load("../../../sex.bias/rscripts/results/esd.RData")
esd.rare.mal <- results

resultsplot <- esd.rare.mal[(esd.rare.mal$h ==0.5),]


p1 <- ggplot(resultsplot, aes(y=freq0, x=as.factor(num.com))) + 
  geom_violin(aes(fill=as.factor(OSR)), stat="ydensity", position="dodge", 
              alpha=0.5, trim=TRUE, scale="area") + facet_grid(h ~ s) + theme_classic() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill=guide_legend(title="OSR")) + xlab("Common sex (females)") + 
  ylab("Frequency of the male beneficial allele") + scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) 

ggsave(filename = "esd.raremal.h05.pdf",
       plot= p1, width=11, height=10, units="in")

resultsplot2 <- resultsplot[(resultsplot$s >= 0.5),]

p2 <- ggplot(resultsplot2, aes(y=freq0, x=as.factor(num.com))) + 
  geom_violin(aes(fill=as.factor(OSR)), stat="ydensity", position="dodge", 
              alpha=0.5, trim=TRUE, scale="area") + facet_grid(h ~ s) + theme_classic() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill=guide_legend(title="OSR")) + xlab("Common sex (females)") + 
  ylab("Frequency of the male beneficial allele") + scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) 

ggsave(filename = "esd.raremal.h05.s>=.5.pdf",
       plot= p2, width=11, height=10, units="in")

resultsplot3 <- resultsplot[(resultsplot$s == 0.1),]
resultsplot4 <- resultsplot[(resultsplot$s == 0.9),]
rplot <- rbind(resultsplot3, resultsplot4)

p3 <- ggplot(rplot, aes(y=freq0, x=as.factor(num.com))) + 
  geom_violin(aes(fill=as.factor(OSR)), stat="ydensity", position="dodge", 
              alpha=0.5, trim=TRUE, scale="area") + facet_grid(h ~ s) + theme_classic() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill=guide_legend(title="OSR")) + xlab("Common sex (females)") + 
  ylab("Frequency of the male beneficial allele") + scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) 

ggsave(filename = "esd.raremal.h05.s=0.1&0.9.pdf",
       plot= p3, width=11, height=10, units="in")

result.h <- esd.rare.mal[(esd.rare.mal$s == 0.5),]
results.1000 <- result.h[(result.h$num.com == 1000),]
results.100 <- result.h[(result.h$num.com == 100),]
result.h <- rbind(results.1000, results.100)

p4 <- ggplot(result.h, aes(y=freq0, x=as.factor(num.com))) + 
  geom_violin(aes(fill=as.factor(OSR)), stat="ydensity", position="dodge", 
              alpha=0.5, trim=TRUE, scale="area") + facet_grid(h ~ s) + theme_light() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill=guide_legend(title="OSR")) + xlab("Common sex (females)") + 
  ylab("Frequency of the female beneficial allele") + 
  ggtitle("Environmental sex determination model") +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) 

ggsave(filename = "esd.raremal.s05.pdf",
       plot= p4, width=11, height=10, units="in")
