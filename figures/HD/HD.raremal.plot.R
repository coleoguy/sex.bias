# Julio Rincones Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com
# P.I. Dr. Heath Blackmon
# coleoguy@gmail.com

# Script to plot allele frequencies for the haplodiploidy model.

library(ggraptR)
library(viridis)

# Set your working directory to sex.bias and load the data.
hd.rare.mal <- read.csv("parsed.HD.rare.male.csv", 
                        header = TRUE, sep = ",", as.is = T, check.names = F,
                        row.names = 1)

resultsplot <- hd.rare.mal[(hd.rare.mal$h == "h=0.5"),]


p1 <- ggplot(resultsplot, aes(y=freq, x=as.factor(osr))) + 
  geom_violin(aes(fill=as.factor(osr)), stat="ydensity", position="dodge", 
              alpha=0.5, trim=TRUE, scale="area") + facet_grid(h ~ s) + theme_classic() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill=guide_legend(title="osr")) + xlab("Common sex (females)") + 
  ylab("Allele frequency for the female beneficial allele") +
  ggtitle("Haplodiploidy model") + 
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) 

ggsave(filename = "hd.raremal.h05.pdf",
       plot= p1, width=11, height=10, units="in")


resultsplot2 <- hd.rare.mal[(hd.rare.mal$h == "h=1"),]

p2 <- ggplot(resultsplot2, aes(y=freq, x=as.factor(osr))) + 
  geom_violin(aes(fill=as.factor(osr)), stat="ydensity", position="dodge", 
              alpha=0.5, trim=TRUE, scale="area") + facet_grid(h ~ s) + theme_classic() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill=guide_legend(title="osr")) + xlab("Common sex (females)") + 
  ylab("Frequency of the female beneficial allele") + scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) 

ggsave(filename = "hd.raremal.h1.pdf",
       plot= p2, width=11, height=10, units="in")

resultsplot3 <- hd.rare.mal[(hd.rare.mal$h == "h=0"),]

p3 <- ggplot(resultsplot3, aes(y=freq, x=as.factor(osr))) + 
  geom_violin(aes(fill=as.factor(osr)), stat="ydensity", position="dodge", 
              alpha=0.5, trim=TRUE, scale="area") + facet_grid(h ~ s) + theme_classic() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill=guide_legend(title="osr")) + xlab("Common sex (females)") + 
  ylab("Frequency of the female beneficial allele") + scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) 

ggsave(filename = "hd.raremal.h0.pdf",
       plot= p3, width=11, height=10, units="in")

result.h <- hd.rare.mal[(hd.rare.mal$s == "s0.5"),]
results.1000 <- result.h[(result.h$comm == "comm.sex1000"),]
results.100 <- result.h[(result.h$comm == "comm.sex100"),]
result.h <- rbind(results.1000, results.100)

p4 <- ggplot(result.h, aes(y=freq, x=as.factor(comm))) + 
  geom_violin(aes(fill=as.factor(osr)), stat="ydensity", position="dodge", 
              alpha=0.5, trim=TRUE, scale="area") + facet_grid(h ~ s) + theme_light() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill=guide_legend(title="osr")) + xlab("Common sex (females)") + 
  ylab("Allele frequency") +
  ggtitle("Female beneficial allele under a Haplodiploidy model") + 
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) 

ggsave(filename = "hd.raremal.s05.pdf",
       plot= p4, width=11, height=10, units="in", path = "../figures/HD")


result.h <- hd.rare.mal[(hd.rare.mal$s == "s0.5"),]
results.1000 <- result.h[(result.h$comm == "comm.sex1000"),]
results.100 <- result.h[(result.h$comm == "comm.sex100"),]
result.h <- rbind(results.1000, results.100)

p5 <- ggplot(result.h, aes(y=freq, x=as.factor(comm))) + 
  geom_violin(aes(fill=as.factor(osr)), stat="ydensity", position="dodge", 
              alpha=0.5, trim=TRUE, scale="area") + facet_grid(h ~ s) + theme_light() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill=guide_legend(title="OSR")) + xlab("Common sex (females)") + 
  ylab("Allele frequency for the female beneficial allele") +
  ggtitle("Haplodiploidy model") + 
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) 

ggsave(filename = "hd.raremal.s05.pdf",
       plot= p5, width=11, height=10, units="in")

