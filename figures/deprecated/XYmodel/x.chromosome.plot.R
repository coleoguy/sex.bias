# Julio Rincones Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com
# P.I. Dr. Heath Blackmon
# coleoguy@gmail.com

# Script to plot allele frequencies for the sexually antagonistic locus when it is present on
# the X chromosome - XY model.
library(ggplot2)
library(viridis)

# set your working directory to sex.bias
resultsplot <- read.csv("Xchrom.rare.male.csv",as.is=T)
resultsplot <- resultsplot[(resultsplot$rd != 0.5),]
resultsplot <- resultsplot[(resultsplot$h != "h99"),]
resultsplot <- resultsplot[(resultsplot$s == "s0.5"),]
resultsplot <- resultsplot[(resultsplot$comm %in% c(500,100)),]
resultsplot <- resultsplot[(resultsplot$rd == .2), ]

resultsplot$h[resultsplot$h == "h1"] <- 0
resultsplot$h[resultsplot$h == "h0.5"] <- .5
resultsplot$h[resultsplot$h == "h0"] <- 1

p4 <- ggplot(resultsplot, aes(y=freq1, x=as.factor(comm))) + 
  geom_violin(aes(fill=as.factor(osr)), stat="ydensity", position="dodge", 
              alpha=0.5, trim=TRUE, scale="area") + facet_grid(h ~ s) + theme_light() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  guides(fill=guide_legend(title="OSR")) + xlab("Common sex (females)") + 
  ylab("Frequency of the female beneficial allele") + 
  ggtitle("Allele frequency on the X-chromosome") +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) 

ggsave(filename = "xchrom.raremal.s05.pdf",
       plot= p4, width=11, height=10, units="in")
