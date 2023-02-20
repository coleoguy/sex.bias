###Used to generate viridis plots for OSR paper
##need to load data from HB.autosomes
library(gridExtra)
library(grid)
library(ggpubr)
library(ggplot2)
library(viridis)

foo.prob.ben$common.num<-as.factor(foo.prob.ben$common.num)
g1<-ggplot(foo.prob.ben, aes(y=value, x=OSR))+geom_line(aes(colour=common.num), size=1)+
  geom_point(aes(shape=common.num, fill=common.num), stat="identity", position="identity", size=3)+
  scale_color_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~ .)+theme_bw()+ theme(text=element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+theme(axis.title.y = element_blank())+
  scale_size(range=c(1,4))+xlab("Operational Sex Ratio")+ylab("Proporion of simulations with fixation")+labs("Number of common sex")+ggtitle("Proportion with fixation")+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
g1
foo.freq$common.num<-as.factor(foo.freq$common.num)
g2<-ggplot(foo.freq, aes(y=value, x=OSR))+geom_line(aes(colour=common.num), size=1)+
  geom_point(aes(shape=common.num, fill=common.num), stat="identity", position="identity", size=3)+
  scale_color_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+scale_fill_viridis_d(end=0.9)+
  facet_grid(h~ .)+theme_bw()+ theme(text=element_text(family="sans", face="plain", color="#000000", size=13, hjust=0.5, vjust=0.5))+theme(axis.title.y = element_blank())+
  scale_size(range=c(1,4))+xlab("Operational Sex Ratio")+ggtitle("Frequency of allele")+ylim(c(0.00,1.00))+
  labs(colour="Number of the\n common sex", shape="Number of the\n common sex", fill="Number of the\n common sex")
g2

final<-ggarrange(g1, g2, ncol=2, legend="bottom", common.legend=TRUE)
final


###for the next plot

gen.avgd$common<-as.factor(gen.avgd$common)
ggplot(gen.avgd, aes(y=fitness, x=OSR)) + 
  geom_line(aes(colour=common), size=1)+geom_point(aes(shape=common, fill=common), stat="identity", position="identity", alpha=1, size=4)  +
  scale_colour_viridis_d(end=0.9)+scale_shape_manual(values=c(21,24,22,23))+ scale_fill_viridis_d(end=0.9)+
  theme_bw() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  scale_size(range=c(1, 3)) + xlab("Operational Sex Ratio") + ylab("fitness difference (common - rare)")+
  labs(colour="number of the\n common sex", shape="number of the\n common sex", fill="number of the\n common sex")+
  theme(legend.position = "right")
str(gen.avgd)
