
rare.fem <- read.csv("../../results/rare.female.250.iter.csv")[,-1]
rare.mal <- read.csv("../../results/rare.male.250.iter.csv")[,-1]
rare.fem.ag <- aggregate(rare.mal, 
                         by=list(Males = rare.fem$males,
                                 OSR = rare.fem$OSR,
                                 Recomb = rare.fem$rd,
                                 h = rare.fem$h,
                                 s = rare.fem$s),
                         FUN = mean)[,1:8]
rare.mal.ag <- aggregate(rare.mal, 
                         by=list(Females = rare.mal$females,
                                 OSR = rare.mal$OSR,
                                 Recomb = rare.mal$rd,
                                 h = rare.mal$h,
                                 s = rare.mal$s),
                         FUN = mean)[,1:8]

rare.mal.ag$Males <- rare.mal.ag$Females * rare.mal.ag$OSR
rare.fem.ag$Females <- rare.fem.ag$Males * rare.fem.ag$OSR
rm(rare.fem,rare.mal)

trial <- rare.mal.ag[rare.mal.ag$h == .5,]
trial <- trial[trial$s == .5,]
trial <- trial[trial$Recomb == .2,]
trial$X <- trial$X - trial$X[28]
trial$Y <- trial$Y - trial$Y[28]

min(trial$X)

plot(x= trial$Females,
     y= trial$Males)


