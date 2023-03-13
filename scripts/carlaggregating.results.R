#Read in files
rare.fem <- read.csv("rare.female.250.iter.csv") #this is expecting column 1 to be row Ns from previously unfixed script
rare.mal <- read.csv("rare.male.250.iter.csv")

#Aggregate data
rare.fem.ag <- aggregate(rare.fem, 
                         by=list(Males = rare.fem$males,
                                 OSR = rare.fem$OSR,
                                 Recomb = rare.fem$rd,
                                 h = rare.fem$h,
                                 s = rare.fem$s),
                         mean)[,1:8]
rare.mal.ag <- aggregate(rare.mal, 
                         by=list(Females = rare.mal$females,
                                 OSR = rare.mal$OSR,
                                 Recomb = rare.mal$rd,
                                 h = rare.mal$h,
                                 s = rare.mal$s),
                         FUN = mean)[,1:8]

#Add number of rare sex to data
rare.mal.ag$Males <- rare.mal.ag$Females * rare.mal.ag$OSR
rare.fem.ag$Females <- rare.fem.ag$Males * rare.fem.ag$OSR

#Clean up variables, remove original data frames
rm(rare.fem,rare.mal)


#Trial subset
trial <- rare.mal.ag[rare.mal.ag$h == .5,]
trial <- trial[trial$s == .5,]
trial <- trial[trial$Recomb == .2,]

#subtracting smallest number from each sample to set it starting at 0
trial$X <- trial$X - trial$X[28]
trial$Y <- trial$Y - trial$Y[28]

min(trial$X)

#Checking that we have seven for each trial
plot(x= trial$Females,
     y= trial$Males)



