# function to return male and female fitness from XY model results

GetSexFitness <- function(foo){
  fit.r <- as.data.frame(matrix(NA, 0, 2))
  colnames(fit.r) <- c("female.fit", "male.fit")
  for(i in 1:nrow(foo)){
    if(i %% 5000 == 0) print(i)
    s <- foo$s[i]
    h <- foo$h.mal.ben[i]
    # p is male benefit allele the frequency recorded in Y and A columns
    if(foo$r[i] == .5){
      px <- py <- foo$freqA[i]
    }else{
      px <- foo$freqX[i]
      py <- foo$freqY[i]
    }
    # get frequencies
    # m for male f for female 1 for p 2 for q
    m11 <- py*px
    m12 <- py*(1-px) + (1-py)*px
    m22 <- (1-py) * (1-px)
    f11 <- px*px
    f12 <- 2 * px * (1-px)
    f22 <- (1-px) * (1-px) 
    # get fitnesses
    # male fitnesses are naturally normalized to 1
    wm11 <- 1
    if(foo$h.mal.ben[i] != 99){
      wm12 <- 1 / (1 + h*s)
    }else{
      wm12 <- 1 / (1 + s)
    }
    wm22 <- 1 / (1 + s)
    # female fitnesses are normalized to 1
    wf11 <- 1 / (1 + s)
    if(foo$h.mal.ben[i] != 99){
      wf12 <- (1 + h * s) / (1 + s)
    }else{
      wf12 <- 1 / (1 + s)
    }
    wf22 <- (1 + s) / (1 + s)
    avg.wm <- m11*wm11 + m12*wm12 + m22*wm22
    avg.wf <- f11*wf11 + f12*wf12 + f22*wf22
    fit.r[i, 1] <- avg.wf
    fit.r[i, 2] <- avg.wm
  }
  final.results <- data.frame(foo, fit.r)
  return(final.results)
}

