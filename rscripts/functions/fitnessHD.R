# Julio Rincones Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com
# P.I. Dr. Heath Blackmon
# coleoguy@gmail.com

# functions to assess the impact of allele frequency change in our 
# haplodiploid model in terms of feminisation or masculinisation of the genome

getFitness <- function(data, dom, s) {
  # data data frame of results
  # dom dominance factor for female benefit allele
  # 
  results <- matrix(, nrow(data), 2)
  colnames(results) <- c("wbar.fem", "wbar.mal")
  for (i in 1:nrow(data)) {
    switch(dom,
           "h=0" = p <- data[i, 1],
           "h=0.5" = p <- data[i, 2],
           "h=1" = p <- data[i, 3])
    
      # genotype frequencies
      # females
      fa11 <- p ^ 2
      fa12 <- 2 * p * (1 - p)
      fa22 <- (1 - p) ^ 2
      # males
      ma1 <- p
      ma2 <- 1 - p

    if (dom == "h=0") {
      # Female fitness functions
      wf.a11 <- 1 + s
      wf.a12 <- 1 
      wf.a22 <- 1
      # Male fitness functions
      wm.a1 <- 1 / (1 + s)
      wm.a2 <- 1
      
    }
    if (dom == "h=0.5") {
      # Female fitness functions
      wf.a11 <- 1 + s
      wf.a12 <- 1 + .5 * s
      wf.a22 <- 1
      # Male fitness functions
      wm.a1 <- 1 / (1 + s)
      wm.a2 <- 1    
    }
    if (dom == "h=1") {
      # Female fitness functions
      wf.a11 <- 1 + s
      wf.a12 <- 1 + s
      wf.a22 <- 1
      # Male fitness functions
      wm.a1 <- 1 / (1 + s)
      wm.a2 <- 1      
    }
    # Standardise to maximum possible fitness
    wbar.fem <- ((fa11 * wf.a11) +
                 (fa12 * wf.a12) +
                 (fa22 * wf.a22)) / wf.a11
    wbar.mal <- ((ma1 * wm.a1) + (ma2 * wm.a2)) / wm.a2
    results[i, 1:2] <- c(wbar.fem, wbar.mal)
  }
  return(as.data.frame(results))
}

# Function to get the fitness difference when males are the common sex

getDifference.mal <- function(w.beta) {
  mean(results$wbar.mal - results$wbar.fem)
}

# Function to get the fitness difference when females are the common sex

getDifference.fem <- function(w.beta) {
  mean(results$wbar.fem - results$wbar.mal)
}


