# Julio Rincones Gamboa
# jgamboa@bio.tamu.edu
# j.a.r.gamboa@gmail.com
# P.I. Dr. Heath Blackmon
# coleoguy@gmail.com

# functions to assess the impact of allele frequency change in our haplodiploid model in terms of
# feminisation or masculinisation of the genome

getFitness <- function(data, dom, h, s) {
  results <- matrix(, nrow(data), 2)
  colnames(results) <- c("wbar.fem", "wbar.mal")
  for (i in 1:nrow(data)) {
    if (dom == "h=0") {
      # we pull the frequency of the female beneficial allele and store it in "p".
      # if you want the frequency of the male beneficial allele this is 1-p
      p <- data[i, 1]
      
      # Female fitness functions
      wf.a11 <- 1 + s
      wf.a12 <- 1 + h * s
      wf.a22 <- 1
      # Male fitness functions
      wm.a1 <- 1 / (1 + s)
      wm.a2 <- 1
      
      # genotype frequencies
      # females
      fa11 <- p ^ 2
      fa12 <- 2 * p * (1 - p)
      fa22 <- (1 - p) ^ 2
      # males
      ma1 <- p
      ma2 <- 1 - p
      
      # Standardise
      wbar.fem <- (((p ^ 2) * wf.a11) +
                     ((2 * p * (1 - p)) * wf.a12) +
                     (1 - p ^ 2 * wf.a22)) / wf.a11
      
      wbar.mal <- ((p * wm.a1) + (1 - p * wm.a2)) / wm.a2
    }
    if (dom == "h=0.5") {
      # we pull the frequency of the female beneficial allele and store it in "p".
      # if you want the frequency of the male beneficial allele this is 1-p
      p <- data[i, 2]
      
      # Female fitness functions
      wf.a11 <- 1 + s
      wf.a12 <- 1 + h * s
      wf.a22 <- 1
      # Male fitness functions
      wm.a1 <- 1 / (1 + s)
      wm.a2 <- 1
      
      # genotype frequencies
      # females
      fa11 <- p ^ 2
      fa12 <- 2 * p * (1 - p)
      fa22 <- (1 - p) ^ 2
      # males
      ma1 <- p
      ma2 <- 1 - p
      
      # Standardise
      wbar.fem <- (((p ^ 2) * wf.a11) +
                     ((2 * p * (1 - p)) * wf.a12) +
                     (1 - p ^ 2 * wf.a22)) / wf.a11
      
      wbar.mal <- ((p * wm.a1) + (1 - p * wm.a2)) / wm.a2
    }
    if (dom == "h=1") {
      # we pull the frequency of the female beneficial allele and store it in "p".
      # if you want the frequency of the male beneficial allele this is 1-p
      p <- data[i, 3]
      
      # Female fitness functions
      wf.a11 <- 1 + s
      wf.a12 <- 1 + h * s
      wf.a22 <- 1
      # Male fitness functions
      wm.a1 <- 1 / (1 + s)
      wm.a2 <- 1
      
      # genotype frequencies
      # females
      fa11 <- p ^ 2
      fa12 <- 2 * p * (1 - p)
      fa22 <- (1 - p) ^ 2
      # males
      ma1 <- p
      ma2 <- 1 - p
      
      # Standardise
      wbar.fem <- (((p ^ 2) * wf.a11) +
                     ((2 * p * (1 - p)) * wf.a12) +
                     (1 - p ^ 2 * wf.a22)) / wf.a11
      
      wbar.mal <- ((p * wm.a1) + (1 - p * wm.a2)) / wm.a2
    }
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


