# This code is taken from the evobiR package
# http://coleoguy.github.io/software.html

library(shiny)
library(doParallel)

# Expected frequencies
#TODO

# setwd("C:/Users/pdglenn/OneDrive/TAMU/OSR/ShinyOSRx/server_ui")
# common = 1000
# s = 0.1
# osr = 1.0
# h = 0.5
# males = 1000
# iter = 3
# gens = 1000
# x = "Sex"
# chromosome = 0.2
# #input$chromosome

no_cores <- detectCores(logical = TRUE)
# number at the end here determines the number of cores left free
cl <- makeCluster(no_cores-2)
registerDoParallel(cl)

source("functions.R")


#TODO - set up to wear you can select rare male or rare female in the simulation, currently common male?

gens = 1000
# Simulated frequencies
get.simulated.results <- function(common,s, osr,h, males, iter, chromosome) {

    females <- ceiling(common * osr)
    if (females < 4) {
      females <- 4
      }
  
    # this sets up the container for a given sim
    resultY <- resultX <- resultA <- c()
    wFem <- wMal <- wDiff <- c()
    
    results <- Aresults <- Xresults <- Yresults <- matrix(,iter,gens)
    
    if (iter > 0) {
      for (k in 1:iter) {
        if (chromosome == 0.2) { # sex
        pop <- makeGenomes(females, males,
                           freqs = c(rep(females / 4, 4),
                                     rep(common / 4, 4)))
        segregating <- TRUE
        c.ite <- 0
        while (segregating) {
          # this gets the allele frequencies
          c.ite <- c.ite + 1
          
          Y <- resultY[c.ite] <- GetFreq(pop, chrom = "Y", allele = 1,
                                         females = females, males = common)
          X <- resultX[c.ite] <- GetFreq(pop, chrom = "X", allele = 2,
                                         females = females, males = common)
          
          # Males
          p2 <- (1 - X) * Y
          q2 <- X * (1 - Y)
          wMal[c.ite] <- ((p2 * (1 + s)) + (((X * Y) + ((1 - X) * (1 - Y))) * (1 + (h * s))) + (q2 * 1))
          
          # Females
          p <- X # p is allele 2, beneficial for females
          q <- 1 - p
          wFem[c.ite] <- ((p^2 * (1 + s)) + (2 * p * q * (1 + ((1 - h) * s))) + (q^2 * 1))
          
          # Difference
          wDiff[c.ite] <- wMal[c.ite] - wFem[c.ite]
          
          pop <- Generation(pop, females = females, males = common,
                            rd = 0.2, h = h, s = s)
          
          # Check if we have run gens generations
          if (c.ite == gens) {
            segregating <- FALSE
          }
        }
        results[k,] <- wDiff
        Xresults[k,] <- resultX
        Yresults[k,] <- resultY
        
        total <- list(results,Xresults,Yresults)
        names(total) <- c("wDiff","Xfreq","Yresults")
        
      } else { # auto
        pop <- makeGenomes(females, males,
                           freqs = c(rep(females / 4, 4),
                                     rep(common / 4, 4)))
        segregating <- TRUE
        c.ite <- 0
        while (segregating) {
          c.ite <- c.ite + 1
          
          A <- resultA[c.ite] <- GetFreq(pop, chrom = "A", allele = 1,
                                         females = females, males = pop)
          
          p <- A # p is allele 1, 1 is beneficial for males
          q <- 1 - p # q is allele 2, 2 is beneficial for females
          
          wMal[c.ite] <- ((p^2 * (1 + s)) + (2 * p * q * (1 + h * s)) + (q^2 * 1)) # /(1+x[8])
          wFem[c.ite] <- (p^2 * 1) + (2 * p * q * (1 + ((1 - h) * s))) + (q^2 * (1 + s))
          wDiff[c.ite] <- wMal[c.ite] - wFem[c.ite]
          
          pop <- Generation(pop, females = females, males = males,
                            rd = 0.5, h = h, s = s)
          
          # Check if we have run gens generations
          if (c.ite == gens) {
            segregating <- FALSE
          }
        }
        results[k,] <- wDiff
        Aresults[k,] <-resultA
        
        total <- list(results,Aresults)
        names(total) <- c("wDiff","Afreq")
      }
        
    }
  }
  return(total)
}


shinyServer(function(input, output) {
  
  # Drift simulations
  data <- reactive({
    set.seed(input$seed.val)
    get.simulated.results(common = as.numeric(input$pop),s = as.numeric(input$s),osr = as.numeric(input$osr),h = as.numeric(input$h),
                                 males = as.numeric(input$pop),iter = as.numeric(input$iter),chromosome = input$chromosome)
  })
  
  
  # Expected calculations
  #TODO
  
  # Plot the results
  output$genePlot <- renderPlot({
    
    if(input$var.plot == 1){ #Wdiff
      # Set up plot area
      plot(0, 0, col = 'white', xlim = c(0, gens), ylim = c(-0.5,0.5),
           xlab = 'Time (generations)', ylab = "Wdiff (common-rare)",
           cex.lab=1.5, cex.axis=1.3, main="Fitness difference")
      mtext("Produced with the package evobiR", 
            side = 1, cex=.8, line=4, adj=1)
    } else { #Allele frequences
      # Set up plot area
      plot(0, 0, col = 'white', xlim = c(0, gens), ylim = c(0,1),
           xlab = 'Time (generations)', ylab = "Allele Frequency",
           cex.lab=1.5, cex.axis=1.3, main="Allele Frequency")
      mtext("Produced with the package evobiR", 
            side = 1, cex=.8, line=4, adj=1)
    }
    

    rep <- as.numeric(input$iter)
    # Plot simulated outcomes
    if (rep > 0) {
      
      if(input$chromosome == 0.5){ #Auto section
      
        for (i in 1:rep) {
          if (input$var.plot == 1) { #Wdiff
            lines(
                1:gens,
                data()[[1]][i,],
                col = rainbow(rep)[i]
              )
          } else if (input$var.plot == 2) { #AutoFreq
            lines(
              1:gens,
              data()[[2]][i,],
              col = rainbow(rep)[i]
            )
            }
          }
        } else if(input$chromosome == 0.2){ #Sex Chr
          for (i in 1:rep) {
          
             if (input$var.plot == 1) { #Wdiff
              lines(
                1:gens,
                data()[[1]][i,],
                col = rainbow(rep)[i]
              )
              
            } else if (input$var.plot == 3) { #Xfreq
              lines(
                1:gens,
                data()[[2]][i,],
                col = rainbow(rep)[i]
              )
              
            } else if (input$var.plot == 4) { #Yfreq
              lines(
                1:gens,
                data()[[3]][i,],
                col = rainbow(rep)[i]
              )
            }
          }
        }
      }
    
    # Plot expected outcome
    #TODO
    
  })

})


stopCluster(cl)

      