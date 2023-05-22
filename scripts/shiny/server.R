# This code is taken from the evobiR package
# http://coleoguy.github.io/software.html

library(shiny)
library(doParallel)

# Expected frequencies
#TODO

# setwd("C:/Users/pdglenn/OneDrive/TAMU/OSR/ShinyOSRx/server_ui")
common = 1000
s = 0.1
osr = 1.0
h = 0.5
males = 1000
iter = 3
gens = 1000
x = "Sex"
chromosome = 0.2
sex = 1
#input$chromosome

no_cores <- detectCores(logical = TRUE)
# number at the end here determines the number of cores left free
cl <- makeCluster(no_cores-2)
registerDoParallel(cl)

source("functions.R")


#TODO - set up to wear you can select rare male or rare female in the simulation, currently common male?

gens = 1000
# Simulated frequencies
get.simulated.results <- function(common,s, osr,h, iter, chromosome, sex) {
  
    # this sets up the container for a given sim
    resultY <- resultX <- resultA <- c()
    wFem <- wMal <- wDiff <- c()
    
    results <- Aresults <- Xresults <- Yresults <- matrix(,iter,gens)
    
    if (iter > 0) {
      for (k in 1:iter) {
        if(sex == 1){
          males <- common
          females <- ceiling(common * osr)
          if (females < 4) {
            females <- 4
          }
        } else if(sex == 2){
          females <- common
          males <- ceiling(common * osr)
          if (males < 4) {
            males <- 4
          }
        }
        if (chromosome == 0.2) { # sex

        pop <- makeGenomes(females = females, males = males,
                           freqs = c(rep(females / 4, 4),
                                     rep(males / 4, 4)))
        segregating <- TRUE
        c.ite <- 0
        while (segregating) {
          # this gets the allele frequencies
          c.ite <- c.ite + 1
          
          Y <- resultY[c.ite] <- GetFreq(pop, chrom = "Y", allele = 1,
                                         females = females, males = males)
          X <- resultX[c.ite] <- GetFreq(pop, chrom = "X", allele = 2,
                                         females = females, males = males)
          
          # Males
          p2 <- (1 - X) * Y
          q2 <- X * (1 - Y)
          wMal[c.ite] <- ((p2 * (1 + s)) + (((X * Y) + ((1 - X) * (1 - Y))) * (1 + (h * s))) + (q2 * 1))
          
          # Females
          p <- X # p is allele 2, beneficial for females
          q <- 1 - p
          wFem[c.ite] <- ((p^2 * (1 + s)) + (2 * p * q * (1 + ((1 - h) * s))) + (q^2 * 1))
          
          # Difference
          #TODO update the Wdiff
          if(sex == 1){ #common male
            wDiff[c.ite] <- wMal[c.ite] - wFem[c.ite]
          } else{ #common female
            wDiff[c.ite] <- wFem[c.ite] - wMal[c.ite]
          }
          
          
          pop <- Generation(pop, females = females, males = males,
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
        pop <- makeGenomes(females = females, males = males,
                           freqs = c(rep(females / 4, 4),
                                     rep(males / 4, 4)))
        segregating <- TRUE
        c.ite <- 0
        while (segregating) {
          c.ite <- c.ite + 1
          
          A <- resultA[c.ite] <- GetFreq(pop, chrom = "A", allele = 1,
                                         females = females, males = males)
          
          p <- A # p is allele 1, 1 is beneficial for males
          q <- 1 - p # q is allele 2, 2 is beneficial for females
          
          wMal[c.ite] <- ((p^2 * (1 + s)) + (2 * p * q * (1 + h * s)) + (q^2 * 1))
          wFem[c.ite] <- (p^2 * 1) + (2 * p * q * (1 + ((1 - h) * s))) + (q^2 * (1 + s))
          
          if(sex == 1){ #common male
            wDiff[c.ite] <- wMal[c.ite] - wFem[c.ite]
          } else{ #common female
            wDiff[c.ite] <- wFem[c.ite] - wMal[c.ite]
          }
          
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
    input$seed.val
    input$pop
    input$s
    input$osr
    input$h
    input$iter
    input$chromosome
    input$sex
    
    get.simulated.results(
      common = as.numeric(input$pop),
      s = as.numeric(input$s),
      osr = as.numeric(input$osr),
      h = as.numeric(input$h),
      iter = as.numeric(input$iter),
      chromosome = input$chromosome,
      sex = as.numeric(input$sex)
    )
  })
  
  # Expected calculations Heath example
  # expected.A <- reactive({
  #   get.expected.results(input$initial.A, input$gen, input$fit.AA, input$fit.Aa,
  #                        input$fit.aa, input$qAa, input$qaA)
  # })
  
  
  #TODO Find the mean?
  
  
  # Expected calculations
  #TODO
  
  # Plot the results
  output$genePlot <- renderPlot({
    
    if(input$var.plot == 1){ #Wdiff
      # Set up plot area
      plot(0, 0, col = 'white', xlim = c(0, gens), ylim = c(-1,1),
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
      #TODO add lines for expected data
      
      # Plot expected outcome from Heath's example
      # if (input$var.plot == 1) 
      #   lines(1:input$gen, expected.A()^2, lwd = lwd.expected)
      }
    
    # Plot expected outcome
    #TODO
    
  })

})


stopCluster(cl)

      