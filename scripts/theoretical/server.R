# This code is taken from the evobiR package
# http://coleoguy.github.io/software.html

library(shiny)
flower <- function(x, y, colx, cex.x, cex.y){
  petalsx <- cex.x * c(0,   1,  0, -1)
  petalsy <- cex.y * c(1,   0, -1,  0)
  for(i in 1:4){
    points(x + petalsx[i], y + petalsy[i], pch = 16, col = colx)
  }
  points(x, y, pch = 16, col = "white", cex=.5)
}
genotype.options <- c('A1 A1', 'A1 A2', 'A2 A2', 'A1', 'A2')

# Expected frequencies
get.expected.results <- function(x, y, wAA, wAa, waa, qAa, qaA) {

  freqA <- vector()
  A <- x

  # genotype frequencies
  AA <- A^2       # p^2
  Aa <- 2*A*(1-A) # 2pq
  aa <- (1-A)^2   # q^2

  # iterate over generations
  for (i in 1:y) {
    w.bar <- AA*wAA + Aa*wAa + aa*waa  # mean fitness

    # genotype frequencies after selection
    AA <- AA * (wAA / w.bar)  
    Aa <- Aa * (wAa / w.bar)
    aa <- aa * (waa / w.bar)

    # frequency of allele A after selection
    freqA[i] <- A <- (AA + .5*Aa)

    # mutation
    if (qAa + qaA != 0)
        A <- A + ((1 - A) * qaA) - (A * qAa)

    # genotype frequencies after reproduction (random mating)
    AA <- A^2
    Aa <- 2*A*(1-A)
    aa <- (1-A)^2
  }

  # frequency of allele A over time
  return(freqA)
}
# Or instead of all the steps above, could jump straight to p'.

# Simulated frequencies
get.simulated.results <- function(fitness, initial.A, pop, gen, var.plot, 
                                  iter, heath, qAa, qaA) {

  results <- matrix(,iter,gen)
  pop2 <- vector()
  if (iter > 0)
  {
    for (k in 1:iter) {
      adults <- c(rep(1, each = round(pop*initial.A^2)), 
                  rep(2, each = round(pop*2*initial.A*{1-initial.A})), 
                  rep(3, each = round(pop*{1-initial.A}^2)))
      plot.val <- vector()
      for (i in 1:gen) {
        A <- (2 * sum(adults == 1) + sum(adults ==2)) / (pop*2)
        if (qAa + qaA != 0) A <- A + ((1 - A) * qaA) - (A * qAa)
        pop2 <- pop
        babies <-  c(rep(1, each = round(pop2*A^2)), 
                     rep(2, each = round(pop2*2*A*{1-A})), 
                     rep(3, each = round(pop2*(1-A)^2)))
        pop.fit <- vector(length = length(babies)) # fitness for each offspring
        pop.fit[babies == 1] <- fitness[1]
        pop.fit[babies == 2] <- fitness[2]
        pop.fit[babies == 3] <- fitness[3]
        # the sample() function is where the "randomness" really happens
        adults <- sample(babies, pop2, replace = T, prob = pop.fit)
        AA <- sum(adults == 1)
        Aa <- sum(adults == 2)
        plot.val[i] <- AA + .5 * Aa
      }
      results[k,] <- plot.val
    }
  }
  return(results)
}

shinyServer(function(input, output) {

  genotypes <- reactive({
    paste('Frequency of', genotype.options[as.numeric(input$var.plot)])
  })
  
  # Drift simulations
  data <- reactive({
    set.seed <- input$seed.val
    get.simulated.results(fitness=c(input$fit.AA, input$fit.Aa, input$fit.aa), 
                                initial.A = input$initial.A, 
                                pop = input$pop, 
                                gen = input$gen, 
                                var.plot = input$var.plot, 
                                iter = input$iter,
                                heath = input$heath, input$qAa, input$qaA)
  })
  # "data" now contains the number of A alleles in the population
  
  # Expected calculations
  expected.A <- reactive({
    get.expected.results(input$initial.A, input$gen, input$fit.AA, input$fit.Aa,
                         input$fit.aa, input$qAa, input$qaA)
  })
  
  # Plot the results
  output$genePlot <- renderPlot({

      # Set up plot area
      plot(0, 0, col = 'white', ylim = c(0, 1), xlim = c(0, input$gen),
           xlab = 'Time (generations)', ylab = genotypes(), 
           cex.lab=1.5, cex.axis=1.3, main="Change in Allele or Genotype Frequency")
      mtext("Produced with the package evobiR", 
            side = 1, cex=.8, line=4, adj=1)

      lwd.expected <- 5
      lwd.simulated <- 3

      # Plot drift outcomes
      if (input$iter > 0) {
        for (i in 1:input$iter) {
          if (input$var.plot == 1) {
            lines(1:input$gen, 
                  (data()[i,1:input$gen]/input$pop)^2,
                  col=rainbow(input$iter)[i],
                  lwd = lwd.simulated)
            
          } else if (input$var.plot == 2) {
            lines(1:input$gen, 
                  2 * (data()[i,1:input$gen]/input$pop) *
                    (1-(data()[i,1:input$gen]/input$pop)),
                  col=rainbow(input$iter)[i], 
                  lwd = lwd.simulated)
            
          } else if (input$var.plot == 3) {
            lines(1:input$gen, 
                  (1-(data()[i,1:input$gen]/input$pop))^2,
                  col=rainbow(input$iter)[i], 
                  lwd = lwd.simulated)
            
          } else if (input$var.plot == 4) {
            lines(1:input$gen, 
                  data()[i,1:input$gen]/input$pop,
                  col=rainbow(input$iter)[i], 
                  lwd = lwd.simulated)
            
          } else {
            lines(1:input$gen, 
                  1 - (data()[i,1:input$gen])/input$pop,
                  col=rainbow(input$iter)[i], 
                  lwd = lwd.simulated)
          }
        }
      }

      # Plot expected outcome
      if (input$var.plot == 1) 
          lines(1:input$gen, expected.A()^2, lwd = lwd.expected)
      if (input$var.plot == 2) 
          lines(1:input$gen, 2 * expected.A() * 
                             (1-expected.A()), lwd = lwd.expected)
      if (input$var.plot == 3) 
          lines(1:input$gen, (1-expected.A())^2, lwd = lwd.expected)
      if (input$var.plot == 4) 
          lines(1:input$gen,   expected.A(), lwd = lwd.expected)
      if (input$var.plot == 5) 
          lines(1:input$gen, 1-expected.A(), lwd = lwd.expected)
  })
  output$popPlot <- renderPlot({
     a1a1 <- (expected.A()^2)[input$gen]
     a1a2 <- (2 * expected.A() * (1-expected.A()))[input$gen]
     a2a2 <- ((1-expected.A())^2)[input$gen]
    
    
    plot(0, 0, col = 'white', ylim = c(0, 100), xlim = c(0, 100),
         xaxt="n", yaxt="n", xlab = "", ylab="", cex.lab=1.5, cex.axis=1.3,
         main="Expected Genotypes")
    pheno.col <- c("#440154FF", "#21908CFF", "#FDE725FF")
    text(90,90,label="A1A1", pos=4)
    text(90,80,label="A1A2", pos=4)
    text(90,70,label="A2A2", pos=4)
    flower(x = 90, y = 90, colx=pheno.col[1], cex.x=1, cex.y=1)
    flower(x = 90, y = 80, colx=pheno.col[2], cex.x=1, cex.y=1)
    flower(x = 90, y = 70, colx=pheno.col[3], cex.x=1, cex.y=1)
    x <- sample(1:85, size=100, replace = T)
    y <- sample(1:100, size=100, replace = T)
    
    flower(x,y,colx=sample(pheno.col, size=100, prob=c(a1a1,a1a2,a2a2), replace=T), 
           cex.x=.8, cex.y=1.3)
    

  })
})
