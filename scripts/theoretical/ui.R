# This interface is modified from the evobiR package
# http://coleoguy.github.io/software.html

shinyUI(fluidPage(

  #titlePanel("", "Change in Allele or Genotype Frequency"),

  fluidRow(
    column(4,

      fluidRow(
        column(6,

          h3("Plot Basics"),
          sliderInput("initial.A", "Initial frequency of A1:", 
                      min = 0, max = 1, value=0.5, step =0.01),
          sliderInput("gen", "Generations:", 
                      min = 10, max = 1000, value = 100, step = 10),
          selectInput("var.plot", "What to plot:",
                      list("A1" = 4, 
                           "A2" = 5, 
                           "A1 A1" = 1,
                           "A1 A2"  = 2,
                           "A2 A2"  = 3)),
         # TODO: checkboxes for what to plot, so can see multiple quantities at once

          h3("Selection"),
          sliderInput("fit.AA", "Fitness of A1 A1:", 
                      min = 0, max = 1.5, value = 1, step = 0.05),
          sliderInput("fit.Aa", "Fitness of A1 A2:", 
                      min = 0, max = 1.5, value = 1, step = 0.05),
          sliderInput("fit.aa", "Fitness of A2 A2:", 
                      min = 0, max = 1.5, value = 1, step = 0.05)

        ),

        column(6,

          h3("Drift"),
          sliderInput("pop", "Population Size:", 
                      min = 10, max = 1000, value = 100, step = 10),
          sliderInput("iter", "Replicates:", 
                      min = 0, max = 50, value = 0, step = 1),
          # TODO: option for starting with a single copy of A1
          h3("Mutation"),
          sliderInput("qAa", "Rate from A1 to A2:", 
                      min = 0, max = 0.1, value = 0, step = 0.0001),
          sliderInput("qaA", "Rate from A2 to A1:", 
                      min = 0, max = 0.1, value = 0, step = 0.0001),
          h5(strong("Re-run simulations:")),
          actionButton("seed.val", 'Refresh')
        )
      )
    ),

    column(8,
           plotOutput("genePlot", width='100%', height='600px'),
           plotOutput("popPlot", width='400px', height='300px')
    )
  )
))
