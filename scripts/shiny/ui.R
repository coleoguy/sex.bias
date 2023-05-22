#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
fluidPage(
  titlePanel("OSR Fitness Differences"),
  column(4,
         
         fluidRow(
           column(6,
                  
                  h3("Plot Basics"),
                  selectInput("chromosome", label = h3("Chromosome selection"), 
                              choices = list("Auto" = "0.5", "Sex" = "0.2"), 
                              selected = "0.5"),
                  selectInput("sex", label = h3("Common Sex"),
                              choices = list("Male" = "1", "Female" = "2"),
                              selected = "Male"),
                  selectInput("pop", label = h3("Common Population size"), 
                              choices = list("50" = "50", "100" = "100", "500" = "500", "1000" = "1000"), 
                              selected = "1000"),

                  
                  h3("Selection"),
                  selectInput("h", label = h3("Dominance factor"), 
                              choices = list("0" = "0", "0.5" = "0.5", "1" = "1"), 
                              selected = "0.5"),
                  
                  selectInput("s", label = h3("Selection Pressure"), 
                              choices = list("0.1" = "0.1", "0.5" = "0.5", "0.9" = "0.9"), 
                              selected = "0.1"),
                  
                  selectInput("osr", label = h3("Operational Sex Ratio"), 
                              choices = list("1" = "1", ".8" = ".8", ".6" = ".6", ".4" = ".4", ".2" = ".2", ".1" = ".1", "0.5" = "0.5"), 
                              selected = "1"),
                  
           ),
           
           column(6,
                  
                  h3("Drift"),
                  sliderInput("iter", "Replicates:", 
                              min = 0, max = 50, value = 0, step = 1),
                  selectInput("var.plot", "What to plot:",
                              list("Wdiff" = 1, 
                                   "XFreq" = 3, 
                                   "Yfreq" = 4,
                                   "AFreq"  = 2)),
                  actionButton("seed.val", 'Refresh'),

           )
         )
  ),
  
  column(8,
         plotOutput("genePlot", width='100%', height='600px')
  )
)

