#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(expm)
library(shiny)
library(ggplot2)
library(tidyverse)
library(patchwork)

source("functions.R")

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    
    # Application title
    titlePanel( div(column(width = 6, h2("Covid19 Variant of Concern Monitor")),
                    column(width = 6, tags$img(src = "__AAU_LEFT_RGB_UK.png", align = "right"))), #, width = "25%", height = "25%"))),
                windowTitle="Covid19 Variant of Concern Monitor"
    ),
    # tags$img(src = "__AAU_LEFT_RGB_UK.png"),
    # titlePanel("Covid19 Variant of Concern Monitor"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("N1",
                        "Population size:",
                        min = 100000,
                        max = 1000000,
                        value = 600000,
                        step = 5000),
            sliderInput("n1",
                        "Number of sequenced individuals per week:",
                        min = 100,
                        max = 20000,
                        value = 10000,
                        step = 100),
            sliderInput("gamma1",
                        "Recovery time parameter (gamma):",
                        min = 0,
                        max = 1,
                        value = 0.5),
            sliderInput("R1",
                        "Reproduction rate (R):",
                        min = 0,
                        max = 4,
                        step = 0.1,
                        value = 1.0),
            
            sliderInput("NumDays",
                        "Number of weeks to run calculation:",
                        min = 1,
                        max = 100,
                        value = 60),
            #sliderInput("MaxI",
            #            "Max number of individuals with VOC:",
            #            min = 1,
            #            max = 1000,
            #            value = 100),
            sliderInput("Threshold",
                        "Desired probability of extinction:",
                        min = 0,
                        max = 1,
                        value = 0.9),
            selectInput("IniProb", "Probability distribution for initial number of VOC:", choices = c("Atom", "Poisson", "Uniform")),
            conditionalPanel(
                "input.IniProb == 'Uniform'",
                sliderInput("Nuni", label = "Uniform Range:", min = 1, 
                            max = 100, value = c(1, 10))
            ),
            conditionalPanel(
                "input.IniProb == 'Poisson'",
                sliderInput("IniMean",
                            "Poisson mean:",
                            min = 0,
                            max = 100,
                            value = 11)
            ),
            conditionalPanel(
                "input.IniProb == 'Atom'",
                sliderInput("IniMean2",
                            "Atom in:",
                            min = 0,
                            max = 100,
                            value = 11)
            )
            
        ),
        

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot", width = "auto", height = "800px")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    output$distPlot <- renderPlot({
                IniMean <- ifelse(input$IniProb == "Atom", input$IniMean2, input$IniMean)
                
                PP <- runCalc(N1=input$N1, N2=input$N1,
                              n1=input$n1, n2=input$n1,
                              gamma1=input$gamma1, gamma2=input$gamma1,
                              R1=input$R1, R2=input$R1,
                              Nlow = input$Nuni[1], Nhigh = input$Nuni[2],
                              MaxI = NA,
                              Y = c(11,rep(0,input$NumDays-1)),
                              NumDays=input$NumDays,
                              IniMean= IniMean,
                              IniProb=input$IniProb)
                
                gg <- data.frame(ts(PP)) %>% 
                    rownames_to_column() %>% 
                    mutate(rowname = as.numeric(rowname)) %>% 
                    pivot_longer(!rowname) %>% 
                    mutate("NumDays" = as.numeric(sub("Series\\.", "", name))) %>% 
                    group_by(NumDays) %>% 
                    mutate(cs = cumsum(value))
                
                ## Pick a nice looking xlim
                max_x_value <- gg %>% group_by(rowname) %>% 
                    mutate(pdf_sum = sum(value)) %>% 
                    filter(pdf_sum < 0.01) %>% 
                    pull(rowname) %>% 
                    min()
                
                
                ## Plot probability density function
                pdf_plot <- gg %>%
                    ggplot(aes(x = rowname, y = value, group = name, color = NumDays)) +
                    coord_cartesian(xlim=c(0,max_x_value)) +
                    geom_line() +
                    xlab("# Cluster 5") +
                    ylab("PDF") +
                    scale_color_viridis_c(name = "Number of weeks")# +
                    #theme(legend.position = "bottom")
                
                ## Plot distribution function
                cdf_plot <- gg %>%  arrange(NumDays, rowname)  %>% 
                    ggplot(aes(x = rowname, y = cs, group = name, color = NumDays)) +
                    coord_cartesian(xlim=c(0,max_x_value)) +
                    geom_line() +
                    xlab("# Cluster 5") +
                    ylab("CDF") +
                    scale_color_viridis_c(name = "Number of weeks")# +
                    #theme(legend.position = "bottom")
                
                prob_days_plot <- gg %>% filter(rowname == 1) %>% 
                    ggplot(aes(x = NumDays, y = cs)) +
                    geom_line(color = "#211a52", size = 2) +
                    xlab("Day") + 
                    ylab("Probability of extinction") +
                    ylim(c(0,1)) +
                    geom_hline(yintercept=input$Threshold, linetype="dashed", color = "grey", size = 1) +
                    annotate("text", label = paste(thres(gg$cs, input$Threshold), "weeks till", input$Threshold, "probability of extinction"),
                             x = 5, y = 1)
                
                ## Combine plots
                #(pdf_plot / cdf_plot) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
                prob_days_plot
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
