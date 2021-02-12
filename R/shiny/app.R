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
                        max = 500000,
                        value = 300000,
                        step = 5000),
            sliderInput("n1",
                        "Number of sequenced individuals:",
                        min = 100,
                        max = 20000,
                        value = 10000,
                        step = 100),
            sliderInput("gamma1",
                        "Recovery time parameter (gamma):",
                        min = 0,
                        max = 1,
                        value = 1/3.4),
            sliderInput("R1",
                        "Reproduction rate (R):",
                        min = 0,
                        max = 4,
                        step = 0.1,
                        value = 1.2),
            
            sliderInput("NumDays",
                        "Number of Days:",
                        min = 20,
                        max = 100,
                        value = 28),
            sliderInput("IniMean",
                        "Initial number of individuals with VOC:",
                        min = 0,
                        max = 100,
                        value = 11),
            sliderInput("MaxI",
                        "Max number of individuals with VOC:",
                        min = 1,
                        max = 1000,
                        value = 100),
            sliderInput("Threshold",
                        "Threshold:",
                        min = 0,
                        max = 1,
                        value = 0.9),
            radioButtons("IniProb",
                        "IniProb:",
                        choices = c("Kronecker delta" = "1",
                                    #"Uniform" = "2",
                                    "Poisson" = "3"),
                        selected = 3),
            radioButtons("UseMeas",
                         "UseMeas",
                         choices = c("Yes","No"))
            
        ),
        

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot", width = "auto", height = "800px")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    thres <- function(x, threshold){
        match(FALSE,x<threshold)
    }
    
    constructQ <- function(MaxI, gamma, beta){
        Q1 = matrix(0, MaxI, MaxI)
        delta = row(Q1) - col(Q1)
        Q1[delta == 1] <- (1:(MaxI-1))*gamma
        Q1[delta == -1] <- (0:(MaxI-2))*beta
        diag(Q1) <- -1*apply(Q1,1,sum)
        EQ1 = expm(Q1)
        return(EQ1)
        
    }
    
    runCalc <- function(  N1 = 300000, N2=N1,
                          n1 = 10000,  n2=n1,
                          gamma1=1/3.4,gamma2=gamma1,
                          R1=1.2, R2=R1,
                          MaxI=100,
                          NumDays=28,
                          Y = rep(0,NumDays),
                          TimeIntervention = Inf,
                          IniMean=12,
                          UseMeas="Yes",
                          IniProb="3",
                          Threshold = 0.9){
        
        
        # Definitions etc.
        beta1 = R1*gamma1; #Birth rate
        #    n1 = matrix(1, nrow = NumDays, ncol =1)*n1
        
        beta2 = R2*gamma2; #Birth rate
        #    n2 = matrix(1, nrow = NumDays, ncol =1)*n2
        
        if(TimeIntervention < NumDays){    
            n = rbind(matrix(n1, nrow = TimeIntervention-1, ncol=1),
                      matrix(n2, nrow = NumDays-TimeIntervention+1, ncol=1))}else
                          n = matrix(n1, nrow = NumDays, ncol=1)
                      
                      
                      EQ1 = constructQ(MaxI, beta1, gamma1)
                      EQ2 = constructQ(MaxI, beta2, gamma2)
                      
                      # Measurements
                      #    Y = matrix(0, nrow = NumDays, ncol = 1); # Assume no Cluster5
                      #    Y[20,1]= 0*1;                            # Include a measurements of Cluster5
                      #    if (UseMeas == "No"){
                      #      Y= Y*NaN;
                      #    }
                      
                      # Initial probability for states
                      P = matrix(0, nrow = MaxI, ncol = 1);
                      switch(IniProb,                      
                             "1" = {P[IniMean,1] <- 1 },                       # Kronecker delta
                             # "2" = {P[Nlow:Nhigh] <- 1/(Nhigh-Nlow+1) },       # Uniform
                             "3" = {P[,1] <- dpois( 1:(MaxI), IniMean) }       # Poisson
                      )
                      
                      
                      # nstates x ndays container, holding day by day state aposteriori probabilities.
                      PP = matrix(0, nrow = MaxI, ncol = NumDays); 
                      
                      
                      for(i in 1:NumDays){
                          # Measurement update
                          if (!is.nan(Y[i])){ # If there is a measurement
                              P = P*dbinom(Y[i]*as.vector(matrix(1,MaxI,1)),
                                           as.vector(n[i]*matrix(1,MaxI,1)),
                                           (1:(MaxI)-1)/N1);
                              P= P/sum(P);
                          }
                          PP[,i] = P;                # PP does not include the initial PDF
                          # Time update;
                          if(i<TimeIntervention){
                              EQ = EQ1}else{
                                  EQ = EQ2
                              }
                          
                          P = t(EQ)%*%P;
                          
                          # Correction for the error due to truncating number of infected to MaxI
                          P = P/sum(P); 
                      }
                      return(PP)
    }
        
        output$distPlot <- renderPlot({
                PP <- runCalc(N1=input$N1, N2=input$N1,
                              n1=input$n1, n2=input$n1,
                              gamma1=input$gamma1, gamma2=input$gamma1,
                              R1=input$R1, R2=input$R1,
                              MaxI = input$MaxI,
                              NumDays=input$NumDays,
                              Threshold = input$Threshold,
                              IniMean=input$IniMean,
                              UseMeas=input$UseMeas,
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
                    scale_color_viridis_c(name = "Number of days")# +
                    #theme(legend.position = "bottom")
                
                ## Plot distribution function
                cdf_plot <- gg %>%  arrange(NumDays, rowname)  %>% 
                    ggplot(aes(x = rowname, y = cs, group = name, color = NumDays)) +
                    coord_cartesian(xlim=c(0,max_x_value)) +
                    geom_line() +
                    xlab("# Cluster 5") +
                    ylab("CDF") +
                    scale_color_viridis_c(name = "Number of days")# +
                    #theme(legend.position = "bottom")
                
                prob_days_plot <- gg %>% filter(rowname == 1) %>% 
                    ggplot(aes(x = NumDays, y = cs)) +
                    geom_line(color = "#221a52", size = 2) +
                    xlab("Day") + 
                    ylab("Probability of extinction") +
                    ylim(c(0,1)) +
                    geom_hline(yintercept=input$Threshold, linetype="dashed", color = "grey", size = 1) +
                    annotate("text", label = paste(thres(gg$cs, input$Threshold), "days till", input$Threshold, "probability of extinction"),
                             x = 5, y = 1)
                
                ## Combine plots
                #(pdf_plot / cdf_plot) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
                prob_days_plot
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
