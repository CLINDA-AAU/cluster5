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
    titlePanel("Cluster5"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("N",
                        "N:",
                        min = 100000,
                        max = 500000,
                        value = 300000),
            sliderInput("n",
                        "n:",
                        min = 5000,
                        max = 20000,
                        value = 10000),
            sliderInput("gamma",
                        "gamma:",
                        min = 0,
                        max = 1,
                        value = 1/3.4),
            sliderInput("R",
                        "Reproduction rate (R):",
                        min = 0,
                        max = 10,
                        step = 0.1,
                        value = 1.2),
            sliderInput("MaxI",
                        "MaxI:",
                        min = 1,
                        max = 1000,
                        value = 100),
            sliderInput("NumDays",
                        "Number of Days:",
                        min = 20,
                        max = 100,
                        value = 28),
            sliderInput("IniMean",
                        "IniMean:",
                        min = 0,
                        max = 100,
                        value = 12),
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
    
    runCalc <- function(N=300000,
                        n=10000,
                        gamma=1/3.4,
                        R=1.2,
                        MaxI=100,
                        NumDays=28,
                        IniMean=12,
                        UseMeas="Yes",
                        IniProb="3"){
        
        
        # Definitions etc.
        beta = R*gamma; #Birth rate
        TimeIntervention = Inf
        n = matrix(1, nrow = NumDays, ncol =1)*n
        
        # Construct Q
        Q1 = matrix(0, MaxI, MaxI)
        delta = row(Q1) - col(Q1)
        Q1[delta == 1] <- (1:(MaxI-1))*gamma
        Q1[delta == -1] <- (0:(MaxI-2))*beta
        diag(Q1) <- -1*apply(Q1,1,sum)
        EQ1 = expm(Q1)
        
        # Measurements
        Y = matrix(0, nrow = NumDays, ncol = 1); # Assume no Cluster5
        Y[20,1]= 0*1;                            # Include a measurements of Cluster5
        if (UseMeas == "No"){
            Y= Y*NaN;
        }
        
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
                             (1:(MaxI)-1)/N);
                P= P/sum(P);
            }
            PP[,i] = P;                # PP does not include the initial PDF
            # Time update;
            if(i<TimeIntervention){
                EQ = EQ1}else{
                    EQ = EQ2}
            
            P= t(EQ)%*%P;
            
            # Correction for the error due to truncating number of infected to MaxI
            P = P/sum(P); 
        }
        return(PP)
    }
        
        output$distPlot <- renderPlot({
                PP <- runCalc(N=input$N,
                              n=input$n,
                              gamma=input$gamma,
                              R=input$R,
                              MaxI = input$MaxI,
                              NumDays=input$NumDays,
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
                    geom_line() +
                    ylab("Probability of no Cluster 5") +
                    xlab("Days since last observation") +
                    ylim(c(0,1))
                
                ## Combine plots
                #(pdf_plot / cdf_plot) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
                prob_days_plot
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
