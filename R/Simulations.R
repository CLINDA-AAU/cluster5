#
# This is an R-script, calculating the results of various intervention strategies, 
# the Danish Mink Case, and Kullback-Leibler Divergence of approximations.
#

library(expm)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(Hmisc)
library(sf)


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

## Simulation of increased WGS testing activity and restrictions

NumDays = 90
nk      = 10000
pk <- c(0.25, 0.5, 0.75, 1.0)
R <- c(0.5, 0.75, 1.0, 1.25, 1.5)
Threshold = c(0.85,0.90,0.95)

for(j in 1:(length(R))){
  for(i in 1:(length(pk))){
  PP <- runCalc(N1 = N,
                n1 = pk[i]*nk,
                gamma1=0.3,
                R1=R[j],
                MaxI=100,
                NumDays=NumDays,
                TimeIntervention=Inf,
                IniMean=12,
                UseMeas="Yes",
                IniProb="3", Threshold = Threshold)
  tmp <- data.frame(Day=1:NumDays, Prob=PP[1,], Ratio = pk[i], R = R[j])
  if(i==1 && j==1){
    res <- tmp}else 
      res <- rbind(res, tmp)
  }
}

res.testing <- res[res$R==1,] 

res.testing$Ratio <- as.factor(as.character(res.testing$Ratio))
  
p.testing <- ggplot(res.testing, aes(x=Day, y=Prob)) + 
  geom_point(res.testing, mapping = aes(x=Day, y=Prob, color=Ratio)) + 
  xlab("Day") + ylab("Probability of extinction") + 
  geom_hline(yintercept=0.85, linetype="dashed", color = "grey", size = 1) + 
  geom_hline(yintercept=0.9, linetype="dashed", color = "grey", size = 1) +
  geom_hline(yintercept=0.95, linetype="dashed", color = "grey", size = 1) +
  scale_color_brewer(palette="YlOrRd")
p.testing$labels$colour <- "WGS Ratio"

res.restrictions <- res[res$Ratio==0.25,] 

res.restrictions$R <- as.factor(as.character(res.restrictions$R))

p.restrictions <- ggplot(res.restrictions, aes(x=Day, y=Prob)) + 
  geom_point(res.restrictions, mapping = aes(x=Day, y=Prob, color=R)) + 
  xlab("Day") + ylab("Probability of extinction") + 
  geom_hline(yintercept=0.85, linetype="dashed", color = "grey", size = 1) + 
  geom_hline(yintercept=0.9, linetype="dashed", color = "grey", size = 1) +
  geom_hline(yintercept=0.95, linetype="dashed", color = "grey", size = 1) +
  scale_color_brewer(palette="YlOrRd")
p.restrictions$labels$colour <- "R"

# Figure

p.testing / p.restrictions + plot_annotation(tag_levels = "A") 
ggsave("Interventions.pdf", width = 18, units="cm")

# Table

thresholds <- res %>% 
  group_by(Ratio, R) %>%
  summarise("Prob < 0.85" = thres(Prob,0.85),
            "Prob < 0.90" = thres(Prob,0.90),
            "Prob < 0.95" = thres(Prob,0.95))

rownames(thresholds) <- NULL

latex(thresholds,file="days.tex", 
      caption = "Days to thresholds for various testing and restriction strategies.", 
      label = "tab:intervention", rowname=NULL)


## Danish mink case

# Map of Denmark

geo_sf <- read_sf("https://dawa.aws.dk/kommuner/?format=geojson")
geo_sf <- geo_sf %>% 
  mutate(grp=case_when(navn %in% c("Hjørring","Frederikshavn") ~ "Cluster 5",
                       navn %in% c("Morsø","Thisted","Brønderslev","Frederikshavn",
                                   "Vesthimmerlands","Læsø","Rebild","Mariagerfjord","Jammerbugt","Aalborg") ~ "NDR",
                       TRUE ~ "Other"))

geo_sf %>% ggplot() +
  geom_sf(aes(fill = grp), size=0.1) +
  scale_fill_manual(values=c("red","blue","white"))+
  labs(title="", fill="Category") +
  theme_void()+
  theme(legend.position=c(0.8,0.8))
ggsave("Denmark.pdf",width = 18, units="cm")

# Generate Danish data

# Read Danish Cluster 5 data
c5 <- read.csv("Data/rn_cluster5_data.csv", sep = ";")
c5[is.na(c5)] <- 0
c5 <- c5[c5$week >= "2020-w35" & c5$week<="2020-w49",]

N   <- 589148                 # Population size (North Denmark Region)
nk <- median(c5$sequenced)

PP <- runCalc(N1 = N,
        n1 = nk,
        gamma1=0.3,
        R1=R[j],
        MaxI=100,
        NumDays=dim(c5)[1],
        Y = c5$cluster5,
        TimeIntervention=Inf,
        IniMean=4,
        UseMeas="Yes",
        IniProb="3", Threshold = Threshold)

plot(PP[1,])

# Intervention analysis

N   <- 589148                 # Population size (North Denmark Region)
xkmm <- 0.006*N               # 0.6% of the population is positive (median)
xk  <- round(0.006*N*0.08)    # 8% of infected with specific variant
xkm <- xkmm - xk              # Number of infected with non-specific variant
pk <- 0.28                    # Ratio of positives WGS tested

nk  <- 17000   # Numer of PCR tests (median)
mk  <- round(pk*nk*(xk+xkm)/N)  # 0.28 of positive sent for WGS tests (median)

NumDays = 70

