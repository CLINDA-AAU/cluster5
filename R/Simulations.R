#
# This is an R-script, calculating the results of various intervention strategies, 
# the Danish Mink Case, and Kullback-Leibler Divergence of approximations.
#

# Libraries

library(expm)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(Hmisc)
library(sf)
library(scales)

# Functions

## AAU colours

AAU_palettes <- list(
  "main" = c("#54616E","#594FBF","#211A52"),
  "second" = c("#bb5b17", "#97701f", "#007fa3", "#a16547", "#0e8563", "#cc445b"),
  "third"  = c("#df8e2e", "#b19336", "#31a9c1", "#cc8b88", "#5caf8d", "#e78293")
)

AAU_pal <- function(palette = "main", reverse = FALSE, ...) {
  pal <- AAU_palettes[[palette]]
  
  if (reverse) pal <- rev(pal)
  
  colorRampPalette(pal, ...)
}

### Testing the colour schemes

x = 1:100
data = data.frame("a" = x,
                  "b" = x*1.25,
                  "c" = x*1.5,
                  "d" = x*1.75,
                  "e" = x*2,
                  "f" = x*3,
                  "g" = x*4) %>% 
  pivot_longer(cols=b:g)

# Gradient over primary color scheme
ggplot(data, aes(x = a, y = value, fill = name)) +
  geom_boxplot() +
  scale_fill_manual(values = AAU_pal("main")(6))

# Gradient for secondary color scheme
ggplot(data, aes(x = a, y = value, fill = name)) +
  geom_boxplot() +
  scale_fill_manual(values = AAU_pal("second")(6))

# Gradient for tertiary color scheme
ggplot(data, aes(x = a, y = value, fill = name)) +
  geom_boxplot() +
  scale_fill_manual(values = AAU_pal("third")(6))


## Days to threhold for extinction

thres <- function(x, threshold){
  match(FALSE,x<threshold)
}


## Construction of rate matrices
 
constructQ <- function(MaxI, gamma, beta){
  Q1 = matrix(0, MaxI, MaxI)
  delta = row(Q1) - col(Q1)
  Q1[delta == 1] <- (1:(MaxI-1))*gamma
  Q1[delta == -1] <- (0:(MaxI-2))*beta
  diag(Q1) <- -1*apply(Q1,1,sum)
  EQ1 = expm(Q1)
  return(EQ1)
  
}


## Calculation of a posteriori probabilities

runCalc <- function(  N1 = 300000, N2=N1,
                      n1 = 10000,  n2=n1,
                      gamma1=1/3.4, gamma2=gamma1,
                      R1=1.2, R2=R1,
                      MaxI=100,
                      TimeIntervention = Inf,
                      NumDays=28,
                      Y = rep(0,NumDays),
                      IniMean=12,
                      IniProb="Kronecker"
                  ){
    # Set the birth rate corresponding to the reproduction rate and recovery time    
    beta1 = R1*gamma1; # Before intervention
    beta2 = R2*gamma2; # After intervention

    if(TimeIntervention < NumDays){    
    n = rbind(matrix(n1, nrow = TimeIntervention-1, ncol=1),
              matrix(n2, nrow = NumDays-TimeIntervention+1, ncol=1))}else
    n = matrix(n1, nrow = NumDays, ncol=1)

    # Transition probability matrix
        
    EQ1 = constructQ(MaxI+1, beta1, gamma1)
    EQ2 = constructQ(MaxI+1, beta2, gamma2)
        
    # Initial probability for states
    P = matrix(0, nrow = MaxI+1, ncol = 1);
    switch(IniProb,                      
           "Kronecker" = {P[IniMean+1,1] <- 1 },                # Kronecker delta
           "Uniform" = {P[Nlow:Nhigh,1] <- 1/(Nhigh-Nlow+1) },  # Uniform
           "Poisson" = {P[,1] <- dpois(0:MaxI, IniMean) }       # Poisson
    )
    
    # nstates x ndays container, holding day by day state a posterior probabilities.
    PP = matrix(0, nrow = MaxI+1, ncol = NumDays); 
    
    for(i in 1:NumDays){
      # Measurement update
      if (!is.nan(Y[i])){
        P = P*dbinom(Y[i]*as.vector(matrix(1,MaxI+1,1)),
                     as.vector(n[i]*matrix(1,MaxI+1,1)),
                     (1:(MaxI+1)-1)/N1)
        P= P/sum(P)
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

# Basic statistics

north.denmark <- data.frame(Municipality = c("Frederikshavn", "Hjørring", "Jammerbugt", 
                                             "Brønderslev", "Thisted", "Vesthimmerland",
                                             "Læsø", "Morsø", "Rebild", "Aalborg"),
                            Inhabitants  = c(61158, 66178, 38611, 
                                             35754, 44908, 37534,
                                             1897, 21474, 28911, 201142),
                            Type = c(rep("Locked",7),rep("non-Locked",3)))

nd.summary <- north.denmark %>% 
  group_by(Type) %>% summarise(Inhabitants = sum(Inhabitants))

N.locked <- as.integer(nd.summary[1,2])
N <- sum(nd.summary["Inhabitants"])

# Results of increased WGS testing activity and restrictions

NumDays = 60
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
                IniProb="Kronecker")
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
  xlab("Week") + ylab("Probability of extinction") + 
  geom_hline(yintercept=0.85, linetype="dashed", color = "grey", size = 1) + 
  geom_hline(yintercept=0.9, linetype="dashed", color = "grey", size = 1) +
  geom_hline(yintercept=0.95, linetype="dashed", color = "grey", size = 1) +
#  scale_color_brewer(palette="YlOrRd")
  scale_color_manual(values=AAU_pal("second")(6))
p.testing$labels$colour <- "WGS Ratio"

res.restrictions <- res[res$Ratio==0.25,] 

res.restrictions$R <- as.factor(as.character(res.restrictions$R))

p.restrictions <- ggplot(res.restrictions, aes(x=Day, y=Prob)) + 
  geom_point(res.restrictions, mapping = aes(x=Day, y=Prob, color=R)) + 
  xlab("Week") + ylab("Probability of extinction") + 
  geom_hline(yintercept=0.85, linetype="dashed", color = "grey", size = 1) + 
  geom_hline(yintercept=0.9, linetype="dashed", color = "grey", size = 1) +
  geom_hline(yintercept=0.95, linetype="dashed", color = "grey", size = 1) +
  scale_color_manual(values=AAU_pal("second")(6))
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
      caption = "Weeks to thresholds for various testing and restriction strategies.", 
      label = "tab:intervention", rowname=NULL)



# Effect of reproduction rate

## Simulation, effect of reproduction rate

NumDays = 60
nk      = 10000
pk <- 0.25
Threshold = c(0.1,0.90,0.95)
R <- seq(0.5, 2.5, by = 0.01)
gamma1 <- 1

for(i in 1:length(R)){
  PP <- runCalc(N1 = N,
                n1 = pk*nk,
                gamma1=gamma1,
                R1=R[i],
                MaxI=100,
                NumDays=NumDays,
                TimeIntervention=Inf,
                IniMean=12,
                IniProb="Kronecker")
  tmp <- data.frame(Day=1:NumDays, Prob=PP[1,], Ratio = pk, R = R[i])
  if(i==1){
    res <- tmp}else 
      res <- rbind(res, tmp)
}

thresholds <- res %>% 
  group_by(R) %>%
  summarise("Prob < 0.85" = thres(Prob,0.85),
            "Prob < 0.90" = thres(Prob,0.90),
            "Prob < 0.95" = thres(Prob,0.95))

thresholds.long = thresholds %>%
  gather(`Prob < 0.85`, `Prob < 0.90`, `Prob < 0.95`, key = "Prob", value = "Week")
thresholds.long

p.reproduction <- ggplot(thresholds.long,aes(x=R, y=Prob)) + 
  geom_line(thresholds.long, mapping = aes(x=R, y=Week, color = Prob)) +
  xlab("R") + ylab("Week") + 
  scale_color_manual(values=AAU_pal("second")(6))
p.reproduction
ggsave("Reproduction.pdf", width = 18, units = "cm")

#
## Danish mink case
#

# Map of Denmark

geo_sf <- read_sf("https://dawa.aws.dk/kommuner/?format=geojson")
geo_sf <- geo_sf %>% 
  mutate(grp=case_when(navn %in% c("Hjørring","Frederikshavn") ~ "Cluster 5",
                       navn %in% c("Thisted","Brønderslev","Frederikshavn",
                                   "Vesthimmerlands","Læsø","Jammerbugt") ~ "Locked down",
                       navn %in% c("Morsø","Thisted","Brønderslev","Frederikshavn",
                                   "Vesthimmerlands","Læsø","Rebild","Mariagerfjord","Jammerbugt","Aalborg") ~ "North Denmark",
                       TRUE ~ "Remaining"))

geo_sf %>% ggplot() +
  geom_sf(aes(fill = grp), size=0.1) +
  scale_fill_manual(values=c(rev(c(AAU_pal("second")(3)[1:3])),"white"))+
  labs(title="", fill="Status week 46, 2020") +
  theme_void()+
  theme(legend.position=c(0.8,0.8))
ggsave("Denmark.pdf",width = 18, units="cm")

# Generate Danish mink case

# Read Danish Cluster 5 data
c5 <- read.csv("Data/rn_cluster5_data.csv", sep = ";")
c5$Week <- as.integer(substring(c5$week,7))
c5[is.na(c5)] <- 0
c5 <- c5[c5$Week >= 35 & c5$Week<=53,]

# Intervention table

intervention <- data.frame(Week = c(46,47,50),
                           Announced = as.Date(c("11/06/20","11/13/20","11/06/20"), "%m/%d/%y"),
                           Effect = as.Date(c("11/09/20","11/16/20","12/07/20"), "%m/%d/%y"),
                           Description = c("Lockdown", "Reopening", "Planned"))

dates <- c("02/27/92", "02/27/92", "01/14/92", "02/28/92", "02/01/92")
as.Date(dates, "%m/%d/%y")

# Analyse Danish mink case

R1 <- 1.2                   # Reproduction rate before intervention
R2 <- 1                   # Reproduction rate after intervention
gamma1 <- 1.0               # Average contamination time 2 weeks

# It was the plan to test the entire population of the region over a 4 weeks period
# i.e., nk = N.locked/4 = 71510 PCR tests per week. Assume a postive pct. of 1.5%, i.e,
# 1072 positive tests. The test capacity was up to 5000 a week, so we assume all 1072
# will be WGS sequence, i.e. n2=1072

# Planned
# If restrictions and testing have been realized

PP <- runCalc(N1 = N.locked,
        n1 = round(c5$cases[1:11]/2),
        n2 =  1072,
        gamma1=gamma1,
        R1=R1,
        R2=R2,
        MaxI=100,
        NumDays=dim(c5)[1],
        Y = c5$cluster5,
        TimeIntervention=12,
        IniMean=4,
        IniProb="Kronecker")

res.mink <- data.frame(Week=c5$Week, Prob=PP[1,], 
                       Strategy="Planned")

# Realized
# Actual data

PP <- runCalc(N1 = N.locked,
              n1 = round(c5$cases/2),
              gamma1=gamma1,
              R1=R1,
              R2=R2,
              MaxI=100,
              NumDays=dim(c5)[1],
              Y = c5$cluster5,
              TimeIntervention=Inf,
              IniMean=4,
              IniProb="Kronecker")

res.mink <- rbind(res.mink, data.frame(Week=c5$Week, Prob=PP[1,], 
                          Strategy="Realized"))

res.mink$Strategy <- factor(res.mink$Strategy)

p.minkcase <- ggplot(res.mink, aes(x=Week, y=Prob)) + 
  geom_point(res.mink, mapping = aes(x=Week, y=Prob, color=Strategy), size=2)+
  geom_vline(xintercept=45, linetype="dashed", color = "grey", size = 1) +
  geom_vline(xintercept=47, linetype="dashed", color = "grey", size = 1) + 
  geom_vline(xintercept=50, linetype="dashed", color = "grey", size = 1) + 
  xlab("Week") + ylab("Probability of extinction") +
  scale_color_manual(values=AAU_pal("second")(6))
p.minkcase
ggsave("Mink.pdf", width = 18, units="cm")

# Table of Cuts

cut.week <- function(x,week){
  x[x$Week == week,]
}

cut.mink <- rbind(cut.week(res.mink,45),
      cut.week(res.mink,47),
      cut.week(res.mink,50))

cut.mink$Prob <- round(cut.mink$Prob,4)

cut.mink.wide <- cut.mink %>% spread(Strategy, Prob)

latex(cut.mink.wide,file="cut.tex", 
      caption = "Planned and realized probability of extinction, before, at the break, and
      at time of planned intervention",  label = "tab:min", rowname=NULL)

