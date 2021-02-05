library(tidyverse)

# Double sampling of specific cluster

cluster.sampling <- function(N, xk, xkm, nk, mk, nsim = 100){

  yk <- rep(0,nsim)
  
  for(i in 1:nsim){
    x <- rep("Neg",N)
    x[1:xk] <- "C"
    x[(xk+1):(xk+xkm)] <- "Pos"
    
    z <- sample(x, nk)
    
    y <- sample(z[z!="Neg"], mk)
    
    yk[i] <- sum(y=="C")
    
  }
  yk
}


# Kullback-Leibler divergence between a discrete density dy 
# and the Poisson distribution with parameter lamda

kl <- function(dy, lambda){
  
  sum(dy[,2]*log(dy[,2]/dpois(dy[,1],lambda)))

}


# Read sequence data from the Danish Cluster 5 case

seqdata <- read.csv("Data/rn_cluster5_data.csv", sep = ";")

N   <- 600000                 # Population size (North Denmark Region)
xkmm <- 0.006*N               # 0.6% of the population is positive (median)
xk  <- round(0.006*N*0.08)    # 8% of infected with specific variant
xkm <- xkmm - xk              # Number of infected with non-specific variant

nk  <- 17000   # Numer of PCR tests (median)
mk  <- round(0.28*nk*(xk+xkm)/N)  # 0.28 of positive sent for WGS tests (median)

# Simulation of the sampling distribution

yk <- cluster.sampling(N, xk, xkm, nk, mk, nsim = 10000)

# Poisson approximation

lam <- mk*xk/(xk+xkm)

# Kullback-Leibler divergence between sampling distribution and plot

KL <- kl(dyk, lam)
plot(dyk[,1],dyk[,2], xlab = "Number", ylab = "Probability")
points(dyk[,1],dpois(dyk[,1],lam), col = "red")
legend("topright", pch = c(1,1), col = c("black", "red"), 
       legend = c("Sample","Poisson"), bty = "n")
title(paste("KLD = ", as.character(round(KL,5))))
