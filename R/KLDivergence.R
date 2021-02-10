library(tidyverse)

# Double sampling function of specific cluster

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


# Kullback-Leibler divergence between a discrete densities dy 
# and dz

kl <- function(dy, dz){
  
  sum(dy[,2]*log(dy[,2]/dz[,2]))

}

# Read sequence data from the Danish Cluster 5 case

seqdata <- read.csv("Data/rn_cluster5_data.csv", sep = ";")

N   <- 589148                 # Population size (North Denmark Region)
xkmm <- 0.006*N               # 0.6% of the population is positive (median)
xk  <- round(0.006*N*0.08)    # 8% of infected with specific variant
xkm <- xkmm - xk              # Number of infected with non-specific variant
pk <- 0.28                    # Ratio of positives WGS tested

nk  <- 17000   # Numer of PCR tests (median)
mk  <- round(pk*nk*(xk+xkm)/N)  # 0.28 of positive sent for WGS tests (median)

# Simulation of the sampling distribution

nsim <- 10000
yk <- cluster.sampling(N, xk, xkm, nk, mk, nsim = nsim)
dyk <- as.data.frame(table(yk))
dyk$yk <- as.integer(levels(dyk$yk))
dyk$Freq <- dyk$Freq/nsim

# Poisson approximation

lam <- mk*xk/(xk+xkm)
dzkpois <- data.frame(yk=dyk$yk, Freq=dpois(dyk[,1],lam))
KLpois <- kl(dyk, dzk)

# Binomial approximation

dzk <- data.frame(yk=dyk$yk, Freq =dbinom(dyk[,1], mk, xk/(xk+xkm)))
KLbinom1 <- kl(dyk, dzk)

# Binomial approximation

dzk <- data.frame(yk=dyk$yk, Freq =dbinom(dyk[,1], pk*nk, xk/N))
KLbinom2 <- kl(dyk, dzk)

# Kullback-Leibler divergence between sampling distribution and plot

plot(dyk[,1],dyk[,2], xlab = "Number", ylab = "Probability")
points(dyk[,1],dpois(dyk[,1],lam), col = "red")
points(dyk[,1],dbinom(dyk[,1], mk, xk/(xk+xkm)), col = "blue")
points(dyk[,1],dbinom(dyk[,1], pk*nk, xk/N), col = "green")
legend("topright", pch = c(1,1), col = c("black", "red","blue", "green"), 
       legend = c("Sample","Poisson", "Binom1", "Binom2"), bty = "n")
title(paste("KLD = ", as.character(round(KLpois,5)),",",as.character(round(KLbinom1,5)),"",as.character(round(KLbinom2,5))))


