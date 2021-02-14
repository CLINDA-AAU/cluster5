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
dyk$Type = "Sample"

# Poisson approximation

lam <- mk*xk/(xk+xkm)
dzkpois <- data.frame(yk=dyk$yk, Freq=dpois(dyk[,1],lam), Type = "Poisson")
KLpois <- kl(dyk, dzk)

# Binomial approximation

dzk.binom1 <- data.frame(yk=dyk$yk, Freq =dbinom(dyk[,1], mk, xk/(xk+xkm)), Type = "Binom1")
KLbinom1 <- kl(dyk, dzk.binom1)

# Binomial approximation

dzk.binom2 <- data.frame(yk=dyk$yk, Freq =dbinom(dyk[,1], pk*nk, xk/N), Type = "Binom2")
KLbinom2 <- kl(dyk, dzk.binom2)

dyk.all <- rbind(dyk, dzkpois, dzk.binom1, dzk.binom2)

# Kullback-Leibler divergence between sampling distribution and plot

p.distr <- ggplot(dyk.all, aes(x=yk,y=Freq)) +
  geom_point(dyk.all, mapping = aes(x=yk,y=Freq,color=Type), size = 2) +
  scale_color_manual(values=AAU_pal("second")(6)) +
  xlab("Number") + ylab("Probability")
p.distr
ggsave("Distribution.pdf", width = 18, units="cm")

kld <- data.frame(Type = c("Poisson", "Binom1", "Binom2"),
                 KLD  = round(c(KLpois, KLbinom1, KLbinom2),4))

latex(kld,file="kld.tex", 
      caption = "Kulback-Leibler divergence between the approximations and the simulated sampling distribution",
      label = "tab:kld", rowname=NULL)


