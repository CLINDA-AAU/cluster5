kl <- function(dy, lam){
  res <- 0
  for(i in 1:length(dy[1,])){
    res <- res + dy[2,i] * log(dy[2,i]/dpois(dy[1,i],lam))
  }
  res
}

N   <- 5000000                # Population size
xkmm <- 0.015*N               # 1.5% of the population is positive
xk  <- round(0.015*0.003*N)   # 0.3% of infected with specific variant
xkm <- xkmm - xk              # Number of infected with non-specific variant

nk  <- 350000   # Numer of PCR tests
mk  <- round(0.2*nk*(xk+xkm)/N)  # 0.2 of positive sent for WGS tests

# Simulation of the sampling distribution

nsim <- 1000
yk <- rep(0,nsim)

for(i in 1:nsim){
print(i)
    x <- rep("Neg",N)
  x[1:xk] <- "C"
  x[(xk+1):(xk+xkm)] <- "Pos"

  z <- sample(x, nk)

  y <- sample(z[z!="Neg"], mk)

  yk[i] <- sum(y=="C")

}

dyk <- t(as.matrix(table(yk)))/nsim
dyk <- rbind(as.integer(colnames(dyk)), as.vector(dyk))


# Poisson approximation

lam <- mk*(xk)/(xk+xkm)

# Kullback-Leibler divergence

KL <- kl(dyk, lam)
plot(dyk[1,],dyk[2,], xlab = "Number", ylab = "Probability")
points(dyk[1,],dpois(dyk[1,],lam), col = "red")
legend("topright", pch = c(1,1), col = c("black", "red"), 
       legend = c("Sample","Poisson"), bty = "n")
text(8,0.2,as.character(round(KL,5)))
