## Double sampling function of specific cluster

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


## Kullback-Leibler divergence between a discrete densities dy and dz

kl <- function(dy, dz){
  
  sum(dy[,2]*log(dy[,2]/dz[,2]))
  
}


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
                      Nlow = 1, Nhigh=10,
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
