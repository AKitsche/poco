#Function to calculate the probability of observing a consistent result
#probability delta_{J}/delta_{ALL} > omega and delta{ALL} > 0
#this is denoted by eta_{B} in:

#Sample size and proportion of Japanese patients in multi-regional trials
#Kimitoshi Ikeda and Frank Bretz
#Pharmaceutical Statistics 9: 207–216 (2010)
library(mvtnorm)
powB <- function(delta=0.125, omega=0.5 , sigma=1 ,alpha=0.025, beta=0.2, p=0.1){
  Zalpha <- qnorm(p=1-alpha)
  Zbeta  <- qnorm(p=1-beta)
  na <- nb <- as.integer((2*(Zalpha+Zbeta)^2*sigma^2)/(delta^2),0)
  #x      <- (delta*(1-omega)*sqrt(na*nb*p))/(sigma*sqrt((na+nb)*(p*omega^2-2*p*omega+1)))
  #etaA   <- dnorm(x=x)
  l1      <- -(delta/sigma)*sqrt((na*nb)/(na+nb))
  l2      <- -(delta*(1-omega)*sqrt(na*nb*p))/(sigma*sqrt(na+nb)*(p*omega^2-2*p*omega+1))
  R       <- diag(c(1,1))
  R[1,2] <-  R[2,1]  <- (sqrt(p)*(1-omega))/(sqrt(p*omega^2-2*p*omega+1))
  etaB <- pmvnorm(mean=c(0,0), corr=R, lower=c(l1,l2), upper=c(Inf, Inf))
  return(etaB=etaB[1])
}


