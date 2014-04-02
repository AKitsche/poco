#Function to calculate the probability of observing a consistent result
#joint probability of detecting a significant difference between the study drug and
#placebo in all patients and showing a consistent result between all and Japanese patients
#this is denoted by eta_{C} in:

#Sample size and proportion of Japanese patients in multi-regional trials
#Kimitoshi Ikeda and Frank Bretz
#Pharmaceutical Statistics 9: 207–216 (2010)
library(mvtnorm)
powC <- function(delta=0.125, omega=0.5 , sigma=1 ,alpha=0.025, beta=0.2, p=0.1){
  Zalpha <- qnorm(p=1-alpha)
  Zbeta  <- qnorm(p=1-beta)
  na <- nb <- as.integer((2*(Zalpha+Zbeta)^2*sigma^2)/(delta^2),0)
  #x      <- (delta*(1-omega)*sqrt(na*nb*p))/(sigma*sqrt((na+nb)*(p*omega^2-2*p*omega+1)))
  #etaA   <- dnorm(x=x)
  m1      <- Zalpha-(delta/sigma)*sqrt((na*nb)/(na+nb))
  m2      <- -(delta*(1-omega)*sqrt(na*nb*p))/(sigma*sqrt(na+nb)*(p*omega^2-2*p*omega+1))
  R       <- diag(c(1,1))
  R[1,2] <-  R[2,1]  <- (sqrt(p)*(1-omega))/(sqrt(p*omega^2-2*p*omega+1))
  etaC <- pmvnorm(mean=c(0,0), corr=R, lower=c(m1,m2), upper=c(Inf, Inf))
  return(etaC=etaC[1])
}
