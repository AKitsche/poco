thetaPowCon <- function(power, n, mu, sd, n.sub=2, TreatMat = "Tukey", SubMat = "GrandMean", thetas = 1, alternative = c("two.sided", "less", "greater"),alpha=0.05){
  
  Theta <- function(theta){
    theta <- as.numeric(theta)
    PowCon(mu=mu, 
           sd=sd,
           n = n,
           n.sub=n.sub,
           TreatMat= TreatMat, 
           SubMat = SubMat,
           thetas=theta, 
           alpha=alpha, 
           alternative=alternative)[[1]]-power
  }
  Thetafinal <- as.numeric(uniroot(Theta, lower=-2, upper=2)$root)
  Power <-    PowCon(mu=mu, 
                     sd=sd,
                     n = n,
                     n.sub=n.sub,
                     TreatMat= TreatMat, 
                     SubMat = SubMat,
                     thetas=Thetafinal, 
                     alpha=alpha, 
                     alternative=alternative) 
  out <- list(power = power,
              n=n,
              NonCentrPar=Power[[3]], 
              crit = Power[[4]], 
              alternative = Power[[5]], 
              CorrMat = Power[[6]], 
              CMat=Power[[7]], 
              DMat=Power[[8]],
              thetas=Thetafinal,
              alpha = Power[[10]],
              n.sub=Power[[11]],
              TreatMat=Power[[12]],
              SubMat=Power[[13]])
  class(out) <- "Powerpoco"
  out
}

