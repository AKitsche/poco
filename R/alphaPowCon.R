alphaPowCon <- function(power, n, mu, sd, n.sub=2, TreatMat = "Tukey", SubMat = "GrandMean", thetas = 1, alternative = c("two.sided", "less", "greater")){
  
  Alpha <- function(alpha){
    alpha <- as.numeric(alpha)
    PowCon(mu=mu, 
           sd=sd,
           n = n,
           n.sub=n.sub,
           TreatMat= TreatMat, 
           SubMat = SubMat,
           thetas=thetas, 
           alpha=alpha, 
           alternative=alternative)[[1]]-power
  }
  Alphafinal <- as.numeric(uniroot(Alpha, lower=0.0001, upper=0.9999)$root)
  Power <-    PowCon(mu=mu, 
                     sd=sd,
                     n = n,
                     n.sub=n.sub,
                     TreatMat= TreatMat, 
                     SubMat = SubMat,
                     thetas=thetas, 
                     alpha=Alphafinal, 
                     alternative=alternative) 
  out <- list(power = power,
              n=n,
              NonCentrPar=Power[[3]], 
              crit = Power[[4]], 
              alternative = Power[[5]], 
              CorrMat = Power[[6]], 
              CMat=Power[[7]], 
              DMat=Power[[8]],
              thetas=Power[[9]],
              alpha = Alphafinal,
              n.sub=Power[[11]],
              TreatMat=Power[[12]],
              SubMat=Power[[13]])
  class(out) <- "Powerpoco"
  out
}


