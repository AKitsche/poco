nPowCon <- function(min.power, mu, sd, n.sub=2, TreatMat = "Tukey", SubMat = "GrandMean", thetas = 1, alternative = c("two.sided", "less", "greater"), alpha = 0.05){
  
  if(length(min.power) != 1 || !is.numeric(min.power) | min.power <= 0 | min.power >= 1) {
    stop("min.power must be a single numeric value between 0 and 1")
  }
  
  samplesize <- function(n){
    n <- as.integer(n)
    PowCon(mu=mu, 
           sd=sd,
           n = n,
           n.sub=n.sub,
           TreatMat= TreatMat, 
           SubMat = SubMat,
           thetas=thetas, 
           alpha=alpha, 
           alternative=alternative)[[1]]-min.power
  }
  nfinal <- as.integer(uniroot(samplesize, lower=2, upper=1000)$root)
  Power <-    PowCon(mu=mu, 
                     sd=sd,
                     n = nfinal,
                     n.sub=n.sub,
                     TreatMat= TreatMat, 
                     SubMat = SubMat,
                     thetas=thetas, 
                     alpha=alpha, 
                     alternative=alternative) 
  out <- list(power = Power[[1]],
              n=nfinal,
              NonCentrPar=Power[[3]], 
              crit = Power[[4]], 
              alternative = Power[[5]], 
              CorrMat = Power[[6]], 
              CMat=Power[[7]], 
              DMat=Power[[8]],
              thetas=Power[[9]],
              alpha = Power[[10]],
              n.sub=Power[[11]],
              TreatMat=Power[[12]],
              SubMat=Power[[13]])
  class(out) <- "Powerpoco"
  out
}
  
  
  

