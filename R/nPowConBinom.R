nPowConBinom <- function(min.power, p, n.sub=2, TreatMat = "Tukey", SubMat = "GrandMean", rhs = 1, alternative = c("two.sided", "less", "greater"), alpha = 0.05){
  
  if(length(min.power) != 1 || !is.numeric(min.power) | min.power <= 0 | min.power >= 1) {
    stop("min.power must be a single numeric value between 0 and 1")
  }
  
  samplesize <- function(n){
    n <- as.integer(n)
    PowConBinom(p=p, 
                n=n, 
                n.sub=n.sub, 
                TreatMat = TreatMat, 
                SubMat = SubMat, 
                rhs = rhs, 
                alternative = alternative, 
                alpha = alpha)[[1]]-min.power
  }
  nfinal <- as.integer(uniroot(samplesize, lower=10, upper=1000)$root)
  Power <-        PowConBinom(p=p, 
                              n=nfinal, 
                              n.sub=n.sub, 
                              TreatMat = TreatMat, 
                              SubMat = SubMat, 
                              rhs = rhs, 
                              alternative = alternative, 
                              alpha = alpha) 
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