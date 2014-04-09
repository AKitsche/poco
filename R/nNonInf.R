#function to calculate the sample size assiciated with simultaneous tests for non-inferirority
#according to 
#Dilba, G., Bretz, F., Hothorn, L.A., Guiard, V. 
#Power and sample size computations in simultaneous tests for non-inferiority based on relative margins
#Statistics in Medicine 2006 (25), 1131-1147
library (mvtnorm)
nNonInf <- function(r, g, psi, min.power, k0, theta.star, alpha, n.start){
  rho <- (psi^2)/sqrt((psi^2+1)*(psi^2+1))
  RHO <- matrix(rep(rho,r*r), nr=r)
  diag(RHO) <- rep(1,r) # correlation matrix (balanced design)
  n <- n.start
  power <- 0
  eps <- 0.00001
  while(power < min.power) {
    nu <- (r+1)*(n-1)
    probq <- function(q) {
      pmvt(rep(-Inf,r),rep(q, r),nu,corr=RHO,delta=rep(0,r),abseps=eps)-(1-alpha)}
    C0 <- uniroot(probq, lower=0, upper=4) $root #computes the critical point c
    theta.vec <- rep(theta.star,g )
    deltaR <- (theta.vec - psi)/(k0*sqrt(1/n + (psi^2)/n)) #non-centrality para.
    RHO.LFC <- matrix(rep(rho,g*g), nr =g)
    diag(RHO.LFC) <- rep(1,g)
    power <- 1-pmvt(rep(-Inf,g),rep(C0,g),nu,corr=RHO.LFC,delta=deltaR,abseps=eps)
    n <- n + 1
  }
  cbind(c(sample.size=round(n-1,0),power=round(power,4)))
}
#nNonInf(r=3, g=1, psi=0.7, min.power=0.8, k0=0.5, theta.star=0.95,
#        alpha=0.05, n.start=2)