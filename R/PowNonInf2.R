#Function to calculate the any-pair power at the
#least favourable configuration (LFC)
#associated with simultaneous tests for non-inferiority
PowNonInf2 <- function(r, h, psi, k0, theta.star, alpha, n){
  eps <- 0.0001
  rho <- (psi^2)/sqrt((psi^2+1)*(psi^2+1))
  RHO <- matrix(rep(rho,r*r), nr=r)
  diag(RHO) <- rep(1,r) # correlation matrix (balanced design)
  nu <- (r+1)*(n-1)#degree of freedom
  probq <- function(q) {
    pmvt(rep(-Inf,r),rep(q, r),nu,corr=RHO,delta=rep(0,r),abseps=eps)-(1-alpha)}
  C0 <- uniroot(probq, lower=0, upper=4) $root #computes the critical point c
  theta.vec <- rep(theta.star,h)
  deltaR <- (theta.vec - psi)/(k0*sqrt(1/n + (psi^2)/n)) #non-centrality para.
  RHO.LFC <- matrix(rep(rho,h*h), nr =h)
  diag(RHO.LFC) <- rep(1,h)
  power <- pmvt(rep(C0,h),rep(Inf,h),nu,corr=RHO.LFC,delta=deltaR,abseps=eps)
  return(list=c(power=power))
}

PowNonInf2(r=3, h=1, psi=0.7, k0=0.2, theta.star=0.8, alpha=0.05, n=100)
