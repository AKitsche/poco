#Function to calculate the power for the ratios of treatment effects
#the R Code is based on the Function presented in:

#Dilba, G., Bretz, F., Hothorn, L.A., Guiard, V. 
#Power and sample size computations in simultaneous tests for non-inferiority based on relative margins 
#Statistics in Medicine 25: 1131-1147 (2006)

#and the methodology presented in:
#Kitsche, A., Hothron, L.A.
#Testing for qualitative interaction using ratios of treatment differences
#Statistics in Medicine 33 (9): 1477-1489 (2014)

#Function requires:
#mu - vector of cell means from the cell means model (length should be twice the levels of the subgrouping factor)
#n - per group sample size
#sd - pooled standard deviation
#n.treat - number of treatment groups
#TreatMat - contrast matrix used for the treatment factor
#SubMat  - contrast matrix used for the subgrouping factor
#thetas - inconsistency margin
#alpha - familywise type I error rate
#alternative - test direction

#Function requires mvtnorm and mratios package
#library(multcomp)
#library(mratios)
#library(mvtnorm)

PowCon <- function(mu, n, sd, n.treat=2, TreatMat = "Tukey", SubMat = "GrandMean", thetas = 1, alternative = c("two.sided", "less", "greater"), alpha = 0.05){
  n.subgroup <- length(mu)/n.treat
  #Definition of numerator and denominator product type interaction contrast matrices for the M ratios 
  C_Treat <- contrMat(n=rep(n, n.treat), type=TreatMat)#definition of the treatment effect as the user defined difference of treatment groups
  C_Subgroup_Numerator   <- contrMatRatio(n=rep(n, n.subgroup), type = SubMat)$numC
  C_Subgroup_Denominator <- contrMatRatio(n=rep(n, n.subgroup), type = SubMat)$denC
  CMat <- kronecker(C_Subgroup_Numerator, C_Treat)#numerator interaction contrast matrix
  DMat <- kronecker(C_Subgroup_Denominator, C_Treat)#denominator interaction contrast matrix
  
  Margin.vec <- rep(thetas, nrow(CMat))
  MM <- diag(1/rep(n, ncol(CMat)))# Diagonal matrix containing reciprocals of the ni's
  ncomp <- nrow(CMat)#Number of comparisons 
  nu <- (length(mu))*(n-1)#degree of freedom
  #  Correlation matrix under H0
  CorrMat.H0 <- matrix(rep(NA,ncomp*ncomp),nrow=ncomp)
  for(i in 1:ncomp) {
    for(j in 1:ncomp) {
      CorrMat.H0[i,j] <- (Margin.vec[i]*DMat[i,] - CMat[i,])%*%MM%*%(Margin.vec[j]*DMat[j,] - CMat[j,])/
        (sqrt((Margin.vec[i]*DMat[i,] - CMat[i,])%*%MM%*%(Margin.vec[i]*DMat[i,] - CMat[i,]))*
           sqrt((Margin.vec[j]*DMat[j,] - CMat[j,])%*%MM%*%(Margin.vec[j]*DMat[j,] - CMat[j,])))
    }
  }
  #Calculating the critical value under H0
  switch(alternative, 
         two.sided = {
    crit <- qmvt(p = 1 - alpha, tail = "both.tails", df = nu, corr = CorrMat.H0)[["quantile"]]},
         less = {
    crit <- qmvt(p= 1-alpha, tail="upper", df=nu, corr=CorrMat.H0)[["quantile"]]}, 
          greater = {
    crit <- qmvt(p = 1 - alpha, tail = "lower", df = nu, corr = CorrMat.H0)[["quantile"]]}
    )
  #Calculating the non-centrality parameter vector tau
  tau <- vector(length=ncomp)
  for(i in 1:ncomp){
    tau[i] <- ((t(mu)%*%CMat[i,]) - (Margin.vec[i]*t(mu)%*%DMat[i,]))/(sd*sqrt((1+Margin.vec[i]^2)/n))}
  #calculating beta
  switch(EXPR = alternative, 
          two.sided = {
    beta <- pmvt(lower = rep(-crit, ncomp), upper = rep(crit, ncomp), delta = tau, df = nu, corr = CorrMat.H0)}, 
          less = {
    beta <- pmvt(lower = rep(crit, ncomp),  upper = rep(Inf,  ncomp), delta = tau, df = nu, corr = CorrMat.H0)}, 
          greater = {
    beta <- pmvt(lower = rep(-Inf, ncomp),  upper = rep(crit, ncomp), delta = tau, df = nu, corr = CorrMat.H0)})

  out <- list(power = 1-beta, 
              NonCentrPar=tau, 
              crit = crit, 
              alternative = alternative, 
              CorrMat = CorrMat.H0, 
              CMat=CMat, 
              DMat=DMat,
              thetas=thetas,
              alpha = alpha)
  return(out)
}
