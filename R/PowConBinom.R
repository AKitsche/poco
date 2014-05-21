PowConBinom <- function(p, n, n.sub=2, TreatMat = "Tukey", SubMat = "GrandMean", rhs = 1, alternative = c("two.sided", "less", "greater"), alpha = 0.05){
  #checks
  if(length(n.sub) != 1 || !is.numeric(n.sub) | is.integer(n.sub)) {
    stop("n.sub must be a single integer value specifying the number of subgroups")
  }
  if(length(n) != 1 || !(is.numeric(n) | is.integer(n))) {
    stop("n must be a single integer value of the sample sizes per treatment-by-subgroup combination")
  }
  if(length(p) < 2 || !(is.numeric(p) | is.integer(p))) {
    stop("p must be a vector of expected proportions (at least length 2, containing integer values)")
  }
  n.subgroup <- n.sub
  n.treat    <- length(p)/n.sub
  #Definition of numerator and denominator product type interaction contrast matrices for the M ratios 
  if(n.sub==1){
    CMat <- contrMatRatio(n=rep(n, n.treat), type = TreatMat)$numC
    DMat <- contrMatRatio(n=rep(n, n.treat), type = TreatMat)$denC
  }else{
    C_Treat <- contrMat(n=rep(n, n.treat), type=TreatMat)#definition of the treatment effect as the user defined difference of treatment groups
    C_Subgroup_Numerator   <- contrMatRatio(n=rep(n, n.subgroup), type = SubMat)$numC
    C_Subgroup_Denominator <- contrMatRatio(n=rep(n, n.subgroup), type = SubMat)$denC
    CMat <- kronecker(C_Subgroup_Numerator, C_Treat)#numerator interaction contrast matrix
    DMat <- kronecker(C_Subgroup_Denominator, C_Treat)#denominator interaction contrast matrix
  }
  if(length(rhs) != 1 & length(rhs) != nrow(CMat)) {
    stop("rhs must either be a single numeric value or a numeric vector which length is equal to the number of comparions")
  }
  #definition of the vector of inconsistency margins
  if(length(rhs) == 1){
    Margin.vec <- rep(rhs, nrow(CMat))
  }else{
    Margin.vec <- rhs
  }
  

  
  #MM <- diag(1/rep(n, ncol(CMat)))# Diagonal matrix containing reciprocals of the ni's
  variances <- (p*(1-p))/n
  MM <- diag(variances)
  nu <- (length(p))*(n-1)#degree of freedom
  ncomp <- nrow(CMat)#Number of comparisons 
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
  alternative <- match.arg(alternative)
#   switch(alternative, 
#          two.sided = {
#            crit <- qmvnorm(p = 1 - alpha, tail = "both.tails", sigma = CorrMat.H0)[["quantile"]]},         less = {
#            crit <- qmvnorm(p= 1-alpha, tail="upper", sigma = CorrMat.H0)[["quantile"]]}, 
#          greater = {
#            crit <- qmvnorm(p = 1 - alpha, tail = "lower", sigma = CorrMat.H0)[["quantile"]]}
#   )
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
  for (i in 1:ncomp){
    #Ratio.Estimate[i] <- (CMat[i,]%*%t(p))/(DMat[i,]%*%t(p))
    tau[i] <- ((CMat[i,] - Margin.vec[i]*DMat[i,])%*%p)/
      sqrt((CMat[i,] - Margin.vec[i]*DMat[i,])%*%MM%*%(CMat[i,] - Margin.vec[i]*DMat[i,]))
  }
  #calculating beta
#   switch(EXPR = alternative, 
#          two.sided = {
#            beta <- pmvnorm(lower = rep(-crit, ncomp), upper = rep(crit, ncomp), mean = tau, sigma = CorrMat.H0)}, 
#          less = {
#            beta <- pmvnorm(lower = rep(crit, ncomp),  upper = rep(Inf,  ncomp), mean = tau, sigma = CorrMat.H0)}, 
#          greater = {
#            beta <- pmvnorm(lower = rep(-Inf, ncomp),  upper = rep(crit, ncomp), mean = tau, sigma = CorrMat.H0)})
#   
switch(EXPR = alternative, 
       two.sided = {
         beta <- pmvt(lower = rep(-crit, ncomp), upper = rep(crit, ncomp), delta = tau, df = nu, corr = CorrMat.H0)}, 
       less = {
         beta <- pmvt(lower = rep(crit, ncomp),  upper = rep(Inf,  ncomp), delta = tau, df = nu, corr = CorrMat.H0)}, 
       greater = {
         beta <- pmvt(lower = rep(-Inf, ncomp),  upper = rep(crit, ncomp), delta = tau, df = nu, corr = CorrMat.H0)})

  out <- list(power = 1-beta,
              n=n,
              NonCentrPar=tau, 
              crit = crit, 
              alternative = alternative, 
              CorrMat = CorrMat.H0, 
              CMat=CMat, 
              DMat=DMat,
              rhs=rhs,
              alpha = alpha,
              n.sub=n.sub,
              TreatMat=TreatMat,
              SubMat=SubMat)
  class(out) <- "Powerpoco"
  out
}
