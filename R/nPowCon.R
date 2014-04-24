nPowCon <- function(min.power, n.start, mu, sd, n.sub=2, TreatMat = "Tukey", SubMat = "GrandMean", thetas = 1, alternative = c("two.sided", "less", "greater"), alpha = 0.05){
  power <- 0
  n <- n.start

  #checks
  if(length(sd) != 1 || !is.numeric(sd)) {
    stop("sd must be a single numeric value")
  }
  if(length(alpha) != 1 || !is.numeric(alpha) | alpha <= 0 | alpha >= 1) {
    stop("alpha must be a single numeric value between 0 and 1")
  }
  if(length(n.sub) != 1 || !is.numeric(n.sub) | is.integer(n.sub)) {
    stop("n.sub must be a single integer value specifying the number of subgroups")
  }
  if(length(n) != 1 || !(is.numeric(n) | is.integer(n))) {
    stop("n must be a single integer value of the sample sizes per treatment-by-subgroup combination")
  }
  if(length(mu) < 2 || !(is.numeric(mu) | is.integer(mu))) {
    stop("mu must be a vector of expected means (at least length 2, containing numeric or integer values)")
  }
  n.subgroup <- n.sub
  n.treat    <- length(mu)/n.sub
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
  if(length(thetas) != 1 & length(thetas) != nrow(CMat)) {
    stop("thetas must either be a single numeric value or a numeric vector which length is equal to the number of comparions")
  }
  #definition of the vector of inconsistency margins
  if(length(thetas) == 1){
    Margin.vec <- rep(thetas, nrow(CMat))
  }else{
    Margin.vec <- thetas
  }
  #loop for iterative search of the sample size
  while(power < min.power){
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
    alternative <- match.arg(alternative)
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
    power <- 1-beta
    n <- n+1
  }
  out <- list(power = 1-beta,
              n=n-1,
              NonCentrPar=tau, 
              crit = crit, 
              alternative = alternative, 
              CorrMat = CorrMat.H0, 
              CMat=CMat, 
              DMat=DMat,
              thetas=thetas,
              alpha = alpha,
              n.sub=n.sub,
              TreatMat=TreatMat,
              SubMat=SubMat)
  class(out) <- "Powerpoco"
  out
}