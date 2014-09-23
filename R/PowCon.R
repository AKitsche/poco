#Function to calculate the power for the ratios of treatment effects
#the R Code is based on the Function presented in:

#Dilba, G., Bretz, F., Hothorn, L.A., Guiard, V. 
#Power and sample size computations in simultaneous tests for non-inferiority based on relative margins 
#Statistics in Medicine 25: 1131-1147 (2006)

#and the methodology presented in:
#Kitsche, A., Hothron, L.A.
#Testing for qualitative interaction using ratios of treatment differences
#Statistics in Medicine 33 (9): 1477-1489 (2014)

#Several power definitions are considered: global power, complete (all-pairs) power and minimal (any-pairs) power
#based on the code from the function powermcpt in the package MCPAN

#Function requires:
#mu - vector of cell means from the cell means model (length should be twice the levels of the subgrouping factor)
#n - per group sample size
#sd - pooled standard deviation
#n.sub - number of subgroups
#TreatMat - contrast matrix used for the treatment factor
#SubMat  - contrast matrix used for the subgrouping factor
#thetas - inconsistency margin
#alpha - familywise type I error rate
#alternative - test direction

#Function requires mvtnorm and mratios package
#library(multcomp)
#library(mratios)
#library(mvtnorm)

PowCon <- function(mu, n, sd, n.sub=2, TreatMat = "Tukey", SubMat = "GrandMean", 
                   thetas = 1, alternative = c("two.sided", "less", "greater"), 
                   alpha = 0.05, type=c("anypair","allpair","global")){
  #function to simulate random numbers from a multivariate normal distribution according to Frank Schaarschmidt in powermcpt 
  rmvtFS <- function(n, sigma = diag(2), df = 1, delta = rep(0, nrow(sigma)), type = c("shifted", "Kshirsagar"), method = c("eigen", "svd", "chol")) {
    if (length(delta) != nrow(sigma)) 
      stop("delta and sigma have non-conforming size")
    if (df == 0) 
      return(rmvnorm(n, mean = delta, sigma = sigma, method = method))
    type <- match.arg(type)
    if (type == "Kshirsagar") 
      return(rmvnorm(n, mean = delta, sigma = sigma, method = method)/sqrt(rchisq(n, df)/df))
    if (type == "shifted") {
      sims <- rmvnorm(n, sigma = sigma, method = method)/sqrt(rchisq(n,  df)/df)
      return(sweep(sims, 2, delta, "+"))
    }
  }
  
  #checks
  if(length(sd) != 1 || !is.numeric(sd)) {
    stop("sd must be a single numeric value")
  }
  #  if(length(alpha) != 1 || !is.numeric(alpha) | alpha <= 0 | alpha >= 1) {
  #    stop("alpha must be a single numeric value between 0 and 1")
  #  }
  if(length(n.sub) != 1 || !is.numeric(n.sub) | is.integer(n.sub)) {
    stop("n.sub must be a single integer value specifying the number of subgroups")
  }
  if(length(n) != length(mu) & length(n) != 1) {
    stop("n must be a single integer value of the sample sizes per treatment-by-subgroup combination")
  }
  if(length(mu) < 2 || !(is.numeric(mu) | is.integer(mu))) {
    stop("mu must be a vector of expected means (at least length 2, containing numeric or integer values)")
  }
  n.subgroup <- n.sub
  n.treat    <- length(mu)/n.sub
  
  #determining the sample sizes 
  if(length(n)==1){
    nTreat <- rep(n, n.treat)
    nSub   <- rep(n, n.subgroup)
  }else{
    nTreat <- vector(length=n.treat)#sample size for each treatment group
    IteratorTreat <- seq(1,length(mu), by=n.treat)
    #IteratorTreat <- c(IteratorTreat,length(p)+1)
    for(i in 1:n.treat){
      #k <- IteratorTreat[i]
      nTreat[i] <- sum(n[IteratorTreat+i-1])
    }
    
    nSub   <- vector(length=n.sub)#sample size for each subgroup
    IteratorSub <- seq(1,length(mu), by=n.treat)
    IteratorSub <- c(IteratorSub, length(mu)+1)
    for(i in 1:n.subgroup){
      k <- IteratorSub[i]
      nSub[i] <- sum(n[k:(IteratorSub[i+1]-1)]) 
    }    
  }
  
  #Definition of numerator and denominator product type interaction contrast matrices for the M ratios 
  if(n.sub==1){
    CMat <- contrMatRatio(n=nTreat, type = TreatMat)$numC
    DMat <- contrMatRatio(n=nTreat, type = TreatMat)$denC
  }else{
    #Definition of numerator and denominator product type interaction contrast matrices for the M ratios 
    C_Treat <- contrMat(n=nTreat, type=TreatMat)#definition of the treatment effect as the user defined difference of treatment groups
    C_Subgroup_Numerator   <- contrMatRatio(n=nSub, type = SubMat)$numC
    CMat <- kronecker(C_Subgroup_Numerator, C_Treat)#numerator interaction contrast matrix
    #denominator contrast matrix
    #identification of required sample sizes to defines DMat for each contrast of C_Treat
    IDn <- matrix(rep(C_Treat, times=n.subgroup), nrow=nrow(C_Treat))
    #IDn - ID matrix to define the required sample sizes
    N <- matrix(rep(0, ncol(CMat)*nrow(C_Treat)), nrow=nrow(C_Treat))
    #N - matrix of sample sizes required for the denominator matrix for each contrast in C_Treat
    for(i in 1:nrow(C_Treat)){
      N[i,which(IDn[i,]!=0)] <- n[which(IDn[i,]!=0)]
    }
    #summarize the sample sizes according to the subgrouping factor
    Nsum   <- rep(NA, times=nrow(C_Treat)*n.sub)
    ID_Sub <- seq(1,length(N)+1, by=n.treat)#identifies the starting value for each subgrouping sample size in N
    for(i in 1:length(Nsum)){
      Nsum[i] <- sum(t(N)[seq(ID_Sub[i],ID_Sub[i+1]-1)])
    }
    NSub <- matrix(Nsum, nrow=nrow(C_Treat),byrow=TRUE)
    C_Subgroup_Denominator   <- matrix(rep(0,length(IDn)),nrow=nrow(C_Treat))
    for(i in 1:nrow(C_Treat)){
      C_Subgroup_Denominator[i,] <- kronecker(contrMatRatio(n=NSub[i,], type = SubMat)$denC[1,],C_Treat[i,])
    }
    DMat <- kronecker(rep(1,nrow(C_Subgroup_Numerator)), C_Subgroup_Denominator)
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
  
  
  if(length(n)==1){
    MM <- diag(1/rep(n, ncol(CMat)))# Diagonal matrix containing reciprocals of the ni's
  }else{
    MM <- diag(1/n)
  }  
  ncomp <- nrow(CMat)#Number of comparisons 
  if(length(n)==1){
    nu <- (length(mu))*(n-1)#degree of freedom
  }else{
    nu <- sum(n-1)
  }
  
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
  if(length(n)==1){
    for(i in 1:ncomp){
      tau[i] <- ((t(mu)%*%CMat[i,]) - (Margin.vec[i]*t(mu)%*%DMat[i,]))/(sd*sqrt((1+Margin.vec[i]^2)/n))}
  }else{
    for(i in 1:ncomp){
      tau[i] <- ((t(mu)%*%CMat[i,]) - (Margin.vec[i]*t(mu)%*%DMat[i,]))/(sd*sqrt((1+Margin.vec[i]^2)/nSub[i]))}
  }
  
  #calculating beta for different power definitions
  ptype <- match.arg(type)
  switch(EXPR = ptype, 
  #calculating beta for global power
  global = {
    switch(EXPR = alternative, 
           two.sided = {
             beta <- pmvt(lower = rep(-crit, ncomp), upper = rep(crit, ncomp), delta = tau, df = nu, corr = CorrMat.H0)
             whichHA <- which(tau != 0)}, 
           less = {
             beta <- pmvt(lower = rep(crit, ncomp),  upper = rep(Inf,  ncomp), delta = tau, df = nu, corr = CorrMat.H0)
             whichHA <- which(tau < 0)}, 
           greater = {
             beta <- pmvt(lower = rep(-Inf, ncomp),  upper = rep(crit, ncomp), delta = tau, df = nu, corr = CorrMat.H0)
             whichHA <- which(tau > 0)})
  },
  #calculating beta for any-pair (minimum) power
  anypair = {
    switch(EXPR = alternative,
           two.sided = {
             whichHA <- which(tau != 0)
             NHA <- length(whichHA)
             if(NHA < 1){
               warning("All contrasts are under their corresponding null hypothesis, any-pair power can not be calculated")
               beta <- 1-alpha
             }else{
               beta <- pmvt(lower = rep(-crit, NHA), upper = rep(crit, NHA), delta = tau[whichHA], df = nu, corr = CorrMat.H0[whichHA,whichHA])
             }
           },
           less = {
             whichHA <- which(tau < 0)
             NHA <- length(whichHA)
             if(NHA < 1){
               warning("All contrasts are under their corresponding null hypothesis, any-pair power can not be calculated")
               beta <- 1-alpha
             }else{
               beta <- pmvt(lower = rep(crit, NHA), upper = rep(Inf, NHA), delta = tau[whichHA], df = nu, corr = CorrMat.H0[whichHA,whichHA])
             }
           },
           greater = {
             whichHA <- which(tau > 0)
             NHA <- length(whichHA)
             if(NHA < 1){
               warning("All contrasts are under their corresponding null hypothesis, any-pair power can not be calculated")
               beta <- 1-alpha
             }else{
               beta <- pmvt(lower = rep(-Inf, NHA), upper = rep(crit, NHA), delta = tau[whichHA], df = nu, corr = CorrMat.H0[whichHA,whichHA])
             }
           })
  },
  allpair = (
    switch(EXPR = alternative, 
           two.sided = {
             whichHA <- which(tau != 0)
             NHA <- length(whichHA)
             if (NHA < 1) {
               warning("All contrasts are under the corresponding null hypotheses, all pair power can not be calculated.")
               beta <- 1 - alpha
             } else {
               nsim <- 10000
               RT <- rmvtFS(n = nsim, delta = tau[whichHA], df = nu, sigma = CorrMat.H0[whichHA, whichHA], method = "svd")
               nreject <- sum(apply(RT, 1, function(x) {
                 min(abs(x))
               }) > abs(crit))
               beta <- 1 - (nreject/nsim)
               simerror <- sqrt(beta * (1 - beta)/nsim)
               attr(beta, which = "simerror") <- simerror
             }
           },
           less = {
             whichHA <- which(tau < 0)
             NHA <- length(whichHA)
             if(NHA < 1){
               warning("All contrasts are under their corresponding null hypothesis, any-pair power can not be calculated")
               beta <- 1-alpha
             }else{
               beta <- beta <- 1-pmvt(lower = rep(-Inf, NHA),  upper = rep(crit,  NHA), delta = tau[whichHA], df = nu, corr = CorrMat.H0[whichHA,whichHA])}
             },
           greater = {
             whichHA <- which(tau > 0)
             NHA <- length(whichHA)
             if(NHA < 1){
               warning("All contrasts are under their corresponding null hypothesis, any-pair power can not be calculated")
               beta <- 1-alpha
             }else{
               beta <- 1-pmvt(lower = rep(crit, NHA),  upper = rep(Inf, NHA), delta = tau[whichHA], df = nu, corr = CorrMat.H0[whichHA,whichHA])
             }
           })
    )
  )
  
  
  out <- list(power = 1-beta,
              n=n,
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
              SubMat=SubMat,
              type=ptype)
  class(out) <- "Powerpoco"
  out
}
