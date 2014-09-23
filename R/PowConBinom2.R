PowConBinom <- function(p, n, n.sub=2, TreatMat = "Tukey", SubMat = "GrandMean", rhs = 1, 
                        alternative = c("two.sided", "less", "greater"), alpha = 0.05,
                        type=c("anypair","allpair","global")){
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
  if(any(p > 1)  | any(p < 0)){
    stop("all elements of p must be numeric values between 0 and 1")
  }
  
  if(length(n.sub) != 1 || !is.numeric(n.sub) | is.integer(n.sub)) {
    stop("n.sub must be a single integer value specifying the number of subgroups")
  }
  if(length(n) != length(p) & length(n) != 1) {
    stop("n must be of the same length as p, correponding to the sample size for each treatment-by-subgroup combination")
  }
  if(length(p) < 2 || !(is.numeric(p) | is.integer(p))) {
    stop("p must be a vector of expected proportions (at least length 2, containing integer values)")
  }
  n.subgroup <- n.sub
  n.treat    <- length(p)/n.sub
  
  #determining the sample sizes 
  if(length(n)==1){
    nTreat <- rep(n, n.treat)
    nSub   <- rep(n, n.subgroup)
  }else{
    nTreat <- vector(length=n.treat)#sample size for each treatment group
    IteratorTreat <- seq(1,length(p), by=n.treat)
    #IteratorTreat <- c(IteratorTreat,length(p)+1)
    for(i in 1:n.treat){
      #k <- IteratorTreat[i]
      nTreat[i] <- sum(n[IteratorTreat+i-1])
    }
    
    nSub   <- vector(length=n.sub)#sample size for each subgroup
    IteratorSub <- seq(1,length(p), by=n.treat)
    IteratorSub <- c(IteratorSub, length(p)+1)
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
    C_Treat <- contrMat(n=nTreat, type=TreatMat)#definition of the treatment effect as the user defined difference of treatment groups
    C_Subgroup_Numerator1   <- contrMatRatio(n=nSub, type = SubMat)$numC
    C_Subgroup_Denominator1   <- contrMatRatio(n=nSub, type = SubMat)$denC
    C_Subgroup_Denominator2 <- contrMatRatio(n=n, type = SubMat)$denC#contrMatRatio(n=n, type = SubMat)$denC
    CMat <- kronecker(C_Subgroup_Numerator1, C_Treat)#numerator interaction contrast matrix
    DMat <- t(C_Subgroup_Denominator2[1,]*t(kronecker(matrix(rep(1,length(nSub)*nrow(C_Subgroup_Numerator1)),nrow=nrow(C_Subgroup_Numerator1)),C_Treat)))#denominator interaction contrast matrix
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
  nu <- sum(n-1)#degree of freedom
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
  #   switch(EXPR = alternative, 
  #          two.sided = {
  #            beta <- pmvnorm(lower = rep(-crit, ncomp), upper = rep(crit, ncomp), mean = tau, sigma = CorrMat.H0)}, 
  #          less = {
  #            beta <- pmvnorm(lower = rep(crit, ncomp),  upper = rep(Inf,  ncomp), mean = tau, sigma = CorrMat.H0)}, 
  #          greater = {
  #            beta <- pmvnorm(lower = rep(-Inf, ncomp),  upper = rep(crit, ncomp), mean = tau, sigma = CorrMat.H0)})
  #   
  #calculating beta for different power definitions
  ptype <- match.arg(type)
  switch(EXPR = ptype,
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
                      RT <- rmvtFS(n = nsim, delta = tau[whichHA], df = nu, sigma = matrix(CorrMat.H0[whichHA, whichHA]), method = "svd")
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
              rhs=rhs,
              alpha = alpha,
              n.sub=n.sub,
              TreatMat=TreatMat,
              SubMat=SubMat,
              type=ptype)
  class(out) <- "Powerpoco"
  out
}
