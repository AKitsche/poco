####################################################################################################################
##																			##
## Program name: prob.def2															##
## 																			##
## Purpose: Calculate the unconditional and conditional probability to claim consistency by Definition 2.		##
##																			##
## Definition 2: Observing region effects that exceed a pre-specified effect size.						##
##																			##
## Requirement: Need to install package "mvtnorm" in R.										##
##																			##
## Date: March 10, 2010																##
##																			##
## Author: Mingyu Li																##
##																			##
####################################################################################################################


## Load "mvtnorm"

library(mvtnorm)

prob.def2 = function(alpha=0.025, beta=0.2, delta=0.25, s=3, f=NULL, u=NULL, b=0.25/3)

{
## Check if length(f)=s and length(u)=s

if (length(f)!=s || length(u)!=s) stop("length(f) and length(u) should be equal to s.") 

## Check if sum(f)=1 and sum(f*u)=1 

if (sum(f)!=1 || sum(f*u)!=1) stop("sum(f) and sum(f*u) should be 1.") 

z.alpha = qnorm(1-alpha)

z.beta = qnorm(1-beta)

## Total number of patients per arm to detect overall treatment effect delta 
## with type I error=alpha, power=1-beta in fixed design.

N = 2*(z.alpha+z.beta)^2/delta^2		

N = ceiling(N)

mean = u*delta

sigma = 2/N/f

## Unconditional probability to claim consistency 

un = prod(pnorm(b, mean = mean, sd = sqrt(sigma), lower.tail = FALSE))


mean.con = c(mean, delta)

sigma.con = matrix(rep(0, (s+1)*(s+1)), (s+1), (s+1))

for (i in 1:s) sigma.con[i,i] = 2/N/f[i]

sigma.con[(s+1),] = 2/N
sigma.con[,(s+1)] = 2/N

lower = c(rep(b, s), z.alpha*sqrt(2/N))

upper=rep(Inf, (s+1))

## Conditional probability to claim consistency

con = as.numeric(pmvnorm(lower=lower, upper=upper, mean=mean.con, sigma=sigma.con)[1])/(1-pnorm(z.alpha-delta/sqrt(2/N)))


return(list(uncond.prob=un, cond.prob=con))

}

