####################################################################################################################
##																			##
## Program name: prob.def3															##
## 																			##
## Purpose: Calculate the unconditional and conditional probability to claim consistency by Definition 3.		##
##																			##
## Definition 3: Demonstrating that region effects exceed a proportion of the overall effect using hypothesis	## 
##		testing																##
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

prob.def3 = function(alpha=0.025, beta=0.2, pi=0, delta=0.25, s=3, f=NULL, u=NULL, eta=0.3)

{
## Check if length(f)=s and length(u)=s

if (length(f)!=s || length(u)!=s) stop("length(f) and length(u) should be equal to s.") 

## Check if sum(f)=1 and sum(f*u)=1 

if (sum(f)!=1 || sum(f*u)!=1) stop("sum(f) and sum(f*u) should be 1.") 

z.alpha = qnorm(1-alpha)

z.beta = qnorm(1-beta)

## Total number of patients per arm to detect overall treatment effect delta 
## with type I error=alpha, power=1-beta in fixed design.

N=2*(z.alpha+z.beta)^2/delta^2		

N=ceiling(N)

## eta corresponds to alpha' in definition 3

z.eta = qnorm(1-eta)

mean = (u-pi)*delta

sigma = matrix(rep((-2*pi+pi^2)*2/N, s*s), s, s)

for (i in 1:s)	sigma[i,i]=sigma[i,i]+2/N/f[i]

lower = z.eta*sqrt(2*(1/f-2*pi+pi^2)/N)

upper = rep(Inf,s)

## Unconditional probability to claim consistency

un = as.numeric(pmvnorm(lower=lower, upper=upper, mean=mean, sigma=sigma)[1])		


mean.con = c(mean, delta)

sigma.con = matrix(rep(0, (s+1)*(s+1)), (s+1), (s+1))

sigma.con[1:s,1:s] = sigma

sigma.con[(s+1),] = (1-pi)*2/N
sigma.con[,(s+1)] = (1-pi)*2/N

sigma.con[(s+1),(s+1)] = 2/N

lower = c(z.eta*sqrt(2*(1/f-2*pi+pi^2)/N), z.alpha*sqrt(2/N))

upper = rep(Inf, (s+1))

## Conditional probability to claim consistency

con = as.numeric(pmvnorm(lower=lower, upper=upper, mean=mean.con, sigma=sigma.con)[1])/(1-pnorm(z.alpha-delta/sqrt(2/N)))


return(list(uncond.prob=un, cond.prob=con))

}





