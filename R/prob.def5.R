####################################################################################################################
##																			##
## Program name: prob.def5															##
## 																			##
## Purpose: Calculate the unconditional and conditional probabilities to claim consistency by Definition 5.		##
##          (By definition 5, unconditional and conditional probabilities are theoretically the same)			##
##																			##
## Definition 5: Lack of significnat difference for any regions from the overall						##
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

prob.def5 = function(alpha=0.025, beta=0.2, eta=0.1/3, delta=0.25, s=3, f=NULL, u=NULL)

{
## Check if length(f)=s and length(u)=s

if (length(f)!=s || length(u)!=s) stop("length(f) and length(u) should be equal to s.") 

## Check if sum(f)=1 and sum(f*u)=1 

if (sum(f)!=1 || sum(f*u)!=1) stop("sum(f) and sum(f*u) should be 1.") 

z.alpha = qnorm(1-alpha)

z.beta = qnorm(1-beta)

## Total number of patients per arm to detect overall treatment effect delta 
## with type I error=alpha, power=1-beta in fixed design

N=2*(z.alpha+z.beta)^2/delta^2	

N=ceiling(N)

mean = (u-1)*delta

sigma = matrix(rep(-2/N, s*s), s, s)

for (i in 1:s)	sigma[i,i]=sigma[i,i]+2/N/f[i]

## eta corresponds to alpha' in definition 5

z.eta = qnorm(1-eta)

lower = -z.eta*sqrt(2*(1-f)/N/f)

upper = rep(Inf,s)

## Unconditional probability to claim consistency

un = as.numeric(pmvnorm(lower=lower, upper=upper, mean=mean, sigma=sigma)[1])	


mean.con = c(mean, delta)

sigma.con = matrix(rep(0, (s+1)*(s+1)), (s+1), (s+1))

sigma.con[1:s,1:s] = sigma

sigma.con[(s+1),] = 0

sigma.con[,(s+1)] = 0

sigma.con[(s+1),(s+1)] = 2/N

lower = c(-z.eta*sqrt(2*(1-f)/N/f), z.alpha*sqrt(2/N))

upper = rep(Inf, (s+1))

## Conditional probability to claim consistency

con = as.numeric(pmvnorm(lower=lower, upper=upper, mean=mean.con, sigma=sigma.con)[1])/(1-pnorm(z.alpha-delta/sqrt(2/N)))


return(list(uncond.prob = un, cond.prob = con))

}


#### Specify the parameters. 
###  Default: alpha=0.025, beta=0.2, eta=0.1/3, delta=0.25, s=3, f=NULL, u=NULL

alpha = 0.025

beta = 0.2

delta = 0.25		#standardized overall treatment effect

s = 3

eta = 0.1/3

## Note: Length of f and u should be the same as s; sum(f)=1; sum(f*u)=1.

f = c(1/3, 1/3, 1/3)

u = c(0.25, 0.55, 2.2)

## Can set the random seed to get identical result every time

set.seed(5000)

## unconditional and conditional probabilities are theoretically the same

prob.def5(alpha, beta, eta, delta, s, f, u)


