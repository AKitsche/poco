####################################################################################################################
##																			##
## Program name: prob.def4															##
## 																			##
## Purpose: Calculate the unconditional and conditional probabilities to claim consistency by Definition 4. 	##
##          (By definition 4, unconditional and conditional probabilities are theoretically the same)			##
##																			##
## Definition 4: Absence of significant treatment-by-region interaction								##
##																			##
## Date: March 10, 2010																##
##																			##
## Author: Mingyu Li																##
##																			##												
####################################################################################################################


prob.def4 = function(alpha=0.025, beta=0.2, delta0=0.25, s=3, f=NULL, u=NULL, eps=0.1)
{
## Check if length(f)=s and length(u)=s

if (length(f)!=s || length(u)!=s) stop("length(f) and length(u) should be equal to s.") 

## Check if sum(f)=1 and sum(f*u)=1 

if (sum(f)!=1 || sum(f*u)!=1) stop("sum(f) and sum(f*u) should be 1.") 

z.alpha = qnorm(1 - alpha)

z.beta = qnorm(1 - beta)

## Total number of patients per arm to detect overall treatment effect delta0 
## with type I error=alpha, power=1-beta in fixed design

N = 2*(z.alpha+z.beta)^2/delta0^2	

N = ceiling(N)

## treatment effect vector for s regions

delta = u*delta0

## counter for conditional probability to claim inconsistency

pow.con = 0	

## counter for unconditional probability to claim inconsistency
					
pow.un = 0	

## counter for overall significant effect probability 					

con = 0							

z.chisq = qchisq(1-eps, df=s-1)

## simulated vector of (\hat{\delta}_1, ..., \hat{\delta}_s)

x = rep(0, s)

for (i in 1:100000)

{

## \hat{\delta}_j for jth region

for (j in 1:s)	x[j] = rnorm(1, mean=delta[j], sd=sqrt(2/N/f[j]))

## \hat{\delta}

y = sum(f*x)

## Test statistics Q

Q = sum((x-y)*(x-y)*N*f)/2

## If Q>z.chisq, unconditional probability counter+1

if (Q>z.chisq) pow.un = pow.un+1

## If overall treatment effect is significant, overall significant effect probability counter+1; if Q>z.chisq too, 
## conditional probability counter+1
 
if (y>z.alpha*sqrt(2/N))
	{
	con = con+1
	if (Q>z.chisq) pow.con = pow.con+1
	}
}

## return unconditional and conditional probability to claim inconsistency 

return(list(uncond.prob = 1-pow.un/100000, cond.prob = 1-pow.con/con))

}


#### Specify the parameters. Default: alpha=0.025, beta=0.2, delta0=0.25, s=3, f=NULL, u=NULL, eps=0.1

alpha = 0.025

beta = 0.2

delta0 = 0.25		#standardized overall treatment effect

s = 3

eps = 0.1

## Note: Length of f and u should be the same as s; sum(f)=1; sum(f*u)=1. 

f = c(1/3, 1/3, 1/3)

u = c(0.25,0.55,2.2)

## Can set the random seed to get identical result every time

set.seed(5000)

## unconditional and conditional probabilities are theoretically the same

prob.def4(alpha, beta, delta0, s, f, u, eps)


