####################################################################################################################
##																		##
## Program name: rand1.def1, using \hat{delta}^{&}											##
## 																			##
## 		 rand2.def1, using \hat{delta}												##
## 																			##
## Purpose: Calculate the unconditional and conditional probabilities to claim consistency by Definition 1 for	## 
##	    the random effect model.														##
##																			##
## Definition 1: Achieving in each region a specified proportion of the observed overall effect.			##
##																			##
## Requirement: Need to install package "mvtnorm" in R.										##
##																			##
## Date: March 16, 2010																##
##																			##
## Author: Mingyu Li																##
##																			##
## Reference: Assessment of Consistency of Treatment Effects in Multi-Regional Clinical Trials, to appear in    	##
##            DIJ, 2010																##
##																			##
####################################################################################################################


## Load "mvtnorm"

library(mvtnorm)

rand1.def1 = function(alpha=0.025, beta=0.2, delta=0.25, s=3, f=c(1,1,1), pi=1/3, tao=0)      ##using \hat{delta}^{&}

{
## Check if length(f)=s

if (length(f)!=s) stop("length(f) should be equal to s.") 

## Check if sum(f)=1

if (sum(f)!=1) stop("sum(f) should be 1.") 

z.alpha = qnorm(1-alpha)

z.beta = qnorm(1-beta)

fun1 = function(x, a)		#a=c(z.alpha, z.beta, delta, tao, f), x=N
	{
	z.alpha = a[1]
	z.beta = a[2]
	delta = a[3]
	tao = a[4]
	f = a[-(1:4)]

	sum(1/(tao^2+2/x/f))-(z.alpha+z.beta)^2/delta^2

	}

## N is the root of eq. (15) of Quan et al. 2010

N = as.numeric(uniroot(fun1, c(1,10000), tol=0.0001, a=c(z.alpha, z.beta, delta, tao, f))[1])

N = ceiling(N)

## vector of (w_1, ..., w_s)

w = 1/(tao^2+2/N/f)


mean = rep((1-pi)*delta,s)

sigma = matrix(rep(pi*(pi-2)/sum(w), s*s), s, s)

for (i in 1:s)	sigma[i,i]=sigma[i,i]+1/w[i]

lower=rep(0, s)

upper=rep(Inf,s)

## Unconditional probability to claim consistency

un = as.numeric(pmvnorm(lower=lower, upper=upper, mean=mean, sigma=sigma)[1])		


mean.con = c(mean, delta)

sigma.con = matrix(rep(0, (s+1)*(s+1)), (s+1), (s+1))

sigma.con[1:s,1:s] = sigma

sigma.con[(s+1),] = (1-pi)/sum(w)
sigma.con[,(s+1)] = (1-pi)/sum(w)

sigma.con[(s+1),(s+1)] = 1/sum(w)

lower = c(rep(0, s), z.alpha*sqrt(1/sum(w)))

upper = rep(Inf, (s+1))

## Conditional probability to claim consistency 

con = as.numeric(pmvnorm(lower=lower, upper=upper, mean=mean.con, sigma=sigma.con)[1])/(1-pnorm(z.alpha-delta/sqrt(1/sum(w))))

return(list(N=N, uncond.prob=un, cond.prob=con))

}

#########################################################

rand2.def1 = function(alpha=0.025, beta=0.2, delta=0.25, s=3, f=c(1,1,1), pi=1/3, tao=0) ## using \hat{delta}

{
## Check if length(f)=s

if (length(f)!=s) stop("length(f) should be equal to s.") 

## Check if sum(f)=1

if (sum(f)!=1) stop("sum(f) should be 1.") 

z.alpha = qnorm(1-alpha)

z.beta = qnorm(1-beta)

## Total number of patients per arm to detect overall treatment effect delta 
## with type I error=alpha, power=1-beta in random effect model

N = 1/(delta^2/2/(z.alpha+z.beta)^2-tao^2*sum(f*f)/2)

N = ceiling(N)

## vector of (w_1, ..., w_s)

w = 1/(tao^2+2/N/f)


mean = rep((1-pi)*delta, s)

c = pi*pi*sum(f*f/w)

sigma = matrix(rep(0, s*s), s, s)

for (i in 1:s)	
	for (j in 1:s) sigma[i,j] = -pi*(f[i]/w[i]+f[j]/w[j])+c

for (i in 1:s) sigma[i,i] = sigma[i,i]+1/w[i]

lower=rep(0, s)

upper=rep(Inf,s)

## Unconditional probability to claim consistency 

un = as.numeric(pmvnorm(lower=lower, upper=upper, mean=mean, sigma=sigma)[1])		


mean.con = c(mean, delta)

sigma.con = matrix(rep(0, (s+1)*(s+1)), (s+1), (s+1))

sigma.con[1:s,1:s] = sigma

sigma.con[(s+1),(1:s)] = f/w-c/pi
sigma.con[(1:s),(s+1)] = f/w-c/pi

sigma.con[(s+1),(s+1)] = c/pi/pi

lower = c(rep(0, s), z.alpha*sqrt(c/pi/pi))

upper = rep(Inf, (s+1))

## Conditional probability to claim consistency 

con = as.numeric(pmvnorm(lower=lower, upper=upper, mean=mean.con, sigma=sigma.con)[1])/(1-pnorm(z.alpha-delta/sqrt(c/pi/pi)))


return(list(N=N, uncond.prob=un, cond.prob=con))

}


#### Specify the parameters. Default: alpha=0.025 beta=0.2 delta=0.25 s=3 pi=1/3 tao=0 f=c(1/3,1/3,1/3)

alpha = 0.025

beta = 0.2

delta = 0.25		#standardized overall treatment effect

s = 3

pi = 1/s

tao = 0
# tao = 0.05
# tao = 0.1
# tao = 0.15

## Note: Length of f should be the same as s; sum(f)=1

f = c(0.1, 0.45, 0.45)

## Can set the random seed to get identical result every time

set.seed(5000)

rand1.def1(alpha, beta, delta, s, f, pi, tao)

rand2.def1(alpha, beta, delta, s, f, pi, tao)

