\name{poco}
\alias{poco}
\alias{poco}
\docType{package}
\title{
Power and smaple size considerations for consistency assessment
}
\description{
Companion of functions for the design of multi-regional trials considering consitency assessment:

- calculate proportion of study patients for one region

- calculate proportion of Japanese patients

- calculate power for concistency assessment

- calculate smaple size for consistency assessment

}
\details{
\tabular{ll}{
Package: \tab poco\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-03-04\cr
License: \tab GPL\cr
}
}
\author{
Andreas Kitsche <kitsche@biostat.uni-hannover.de>
}
\references{
R Functions for Sample Size and Probability Calculations for Assessing Consistency of Treatment Effects in Multi-Regional Clinical Trials, Li, M., Quan, H. Chen, J., Tanaka, Y., Ouyang, P., Luo, X., Li, G. Journal of Statistical Software, 2012 (47)

Ikeda, K. ,Bretz, F. Sample size and proportion of Japanese patients in multi-regional trials
Pharmaceutical Statistics 9: 207-216 (2010)
}

\keyword{ package }
\seealso{
}
\examples{
####Specify the parameters. Default: 
####alpha=0.025, beta=0.01, delta=0.005, 
####sigma=0.013, u=c(1,1,1,1), pi=1/4, set.pow=0.8

alpha <- 0.025
beta <- 0.01
delta <- 0.005
sigma <- 0.013
pi <- 1/4
u <- rep(1, 4)
set.pow <- 0.8  

#For some higher powers, there may not be a solution, 
#may need to increase total sample size

f11.def1(alpha, beta, delta, sigma, u, pi, set.pow)
}
