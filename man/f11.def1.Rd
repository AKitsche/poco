\name{f11.def1}
\alias{f11.def1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculate proportion of study patients for one region
}
\description{
Purpose: To calculate the minimum proportion of study patients for one region (f1) to have an specified chance of showing consistency by Definition 1.

Definition 1: Achieving in each region a specified proportion of the observed overall effect.

Under the restriction f1 < f2 = f3 = f4
}
\usage{
f11.def1(alpha = 0.025, beta = 0.01, delta = 0.005, 
         sigma = 0.013, u = c(1, 1, 1, 1), pi = 1/4, set.pow = 0.8)
}
\arguments{
  \item{alpha}{
One-sided signifcance level for testing the overall treatment effect.
}
  \item{beta}{
1-\eqn{\beta} power for detecting an overall treatment effect of delta.
}
  \item{delta}{
The standardized overall treatment effect for calculating the (1-\eqn{\beta}) overall power.
}
  \item{sigma}{
varianz
}
  \item{u}{
Vector of \eqn{(f_{1},f_{2},\ldots,f_{s})}, where \eqn{f_{i}} is the proportion of patients within region \eqn{i} and \eqn{\sum_{i=1}^{s}f_{i}=1}
}
  \item{pi}{
A threshold that each region should achieve at least pi of the observed overall effect
in order to claim consistency.
}
  \item{set.pow}{
probability of demonstrating consistency, e.g. 80\%
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
R Functions for Sample Size and Probability Calculations for Assessing Consistency of Treatment Effects in Multi-Regional Clinical Trials, Li, M., Quan, H. Chen, J., Tanaka, Y., Ouyang, P., Luo, X., Li, G. Journal of Statistical Software, 2012 (47)
}
\author{
Mingyu Li
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#### Specify the parameters. Default: 
####alpha=0.025, beta=0.01, delta=0.005, sigma=0.013, u=c(1,1,1,1), pi=1/4, set.pow=0.8

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

#At the design stage, we would like to determine the minimum 
#required proportion of sample size for a particular region, 
#say Region 1, so that there was an 80 percent probability of 
#demonstrating consistency. If the sample sizes of the other 
#3 regions were the same (f1 < f2 = f3 = f4) and conditional 
#probability was the concern, then using Defnition 1 with pi=1/4
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
