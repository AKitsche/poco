\name{prob.def1}
\alias{prob.def1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the probability to claim consistency
}
\description{
Purpose: Calculate the unconditional and conditional probabilities to claim consistency by Definition 1.

Definition 1: Achieving in each region a specified proportion of the observed overall effect. 
}
\usage{
prob.def1(alpha = 0.025, beta = 0.2, delta = 0.25, 
          s = 3, f = NULL, u = NULL, pi = 1/3)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{s}{
Number of regions
}
  \item{f}{
Vector of \eqn{(f_{1},f_{2},\ldots,f_{s})}, where \eqn{f_{i}} is the proportion of patients within region \eqn{i} and \eqn{\sum_{i=1}^{s}f_{i}=1}
}
  \item{u}{
Vector of \eqn{(u_{1},u_{2},\ldots,u_{s})}, where \eqn{u_{i}} is the ratio of the treatment effect for region \eqn{i} to the overall effect, \eqn{u_{i}=\delta_{i}/\delta}, and \eqn{\sum_{i=1}^{s}f_{i}u_{i}=1}
}
  \item{pi}{
A threshold that each region should achieve at least pi of the observed overall eect
in order to claim consistency.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{uncond.prob}{unconditional probability to claim consistency}
\item{cond.prob}{conditional probability to claim consistency given a signifcant overall treat-
ment effect}

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
#alpha=0.025, beta=0.2, delta=0.25, s=3, pi=1/3, f=NULL, u=NULL

alpha <- 0.025
beta <- 0.2
delta <- 0.25  	#standardized overall treatment effect
s <- 3
pi <- 1/s

## Note: Length of f and u should be the same as s; sum(f)=1; sum(f*u)=1. 
f <- c(1/3, 1/3, 1/3)
u <- c(1, 1, 1)

## Can set the random seed to get identical result every time
set.seed(5000)
prob.def1(alpha, beta, delta, s, f, u, pi)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
