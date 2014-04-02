\name{f11.def2}
\alias{f11.def2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculate proportion of study patients for one region
}
\description{
Purpose: To calculate the minimum proportion of study patients for one region (f1) to have an specified chance of showing consistency by Definition 1.

Definition 2: Observing region effects that exceed a pre-specified effect size.  

Under the restriction f1 < f2 = f3 = f4
}
\usage{
f11.def2(alpha = 0.025, beta = 0.01, delta = 0.005, sigma = 0.013, 
         u = c(1, 1, 1, 1), b = 0.005/4, set.pow = 0.8)
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
  \item{sigma}{
varianz
}
  \item{u}{
Vector of \eqn{(f_{1},f_{2},\ldots,f_{s})}, where \eqn{f_{i}} is the proportion of patients within region \eqn{i} and \eqn{\sum_{i=1}^{s}f_{i}=1}
}
  \item{b}{
pre-specified effect size
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
#alpha=0.025, beta=0.01, delta=0.005, sigma=0.013,
#u=c(1,1,1,1), b=0.005/4, set.pow=0.8

alpha <- 0.025
beta <- 0.01
delta <- 0.005
sigma <- 0.013
b <- 0.005/4
u <- rep(1, 4)
set.pow <- 0.8    

#For some higher powers, there may not be a solution, 
#may need to increase total sample size

f11.def2(alpha, beta, delta, sigma, u, b, set.pow)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
