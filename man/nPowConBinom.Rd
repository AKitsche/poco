\name{nPowConBinom}
\alias{nPowConBinom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculate the sample size for multiple ratios of treatment differences for binomial response variables
}
\description{
This function calculates the sample size associated with simltaneous tests for the ratios of treatment differences with a user defined inconsistency margin for the assessment of treatment-by-subgroup interactions for binomial response variables. The required size n is determined iteratively by starting with a given sample size and search until the power condition is satisfied.
}
\usage{
nPowConBinom(min.power, p, n.sub = 2, 
             TreatMat = "Tukey", SubMat = "GrandMean", 
             rhs = 1, alternative = c("two.sided", "less", "greater"), 
             alpha = 0.05,type="anypair","allpair","global")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{min.power}{
any pair power
}
  \item{p}{
vector of success probabilities for each treatment-by-subset combination
}
  \item{n.sub}{
a single numeric value, number of levels of the subgrouping factor; if n.sub=1 the argument SubMat is ignored
}
  \item{TreatMat}{
a single character string specifying the type of contrast for the treatment factor; naming one of the contrast types available in contrMat(multcomp)
}
  \item{SubMat}{
a single character string specifying the type of contrast for the subgrouping factor; naming one of the contrast types available in contrMat(multcomp); if n.sub=1 the argument SubMat is ignored
}
  \item{rhs}{
a single numeric value or a numeric vector of length, inconsistency margin(s)
}
  \item{alternative}{
a single character string, specifying the direction of the alternative hypothesis, one of "two.sided","less","greater"
}
  \item{alpha}{
a single numeric value, familywise type I error to be controlled
}
  \item{type}{
  a single character string, naming the type of rejection probability to be computed;  "anypair" for the rejection probability that at least one of the alternative hypothses is true, "allpair" for the probability that all elementary alternatives are rejected,
  "global" for the global rejection probability.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
nPowCon return a list containing the following components:
\item{power}{any pair power}
\item{n}{number of samples per treatment-by-subgroup combination}
\item{NonCentrPar}{Non-centrality parameter}
\item{crit}{quantile from the multivariate t-distribution}
\item{alternative}{direction of the alternative hypothesis}
\item{CorrMat}{correlation matrix of the test statistics}
\item{CMat}{Numeration interaction contrast matrix}
\item{DMat}{Denominator interaction contrast matrix}
\item{thetas}{inconsistency margin}
\item{alpha}{global significance level}
\item{n.sub}{number of levels of the subgrouping factor}
\item{TreatMat}{type of contrast for the treatment factor}
\item{SubMat}{type of contrast for the subgrouping factor}
\item{type}{power definition used}
}
\references{
Kitsche A. Detecting qualitative interactions in clinical trials with binary responses. Pharmaceutical statistics (2014); early view

Kitsche A. Power and sample size computations for simultaneous consistency assessment of treatment effects [manuscript in preparation]
}
\author{
Andreas Kitsche
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
nPowConBinom(min.power=0.8, 
             p=c(0.2,0.3,0.2,0.8), 
             n.sub=2, 
             TreatMat = "Tukey", 
             SubMat = "GrandMean", 
             rhs = 1, 
             alternative = "two.sided",
             alpha = 0.05,
             tyoe="anypair")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
