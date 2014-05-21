\name{PowConBinom}
\alias{PowConBinom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculate any pair power for multiple ratios of treatment differences for binomial distributed endpoint
}
\description{
This function calculates the minimal power (any-pair power) associated with simltaneous tests for the ratios of treatment differences with a user defined inconsistency margin for the assessment of treatment-by-subgroup interactions.
}
\usage{
PowConBinom(p, n, sd, n.treat = 2, TreatMat = "Tukey", SubMat = "GrandMean",
       rhs = 1, alternative = c("two.sided", "less", "greater"), alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu}{
vector of cell means from the cell means model corresponding to a two-way layout
}
  \item{n}{
a single numeric value, number samples per treatment-by-subgroup combination
}
  \item{sd}{
a single numeric value, pooled standard deviation
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
a single numeric value or a numeric vector of length,  inconsistency margin(s)
}
  \item{alternative}{
a single character string, specifying the direction of the alternative hypothesis, one of "two.sided","less","greater"
}
  \item{alpha}{
a single numeric value, familywise type I error to be controlled
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
PowCon return a list containing the following components:
\item{power}{any pair power}
\item{n}{number of samples per treatment-by-subgroup combination}
\item{NonCentrPar}{Non-centrality parameter}
\item{crit}{quantile from the multivariate t-distribution}
\item{alternative}{direction of the alternative hypothesis}
\item{CorrMat}{correlation matrix of the test statistics}
\item{CMat}{Numeration interaction contrast matrix}
\item{DMat}{Denominator interaction contrast matrix}
\item{rhs}{inconsistency margin}
\item{alpha}{global significance level}
\item{n.sub}{number of levels of the subgrouping factor}
\item{TreatMat}{type of contrast for the treatment factor}
\item{SubMat}{type of contrast for the subgrouping factor}
}
\references{
Kitsche, A. Detecting qualitative interactions in clinical trials with binary responses. Pharmaceutical Statistics, under revision.
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
#Power calculation for consistency assessment of the 
#Metoprolol Controlled-Release Randomized Intervention Trial in Heart Failure

#generating the data set
data(MetoCRXL2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
