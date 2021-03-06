\name{alphaPowCon}
\alias{alphaPowCon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculate the alpha level for multiple ratios of treatment differences
}
\description{
This function calculates the alpha level  associated with simltaneous tests for the ratios of treatment differences with a user defined inconsistency margin for the assessment of treatment-by-subgroup interactions. The alpha level is determined iteratively until the power condition is satisfied.
}
\usage{
alphaPowCon(power, n, mu, sd, n.sub = 2, 
            TreatMat = "Tukey", SubMat = "GrandMean", 
            thetas = 1, alternative = c("two.sided", "less", "greater"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{power}{
any pair power
}
  \item{n}{
a single numeric value, number samples per treatment-by-subgroup combination
}
  \item{mu}{
numeric vector of cell means from the cell means model corresponding to a two-way layout
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
  \item{thetas}{
a single numeric value or a numeric vector of length, inconsistency margin(s)
}
  \item{alternative}{
a single character string, specifying the direction of the alternative hypothesis, one of "two.sided","less","greater"
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
alphaPowCon return a list containing the following components:
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
}
\references{
Kitsche, A., Hothorn, L.A. Testing for qualitative interaction using ratios of treatment differences. Statistics in Medicine 2014; 13(9):1477-1489.
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
#analysis of the consistency in a reproduced pivotal study
alphaPowCon(power=0.8,
       mu=c(10,15,10,14), 
       sd=sqrt(5),
       n = 30,
       n.sub=2,
       TreatMat= "Tukey", 
       SubMat = "Tukey",
       thetas=0.5, 
       alternative="greater")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
