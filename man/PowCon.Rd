\name{PowCon}
\alias{PowCon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculate any pair power for multiple ratios of treatment differences
}
\description{
This function calculates the minimal power (any-pair power) associated with simltaneous tests for the ratios of treatment differences with a user defined inconsistency margin for the assessment of treatment-by-subgroup interactions.
}
\usage{
PowCon(mu, n, sd, n.treat = 2, TreatMat = "Tukey", SubMat = "GrandMean",
       thetas = 1, alternative = c("two.sided", "less", "greater"), alpha = 0.05)
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
  \item{n.treat}{
a single numeric value, number of levels of the subgrouping factor
}
  \item{TreatMat}{
a single character string specifying the type of contrast for the treatment factor; naming one of the contrast types available in contrMat(multcomp)
}
  \item{SubMat}{
a single character string specifying the type of contrast for the subgrouping factor; naming one of the contrast types available in contrMat(multcomp)
}
  \item{thetas}{
a single numeric value, inconsistency margin
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
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
#according to Kitsche and Hothorn (2014) 
#analyis of a multi-centre clinical trial
PowCon(mu = c(10,8,10,12,10,12,10,12,10,12,10,12,10,12,10,12,10,12,10,12),
       n=30,
       sd=sqrt(5),
       n.treat=2,
       TreatMat="Tukey",
       SubMat="GrandMean",
       thetas=0,
       alternative="less",
       alpha=0.05)
       
#analysis of the consistency in a reproduced pivotal study
PowCon(mu=c(10,15,10,14), 
       sd=sqrt(5),
       n = 30,
       TreatMat= "Tukey", 
       SubMat = "Tukey",
       thetas=0.5, 
       alpha=0.05, 
       alternative="greater")
       
#analysis of treatment-by-subgroup interaction (e.g. in a multi-regional trial)
PowCon(mu=c(10,15,10,15,10,15,10,15,10,11), 
       sd=sqrt(5),
       n = 10,
       TreatMat= "Tukey", 
       SubMat = "GrandMean",
       thetas=1, 
       alpha=0.05, 
       alternative="two.sided")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
