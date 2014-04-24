\name{nPowCon}
\alias{nPowCon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculate the sample size for multiple ratios of treatment differences
}
\description{
This function calculates the sample sizeassociated with simltaneous tests for the ratios of treatment differences with a user defined inconsistency margin for the assessment of treatment-by-subgroup interactions. The required size n is determined iteratively by starting with a given sample size and search until the power condition is satisfied.
}
\usage{
nPowCon(min.power, n.start=1, mu, sd, n.sub = 2, 
        TreatMat = "Tukey", SubMat = "GrandMean", 
        thetas = 1, alternative = c("two.sided", "less", "greater"), 
        alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{min.power}{
any pair power
}
  \item{n.start}{
a single interger value, starting sample size
}
  \item{mu}{
vector of cell means from the cell means model corresponding to a two-way layout
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
#analysis of the consistency in a reproduced pivotal study
nPowCon(min.power=0.8,
        n.start=10,
        mu=c(10,15,10,14), 
        sd=sqrt(5),
        n.sub=2,
        TreatMat= "Tukey", 
        SubMat = "Tukey",
        thetas=0.5, 
        alpha=0.05, 
        alternative="greater")
        
#analysis of treatment-by-subgroup interaction (e.g. in a multi-regional trial)
nPowCon(min.power=0.9,
       n.start=10,
       mu=c(10,15,10,15,10,15,10,15,10,11), 
       sd=sqrt(5),
       n.sub=5,
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
