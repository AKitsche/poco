\name{PowNonInf}
\alias{PowNonInf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculate the any-pair power at the least favourable configuration (LFC)
associated with simultaneous tests for non-inferiority
}
\description{
This function calculates the minimal power (any-pair power) at the least favourable configuration (LFC) associated with simultaneous tests for non-inferiority for the ratio-based inference when large response values indicate better treatment effects
}
\usage{
PowNonInf(r, g, psi, min.power, k0, theta.star, alpha, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
number of comparisons (number of treatment groups)
}
  \item{g}{
interger value for an a-priori known lower bound for the m true alternatives (if there is no prior information about m set g=1)
}
  \item{psi}{
relative non-inferiority margin
}
  \item{k0}{
coefficient of variation of the control group
}
  \item{theta.star}{
smallest clinically irrelevant percentage of the control to be detected
}
  \item{alpha}{
familywise type I error rate
}
  \item{n}{
per group sample size
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
Dilba, G., Bretz, F., Hothorn, L.A., Guiard, V. Power and sample size computations in simultaneous tests for non-inferiority based on relative margins Statistics in Medicine 25: 1131-1147 (2006)
}
\author{
Dilba, G.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
PowNonInf(r=3, g=1, psi=0.8, k0=0.2, 
          theta.star=0.8, alpha=0.05, n=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
