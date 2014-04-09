\name{nNonInf2}
\alias{nNonInf2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculate the sample size associated with simultaneous tests for non-inferiority
}
\description{
This function calculates the sample size associated with simultaneous tests for non-inferiority for the ratio-based inference when large response values indicate better treatment effects based usiong the complete power definition.
}
\usage{
nNonInf2(r, h, psi, comp.power, k0, theta.star, alpha, n.start)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
number of comparisons
}
  \item{h}{
interger value for an a-priori known upper bound for the m true alternatives (if there is no prior information abut m set g=r)
}
  \item{psi}{
relative non-inferiority margin
}
  \item{comp.power}{
minimal power
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
  \item{n.start}{
starting value for the sample size
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{sample.size}{Description of 'comp1'}
\item{power}{Description of 'comp2'}
}
\references{
Dilba, G., Bretz, F., Hothorn, L.A., Guiard, V. Power and sample size computations in simultaneous tests for non-inferiority based on relative margins
Statistics in Medicine 25: 1131-1147 (2006)
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#Example from Dilba et al. (2006):
#r =3, g=1, psi=0.70, 1-beta=0.80, k0=0.5, theta_star=0.95 and alpha==0.05.
nNonInf(r=3, h=3, psi=0.7, comp.power=0.8, k0=0.5, theta.star=0.95,
        alpha=0.05, n.start=2)}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
