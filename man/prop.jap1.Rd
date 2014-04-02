\name{prop.jap1}
\alias{prop.jap1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculate proportion of Japanese patients in a multi-regional trial
}
\description{
This function calculates the proportion of Japanese patients in a multi-regional clinical trial. The methodlogy is in accordence to Method 1 of the MHLW Guidance. It assumes that \eqn{\delta}(Japan)/\eqn{\delta}(ALL) > \eqn{\omega} and \eqn{\delta}(ALL) > 0.
}
\usage{
prop.jap1(omega = 0.5, beta = 0.8, alpha = 0.025, gamma = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{omega}{
pre-specified consistency threshold
}
  \item{beta}{
type II error
}
  \item{alpha}{
type I error
}
  \item{gamma}{
Pr(\eqn{\delta}(Japan)>\eqn{\omega}\eqn{\delta}(ALL)=1-\eqn{\gamma}, 0<\eqn{\gamma}<1
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{Proportion}{proportion of Japanese patients}
}
\references{
Ikeda, K. ,Bretz, F. Sample size and proportion of Japanese patients in multi-regional trials
Pharmaceutical Statistics 9: 207-216 (2010)
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
#According to Ikeda and Bretz the proportions of japanese patients for
#alpha=0.025, omega=0.5, gamma=0.2 and beta=0.2
powB(delta=0.125, omega=0.5 , na=100 ,nb=100 , sigma=1 ,alpha=0.025, p=0.1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
