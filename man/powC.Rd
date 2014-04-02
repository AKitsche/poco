\name{powC}
\alias{powC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculate the probability of showing a consistent result between all and Japanese patients
}
\description{
Function to calculate the joint probability of detecting a significant difference between the study drug and placebo in all patients and showing a consistent result between all and Japanese patients
}
\usage{
powC(delta = 0.125, omega = 0.5, sigma = 1, 
     alpha = 0.025, beta = 0.2, p = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
powC(delta = 0.125, omega = 0.5, na = 100, nb = 100, 
     sigma = 1, alpha = 0.025, p = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{
overall treatment effect
}
  \item{omega}{
consistency margin
}
  \item{sigma}{
square root of the common variance
}
  \item{alpha}{
type I error rate
}
  \item{beta}{
type II error rate  
}
  \item{p}{
proportion of Japanese patients
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{etaC}{probability of showing a consistent result between all and Japanese patients}
}
\references{
Ikeda, K. ,Bretz, F. Sample size and proportion of Japanese patients in multi-regional trials
Pharmaceutical Statistics 9: 207-216 (2010)
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
#According to Ikeda and Bretz the proportions of japanese patients for
#delta=0.1, omega=0.5, sigma=1 ,alpha=0.025, p=0.1, beta=0.1
powC(delta=0.1, omega=0.5, sigma=1 ,alpha=0.025, p=0.1, beta=0.1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
