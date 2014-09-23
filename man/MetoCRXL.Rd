\name{MetoCRXL}
\alias{MetoCRXL}
\docType{data}
\title{
Metoprolol Controlled-Release Randomized Intervention Trial in Heart Failure
}
\description{
Data set with the number of successes and failures for each region in the Metoprolol Controlled-Release Randomized Intervention Trial in Heart Failure.
}
\usage{data("MetoCRXL")}
\format{
  A data frame with 24 observations on the following 4 variables.
  \describe{
    \item{\code{Country}}{a factor with levels \code{Belgium} \code{Czech Republic} \code{Denmark/Finland} \code{Germany} \code{Hungary} \code{Iceland} \code{Norway} \code{Poland} \code{Sweden} \code{The Netherland/Switzerlnd} \code{UK} \code{USA}}
    \item{\code{Treatment}}{a factor with levels \code{Meto CR/XL} \code{Placebo}}
    \item{\code{SampleSize}}{a numeric vector, number of patients}
    \item{\code{Events}}{a numeric vector, number of successes}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
MERIT-HF Study Group. Effect of metoprolol CR/XL in chronic heart failure: Metoprolol CR/XL Randomised Intervention Trial In-congestive Heart Failure (MERIT-HF). Lancet (1999); 353(9169):2001-2007.
}
\examples{
data(MetoCRXL)
}
\keyword{datasets}
