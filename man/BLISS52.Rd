\name{BLISS52}
\alias{BLISS52}
\docType{data}
\title{
A randomised, placebo-controlled, phase 3 trial compared belimumab
(trade name Benlysta) to placebo in patients with systemic lupus erythematosus (SLE)
}
\description{
BLISS-52 randomized patients from Asia, South America, and Eastern Europe; patients were to be treated for 52 weeks. The primary outcome was a binary outcome measured at 52 weeks with success defined by an SLE Responder Index (SRI); participating three regions in 13 countries in Latin America (Argentina, Brazil, Chile, Colombia, and Peru), Asia-Pacific (Australia, Hong Kong, India, Korea, Philippines, and Taiwan), and eastern Europe (Romania and Russia)
}
\usage{data("BLISS52")}
\format{
  A data frame with 9 observations on the following 4 variables.
  \describe{
    \item{\code{Region}}{a factor with levels \code{AsiaPacific}, \code{EasternEurope}, \code{LatinAmerica}}
    \item{\code{Treatment}}{a factor with levels \code{10mg} \code{1mg} \code{Placebo}}
    \item{\code{SampleSize}}{a numeric vector, determing the sample size}
    \item{\code{Events}}{a numeric vector, determining the number of events}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Efficacy and safety of belimumab in patients with active systemic lupus erythematosus: a randomised, placebo-controlled, phase 3 trial. Lancet 2011; 377(9767):721-731.
}
\examples{
data(BLISS52)
## maybe str(BLISS52) ; plot(BLISS52) ...
}
\keyword{datasets}
