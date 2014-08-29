\name{MetoCRXL2}
\alias{MetoCRXL2}
\docType{data}
\title{
Metoprolol Controlled-Release Randomized Intervention Trial in Heart Failure
}
\description{
The large-scale randomised, double-blind, placebo-controlled trial was conducted to investigate the treatment effect of adding once-daily doses of metoprolol controlled-release/extended- release (MetoCR/XL) to the optimum standard therapy in terms of lowering mortality in patients with symptomatic heart failure. A total number of 3991 patients were randomised into the placebo or the Meto CR/XL group in 14 countries. The data from Finland were combined with the data from Denmark, and the data from the Netherlands were combined with the data from Switzerland because no event was observed in the Meto CR/XL group in Finland and Switzerland. From Table II, a decreasing overall treatment effect is observable, whereas in two regions, Iceland and USA, the treatment effect increases.
}
\usage{data("MetoCRXL2")}
\format{
  A data frame with 3991 observations on the following 4 variables.
  \describe{
    \item{\code{Region}}{a factor with levels 12 levels determining the countries/regions}
    \item{\code{Treatment}}{a factor with levels \code{Meto CR/XL} and \code{Placebo}}
    \item{\code{Success}}{a numeric vector, determining the presence of a desease}
    \item{\code{RegionTreat}}{a factor with 24 levels, determing each treatment-by-region group}
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
data(MetoCRXL2)
library(MCPAN)
MCPAN_Est <- binomest(Success ~ RegionTreat,
                      data=MetoCRXL2, 
                      success="1", 
                      method="Wald")
MCPAN_Est$estp
MCPAN_Est$n
#power calculations
PowConBinom(p=MCPAN_Est$estp, 
            n=MCPAN_Est$n, 
            n.sub = 12, 
            TreatMat = "Tukey", 
            SubMat = "GrandMean",
            rhs = 0, 
            alternative = "less", 
            alpha = 0.05)
}
\keyword{datasets}
