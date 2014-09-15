\name{PURSUIT2}
\alias{PURSUIT2}
\docType{data}
\title{
The Platelet Glycoprotein IIb/IIIa in Unstable Angina: Receptor Suppression Using Integrilin Therapy (PURSUIT) study 
}
\description{
The study was an international multi-center double-blind  #randomized placebo-controlled study evaluating the 180/2.0 regimen of INTEGRILIN in 10,948 patients. INTEGRILIN was added to aspirin and heparin and continued until hospital discharge or for up to 72 hours (or up to 96 hours if a coronary intervention was performed towards the end of the 72 hour period). The data set contains the information from four regions: eastern europe, western europe, north america and latin america.
}
\usage{data("PURSUIT2")}
\format{
  A data frame with 9461 observations on the following 4 variables.
  \describe{
    \item{\code{Region}}{a factor with levels \code{EE}, \code{LA}, \code{WE}, \code{NAm}, determing the region: EE - Eastern Europe, LA - Latin America, WE - Western Europe, NA - North America}
    \item{\code{Treat}}{a factor with levels \code{High} \code{Placebo}, determining the treatment group}
    \item{\code{X}}{a numeric vector, indicating the binary response variable}
    \item{\code{RegionTreat}}{a factor with levels \code{EE:High} \code{EE:Placebo} \code{LA:High} \code{LA:Placebo} \code{NAm:High} \code{NAm:Placebo} \code{WE:High} \code{WE:Placebo}}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
http://www.stat.ncsu.edu/courses/sibs/datasets/pursuit/
}
\references{
The PURSUIT Trial Investigators. Inhibition of platelet glycoprotein IIb/IIIa with eptifibatide in patients with acute coronary syndromes without persistent ST-segment elevation. N Engl J Med. 1998; 339:436-443.
}
\examples{
data(PURSUIT2)
library(MCPAN)
MCPAN_Est <- binomest(X ~ RegionTreat,
                      data=PURSUIT2, 
                      success="1", 
                      method="Wald")
MCPAN_Est$estp
MCPAN_Est$n
#power calculations
PowConBinom(p=MCPAN_Est$estp, 
            n=MCPAN_Est$n, 
            n.sub = 4, 
            TreatMat = "Tukey", 
            SubMat = "GrandMean",
            rhs = 0.5, 
            alternative = "less", 
            alpha = 0.05,
            type="anypair")
}
\keyword{datasets}
