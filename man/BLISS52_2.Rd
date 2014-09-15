\name{BLISS52_2}
\alias{BLISS52_2}
\docType{data}
\title{
A randomised, placebo-controlled, phase 3 trial compared belimumab
(trade name Benlysta) to placebo in patients with systemic lupus erythematosus (SLE)
}
\description{
BLISS-52 randomized patients from Asia, South America, and Eastern Europe; patients were to be treated for 52 weeks. The primary outcome was a binary outcome measured at 52 weeks with success defined by an SLE Responder Index (SRI); participating three regions in 13 countries in Latin America (Argentina, Brazil, Chile, Colombia, and Peru), Asia-Pacific (Australia, Hong Kong, India, Korea, Philippines, and Taiwan), and eastern Europe (Romania and Russia)
}
\usage{data("BLISS52_2")}
\format{
  A data frame with 865 observations on the following 4 variables.
  \describe{
    \item{\code{Region}}{a factor with levels \code{AsiaPacific}, \code{EasternEurope}, \code{LatinAmerica}}
    \item{\code{Treatment}}{a factor with levels \code{10mg} \code{1mg} \code{Placebo}}
    \item{\code{Success}}{a numeric vector, indicating the binary response variable}
    \item{\code{RegionTreat}}{a factor with levels \code{AsiaPacific:10mg} \code{AsiaPacific:1mg} \code{AsiaPacific:Placebo} \code{EasternEurope:10mg} \code{EasternEurope:1mg} \code{EasternEurope:Placebo} \code{LatinAmerica:10mg} \code{LatinAmerica:1mg} \code{LatinAmerica:Placebo}}
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
data(BLISS52_2)
BLISS52_2$Treatment <- factor(BLISS52_2$Treatment, levels=c("Placebo", "1mg", "10mg"))
library(MCPAN)
MCPAN_Est <- binomest(Success ~ RegionTreat,
                      data=BLISS52_2, 
                      success="1", 
                      method="Wald")
MCPAN_Est$estp
MCPAN_Est$n
#power calculations
PowConBinom(p=MCPAN_Est$estp, 
            n=MCPAN_Est$n, 
            n.sub = 3, 
            TreatMat = "Dunnett", 
            SubMat = "GrandMean",
            rhs = 0.5, 
            alternative = "less", 
            alpha = 0.05,
            type="anypair")
}
\keyword{datasets}
