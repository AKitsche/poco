\name{Trastuzumab2}
\alias{Trastuzumab2}
\docType{data}
\title{
Multi-regional clinical trial considering treatment of
trastuzumab (Herceptin, Roche) after adjuvant chemotherapy in HER2-positive breast cancer
}
\description{
The international, randomized clinical trial compared two treatment groups, one or two years of trastuzumab being given every three weeks, with observation in patients with HER2-positive breast cancer. The primary endpoint was disease-free survival, defined as time from randomization to the first occurrence of an event, where an event is given as recurrence of breast cancer, contralateral breast cancer, second non-breast malignant disease, or death. The endpoint under consideration is the occurrence of an event because the disease-free survival times were not available. Five regions were pre-specified, namely Central and South America, Eastern Europe, Asia Pacific including Japan and Others (Others includeWestern and Northern Europe, Canada, South Africa, Australia, New Zealand). In the report of the interim analysis, only the results of the treatment group with trastuzumab treatment for one year and the observation group were presented.
}
\usage{data("Trastuzumab2")}
\format{
  A data frame with 3473 observations on the following 4 variables.
  \describe{
    \item{\code{Region}}{a factor with levels \code{Asia Pacific, Japan} \code{Central and South America} \code{Eastern Europe} \code{Japan only} \code{Others}, determining the countries/regions}
    \item{\code{Treatment}}{a factor with levels \code{Observation} \code{Trastuzumab}}
    \item{\code{Success}}{a numeric vector, determining the presence of a desease}
    \item{\code{RegionTreat}}{a factor with 10 levels, determing each treatment-by-region group }
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Romond, E. H., Perez, E. A., Bryant, J., Suman, V. J., Geyer Jr., C. E., Davidson, N. E., Tan-Chiu, E., Martino, S., Paik, S.,Kaufman, P. A., Swain, S. M., Pisansky, T. M., Fehrenbacher, L., Kutteh, L. A., Vogel, V. G., Visscher, D. W., Yothers, G., Jenkins, R. B., Brown, A. M., Dakhil, S. R., Mamounas, E. P., Lingle,W. L., Klein, P. M., Ingle, J. N., andWolmark, N. (2005). Trastuzumab plus adjuvant chemotherapy for operable HER2-positive breast cancer. New England Journal of Medicine, 353(16):1673-1684.
}
\examples{
data(Trastuzumab2)
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
            alpha = 0.05,
            type="anypair")
}
\keyword{datasets}
