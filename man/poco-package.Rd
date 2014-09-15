\name{poco}
\alias{poco}
\alias{poco}
\docType{package}
\title{
Power and smaple size calculations for consistency assessment
}
\description{
In this package provides power and sample size calculations for the analysis of treatment-by-subgroup interactions for continuous and binary response variables:

- calculate the power to detect an interaction effect (also qualitative interactions)

- calculate sample size to detect an interaction effect (also qualitative interactions)

- calculate power for consistency assessment

- calculate sample size for consistency assessment

}
\details{
\tabular{ll}{
Package: \tab poco\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2014-08-28\cr
License: \tab GPL\cr
}
}
\author{
Andreas Kitsche <kitsche@biostat.uni-hannover.de>
}
\references{
Kitsche A, Power and sample size computations for simultaneous consistency assessment of treatment effects [manuscript in preparation]

Kitsche A, Hothorn LA. Testing for qualitative interaction using ratios of treatment differences. Statistics in Medicine (2014); 13(9):1477-1489

Kitsche A. Detecting qualitative interactions in clinical trials with binary responses. Pharmaceutical statistics (2014); early view
}

\keyword{ package }
\seealso{
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
            alpha = 0.05,
            type="anypair")
}
