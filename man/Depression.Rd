\name{Depression}
\alias{Depression}
\docType{data}
\title{
Multi-centre depression trial
}
\description{
In themulti-centre depression trial, two groups of patients, one treatment and one placebo group, were compared. The primary endpoint was the change from the baseline to the end of the 9-week acute treatment phase in the 17-item Hamilton depression rating scale total score (HAMD17 score). The scores range from -2 to 28, and therefore, we assume that this end-point is approximately normally distributed. The experiment was conducted at five centres.
}
\usage{data("Depression")}
\format{
  A data frame with 100 observations on the following 3 variables.
  \describe{
    \item{\code{Centre}}{a numeric vector,  determining the centre}
    \item{\code{Group}}{a factor with levels \code{D} \code{P}, determining the treatment group}
    \item{\code{Score}}{a numeric vector, determining the HAMD17 depression score}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
Dmitrienko A, Molenberghs G, Chuang-Stein C, Offen W. Analysis of Clinical Trials using SAS: A Practical Guide.SAS Institute Inc: Cary, NC, 2005.
}
\references{
Kitsche, A., Hothorn, L.A. Testing for qualitative interaction using ratios of treatment differences. Statistics in Medicine 2014; 13(9):1477-1489.

Kitsche, A. Power and sample size computations for simultaneous consistency assessment of treatment effects [manuscript in preparation]
}
\examples{
data(Depression)
## maybe str(Depression) ; plot(Depression) ...
}
\keyword{datasets}
