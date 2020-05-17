\name{wdemographics}
\alias{worlddemographics}
\docType{data}
\title{
Birth rates and death rates for 229 countries in 2020
}
\description{
This data set contains the birth rates and death rates of 229 countries
in 2020 ranked in decreasing order.
}
\usage{wdemographics}
\format{
  Data frame with 229 observations and 3 variables.
  \describe{
    \item{\code{Country}}{Character. Name of the country.}
    \item{\code{Birth.Rate.Pop}}{Numeric. Birth rate per 1,000 habitants.}
    \item{\code{Death.Rate.Pop}}{Numeric. Death rate per 1,000 habitants.}
  }
}
\references{
Source: The World Factbook of the Central Intelligence Agency of the
United States (\url{https://www.cia.gov/library/publications/the-world-factbook/geos/xx.html}).

See also Figure 1.1 in "A consistent test of independence between random vectors", G. Boglioni Beaulieu,
 Master thesis, Universite de Montreal, (2016).
}
\author{
Geenens G., Lafaye de Micheaux P.
}
\examples{
data(wdemographics)
}
\keyword{datasets}
