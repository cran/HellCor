\name{Chagos}
\alias{Chagos}
\docType{data}
\title{
Coral fish, seabirds and reef productivity study by Graham et al
}
\description{
This data set was first studied by Graham etal (2018). It contains several
variables describing the fish and seabirds populations on twelve islands
of the Chagos Archipelago (British Indian Ocean Territory).
}
\usage{data(Chagos)}
\format{
  A data frame with 12 observations on the following 7 variables.
  \describe{
    \item{\code{Atoll}}{Factor w 3 levels. Atoll containing several islands.}
    \item{\code{Island}}{Factor w 12 levels. Name of the island.}
    \item{\code{Treatment}}{Factor w/ 2 levels. Presence or absence of
      rats on the island.}
    \item{\code{Seabirds_ha}}{numeric vector of the number of birds per hectare.}
    \item{\code{kg_N_ha_yr}}{numeric vector of nitrogen input (in kg) by
      seabirds per hectare over a one year period.}
    \item{\code{Number_of_fishes}}{integer vector of the number of
      fishes recorded during underwater visual surveys along the
reef crest of each island on the lagoonal side of each atoll.}
    \item{\code{biomass}}{numeric vector. Fish counts were converted into biomass using published length-weight
relationships from FishBase.}
  }
}
\references{
Graham N. A. J., Wilson S. K., Carr P., Hoey A. S., Jennings S., MacNeil
M. A. (2018). Seabirds enhance coral reef productivity and
functioning in the absence of invasive rats. \emph{Nature} 559, 250--253.
}
\author{
Geenens G., Lafaye de Micheaux P.
}
\examples{
data(Chagos)
}
\keyword{datasets}
