\name{HellCor}
\alias{HellCor}

\title{
The Hellinger Correlation
}
\description{
Empirical value of the Hellinger correlation between two continuous random variables X and Y.}
\usage{
HellCor(x, y, Kmax = 20L, Lmax = 20L, K = 0L, L = 0L,
        alpha = 6.0, pval.comp = FALSE, conf.level = NULL,
        B1 = 200, B2 = 200, C.version = TRUE)
}

\arguments{

 \item{x}{Sample of X-values.}

  \item{y}{Sample of Y-values.}

  \item{Kmax}{Maximal number of terms to consider in the expansion of
  the square root of the copula density in the basis of Legendre polynomials.}

  \item{Lmax}{Maximal number of terms to consider in the expansion of
  the square root of the copula density in the basis of Legendre polynomials.}

  \item{K}{If K and L are not equal to zero, and Kmax = Lmax = 0, the
  number of terms is set to K in the expansion of
  the square root of the copula density in the basis of Legendre polynomials.}

  \item{L}{If K and L are not equal to zero, and Kmax = Lmax = 0, the
  number of terms is set to L in the expansion of
  the square root of the copula density in the basis of Legendre polynomials.}

  \item{alpha}{Parameter of the Beta(alpha,alpha) distribution
    used to fix boundary issues through transformation when estimating the Hellinger correlation.}

  \item{pval.comp}{Logical. Should we compute the p-value?}

  \item{conf.level}{If not \code{NULL}, a number in (0, 1) that sets the
  confidence level for the computed confidence interval.}

  \item{B1}{Numeric. Number of Bootstrap replicates for the outer loop
  in the double bootstrap procedure used to obtain the confidence
  interval. Only used if \code{conf.level} is not \code{NULL}. For some
  examples, one might need to increase this value to 1,000.}

  \item{B2}{Numeric. Number of Bootstrap replicates for the inner loop
  in the double bootstrap procedure used to obtain the confidence
  interval. Only used if \code{conf.level} is not \code{NULL}.}

\item{C.version}{Logical. If \code{FALSE}, the R version is used instead
  of the C version. Not really useful to change to \code{FALSE} unless
  one wants to understand the algorithm used.}
}

\details{
When \code{Kmax = Lmax = K = L = 0}, the value returned in \code{Hcor} is the
unnormalized version of the empirical Hellinger correlation.
  }

\value{
  List with the following components:
  \item{Hcor}{Value of the empirical Hellinger correlation.}
  \item{p.value}{The p-value associated to the null hypothesis that the Hellinger correlation is
  equal to 0 (computed using a Monte-Carlo simulation). Will be
  \code{NA} if \code{pval.comp} is not set to \code{TRUE}.}
  \item{conf.int}{Confidence interval for the population Hellinger
  correlation (computed using a double Bootstrap). Will be
  \code{NA} if \code{conf.level} is not set to a value in (0,1).}
  \item{Khat}{Value of \eqn{K} selected by cross-validation, if
  computed, otherwise \code{NA}.}
  \item{Lhat}{Value of \eqn{L} selected by cross-validation, if computed, otherwise \code{NA}.}
}
\references{
Geenens G., Lafaye de Micheaux P., (2018). The Hellinger correlation. \emph{}  (), --.

}
\author{
Geenens G., Lafaye de Micheaux P.
}


\examples{

# We illustrate the application of our new measure using data extracted
#  from a study on the coral reef published in
#  [Graham, N.A.J. and Wilson, S.K. and Carr, P. and Hoey, A.S. and
#   Jennings, S. and MacNeil, M.A. (2018), Seabirds enhance coral reef
#   productivity and functioning in the absence of invasive rats, Nature,
#   559, 250--253.]. 
# The two variables we consider are the number of fishes and the nitrogen 
# input by seabirds per hectare recorded on n = 12 islands in the Chagos
# Archipelago. Nitrogen input is an indirect measure of the abundance of
# seabirds. It is worthwhile to notice that nitrogen is absorbed by algae 
# and that herbivorous fishes eat these algae. Fishes and birds living in 
# two different worlds, it might seem odd to suspect a dependence between 
# these two variables. Since our measure is tailored to capture dependence 
# induced by a lurking variable, we can suspect the existence of such a
# hidden variable. This is indeed what researchers in (Graham et al., 2018) 
# were able to show by finding that the presence and abundance of rats on 
# an island had dramatic effects on the number of seabirds (the rats eating 
# their eggs) and consequently on the input of nitrogen in seabirds' guano. 
# In turn, this would diminish the abundance of algae and fishes eating
# these algae.

data(Chagos)
n <- nrow(Chagos)

par.save <- par()$mfrow
par(mfrow = c(1, 2))
plot(Chagos$Seabirds_ha, Chagos$Number_of_fishes, main = "Original data",
     xlab = "Density of seabirds", ylab = "Density of fishes")
plot(rank(Chagos$Seabirds_ha) / (n + 1), rank(Chagos$Number_of_fishes) /
     (n + 1), main = "Rank-Copula transformed data",
     xlab = "Density of seabirds", ylab = "Density of fishes")
par(mfrow = par.save)

set.seed(1)
# Empirical Hellinger correlation
HellCor(Chagos$Seabirds_ha, Chagos$Number_of_fishes, pval.comp = TRUE)
# Pearson correlation
cor.test(Chagos$Seabirds_ha, Chagos$Number_of_fishes)
# Distance correlation
dcor.test(Chagos$Seabirds_ha, Chagos$Number_of_fishes, R = 200)


set.seed(1)
# Empirical Hellinger correlation
HellCor(Chagos$kg_N_ha_yr, Chagos$Seabirds_ha, pval.comp = TRUE)
# Pearson correlation
cor.test(Chagos$kg_N_ha_yr, Chagos$Seabirds_ha)
# Distance correlation
dcor.test(Chagos$kg_N_ha_yr, Chagos$Seabirds_ha, R = 200)

set.seed(1)
# Empirical Hellinger correlation
HellCor(Chagos$kg_N_ha_yr, Chagos$Number_of_fishes, pval.comp = TRUE)
# Pearson correlation
cor.test(Chagos$kg_N_ha_yr, Chagos$Number_of_fishes)
# Distance correlation
dcor.test(Chagos$kg_N_ha_yr, Chagos$Number_of_fishes, R = 200)

t.test(Chagos$kg_N_ha_yr ~ Chagos$Treatment)

######################################################################
# Geenens G., Lafaye de Micheaux P., (2020). The Hellinger correlation
######################################################################
# Figure 5.2

\dontrun{
n <- 500

set.seed(1)
par(mfrow = c(3, 5))

XX <- .datagenW.corrected(n)
plot(XX, xlab = expression(X[1]), ylab = expression(X[2]))
etahat <- HellCor(XX[,1], XX[,2])$Hcor
title(main = substitute(W~~-~~hat(eta)==etahat, list(etahat =
round(etahat, 3))))

XX <- .datagenDiamond(n)
plot(XX, xlab = expression(X[1]), ylab = expression(X[2]))
etahat <- HellCor(XX[,1], XX[,2])$Hcor
title(main = substitute(Diamond~~-~~hat(eta)==etahat, list(etahat =
round(etahat, 3))))

XX <- .datagenParabola.corrected(n)
plot(XX, xlab = expression(X[1]), ylab = expression(X[2]))
etahat <- HellCor(XX[,1], XX[,2])$Hcor
title(main = substitute(Parabola~~-~~hat(eta)==etahat, list(etahat =
round(etahat, 3))))

XX <- .datagen2Parabolas.corrected(n)
plot(XX, xlab = expression(X[1]), ylab = expression(X[2]))
etahat<-HellCor(XX[,1], XX[,2])$Hcor
title(main = substitute(Two~~parabolae~~-~~hat(eta)==etahat, list(etahat
= round(etahat, 3))))

XX <- .datagenCircle.corrected(n)
plot(XX, xlab = expression(X[1]), ylab = expression(X[2]))
etahat <- HellCor(XX[,1], XX[,2])$Hcor
title(main = substitute(Circle~~-~~hat(eta)==etahat, list(etahat =
round(etahat, 3))))

XX <- .datagen4indclouds(n)
plot(XX, xlab = expression(X[1]), ylab = expression(X[2]))
etahat <- HellCor(XX[,1], XX[,2])$Hcor
title(main = substitute(4~~clouds~~-~~hat(eta)==etahat, list(etahat =
round(etahat, 3))))

XX <- .datagenCubic(n)
plot(XX, xlab = expression(X[1]), ylab = expression(X[2]))
etahat <- HellCor(XX[,1], XX[,2])$Hcor
title(main = substitute(Cubic~~-~~hat(eta)==etahat, list(etahat =
round(etahat, 3))))

XX <- .datagenSine(n)
plot(XX, xlab = expression(X[1]), ylab = expression(X[2]))
etahat <- HellCor(XX[,1], XX[,2])$Hcor
title(main = substitute(Sine~~-~~hat(eta)==etahat, list(etahat =
round(etahat, 3))))

XX <- .datagenWedge(n)
plot(XX, xlab = expression(X[1]), ylab = expression(X[2]))
etahat <- HellCor(XX[,1], XX[,2])$Hcor
title(main = substitute(Wedge~~-~~hat(eta)==etahat, list(etahat =
round(etahat, 3))))

XX <- .datagenCross(n)
plot(XX, xlab = expression(X[1]), ylab = expression(X[2]))
etahat <- HellCor(XX[,1], XX[,2])$Hcor
title(main = substitute(Cross~~-~~hat(eta)==etahat, list(etahat =
round(etahat, 3))))

XX <- .datagenSpiral(n, sigma = 0.05)
plot(XX, xlab = expression(X[1]), ylab = expression(X[2]))
etahat <- HellCor(XX[,1], XX[,2])$Hcor
title(main = substitute(Spiral~~-~~hat(eta)==etahat, list(etahat =
round(etahat, 3))))

XX <- .datagen4Circles(n)
plot(XX, xlab = expression(X[1]), ylab = expression(X[2]))
etahat <- HellCor(XX[,1], XX[,2])$Hcor
title(main = substitute(Circles~~-~~hat(eta)==etahat, list(etahat =
round(etahat, 3))))

XX <- .datagenHeavisine(n)
plot(XX, xlab = expression(X[1]), ylab = expression(X[2]))
etahat <- HellCor(XX[,1], XX[,2])$Hcor
title(main = substitute(Heavisine~~-~~hat(eta)==etahat, list(etahat =
round(etahat, 3))))

XX <- .datagenDoppler(n)
plot(XX, xlab = expression(X[1]), ylab = expression(X[2]))
etahat<-HellCor(XX[,1], XX[,2])$Hcor
title(main = substitute(Doppler~~-~~hat(eta)==etahat, list(etahat =
round(etahat, 3))))

XX <- .datagen5Clouds(n)
plot(XX, xlab = expression(X[1]), ylab = expression(X[2]))
etahat <- HellCor(XX[,1], XX[,2])$Hcor
title(main = substitute(5~~clouds~~-~~hat(eta)==etahat, list(etahat =
round(etahat, 3))))
}

############################
# Birth rate vs Death rate #
############################

data(worlddemographics)
x <- wdemographics$Death.Rate.Pop
y <- wdemographics$Birth.Rate.Pop
plot(x, y, xlab = "DEATHS/1,000 POPULATION", ylab =
"BIRTHS/1,000 POPULATION", main =
"Birth rate vs Death rate for 229 countries and territories (in 2020)", col = "orangered", pch = 20)
text(x, y, labels = wdemographics$Country, pos = 3, cex = 0.4, offset = 0.2)
\dontrun{
 HellCor(x, y, pval.comp = TRUE)
 cor.test(x, y) # Pearson
 dcor.test(x, y, R = 200)
}
}
\keyword{package}
