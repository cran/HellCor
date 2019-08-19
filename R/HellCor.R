HellCor <- function(x, y, KLmax = 20L, alpha = 6.0, pvalue = FALSE) {
  
  if (!is.integer(KLmax)) stop("'KLmax' should be an integer.")
  if (length(x) != length(y)) stop("'x' and 'y' should have the same length.")
  
  res <- .C("hellcorC", as.double(cbind(x, y)), as.integer(2 * length(x)), statistic = as.double(0.0), as.integer(pvalue), pvalue = as.double(0.0), as.integer(KLmax), as.double(alpha), PACKAGE = "HellCor")
  
  if (pvalue) {
    n <- length(x)
    twon <- 2 * n
    y <- runif(n)
    pval <- mean(replicate(10000, .C("hellcorC", as.double(cbind(runif(n), y)), as.integer(twon), statistic = as.double(0.0), as.integer(0), as.double(0.0), as.integer(2), as.double(alpha), PACKAGE = "HellCor")$statistic) > res$statistic)
  } else {pval <- NA}
  
  return(list(Hcor = res$statistic, pvalue = pval))
  
}

.datagenW.corrected <- function(n) {
  x <- runif(n, -1, 1)
  u <- x + runif(n) / 3
  v <- 4 * ((x ^ 2 - 0.5) ^ 2 + runif(n) / 500)
  return (cbind(u, v))
}
.datagenDiamond <- function(n) {
  x <- runif(n, min = -1, max = 1)
  y <- runif(n, min = -1, max = 1)
  return (sqrt(2) * (cbind(x - y, x + y)) / 2)
}
.datagenParabola.corrected <- function(n) {
  x <- runif(n, -1, 1)
  y <- (x ^ 2 + runif(n)) / 2
  return (cbind(x, y))
}
.datagen2Parabolas.corrected <- function(n) {
  x <- runif(n, -1, 1)
  y <- (x ^ 2 + runif(n) / 2) * (sample(c(-1, 1), size = n, replace = TRUE))
  return (cbind(x, y))
}
.datagenCircle.corrected <- function(n) {
  x <- runif(n, -pi, pi)
  u <- sin(x) + rnorm(n) / 8
  v <- cos(x) + rnorm(n) / 8
  return (cbind(u, v))
}
.datagen4indclouds <- function(n) {
  dx <- rnorm(n) / 3
  dy <- rnorm(n) / 3
  cx <- sample(c(-1, 1), size = n, replace = TRUE)
  cy <- sample(c(-1, 1), size = n, replace = TRUE)
  u <- cx + dx
  v <- cy + dy 
  return (cbind(u, v))
}
.datagenCubic <- function(n) {
  x <- runif(n, -1.3, 1.1)
  x2 <- x + runif(n, -1.3, 1.1) / 2
  y <- 4 * x2 ^ 3 + x2 ^ 2 - 4 * x2 
  return (cbind(x, y))
}
.datagenSine <- function(n) {
  x <- runif(n, 0, 1)
  y <- sin(8 * pi * x) + runif(n, -1, 1) / 3
  return (cbind(x, y))
}
.datagenWedge <- function(n) {
  x <- runif(n, 0, 1)
  y <- x * (sample(c(-1, 1), size = n, replace = TRUE)) + runif(n, -1, 1) / 2
  return (cbind(x, y))
}
.datagenCross <- function(n) {
  x <- runif(n, 0, 1)
  s <- sample(c(-1, 1), size = n, replace = TRUE)
  y <- 2 * s * x - (s - 1) + runif(n, -1, 1) / 2
  return (cbind(x, y))
}
.datagenSpiral <- function(n,sigma=1/8) {
  theta <- runif(n, 0, 6 * pi)
  r <- (1/3) * theta
  x <- r * cos(theta) + rnorm(n)*sigma
  y <- r * sin(theta) + rnorm(n)*sigma
  return (cbind(x, y))
}
.datagen4Circles <- function(n) {
  theta <- runif(n, -pi, pi)
  s1 <- sample(c(-1, 1), size = n, replace = TRUE)
  s2 <- sample(c(-1, 1), size = n, replace = TRUE)
  y.center <- s1 * 1
  x.center <- s2 * 1
  radius <- 1.4
  y <- radius * sin(theta) + y.center
  x <- radius * cos(theta) + x.center
  return (cbind(x, y))
}
.datagenHeavysine <- function(n) {
  x <- runif(n, 0, 1)
  y <- (4 * sin(4 * pi * x) - sign(x - 0.3) - sign (0.72 - x))/3 + runif(n, -1, 1) / 2
  return (cbind(x, y))
}
.datagenDoppler <- function(n) {
  x <- runif(n, 0, 1)
  y <- sqrt(x * (1 - x)) * sin(2.1 * pi / (x + 0.05)) + runif(n, -1, 1) / 2
  return (cbind(x, y))
}
.datagen5Clouds <- function(n) {
  dx <- rnorm(n) / 5
  dy <- rnorm(n) / 5
  cx <- sample(c(-1, 1), size = n, replace = TRUE)
  cy <- sample(c(-1, 1), size = n, replace = TRUE)
  s <- sample(c(0, 1, 1, 1, 1), size = n, replace = TRUE)
  u <- s * cx + dx
  v <- s * cy + dy 
  return (cbind(u, v))
}
