HellCor <- function(x, y = NULL, Kmax = 20L, Lmax = 20L, K = 0L, L = 0L, alpha = 6.0,
                    pval.comp = FALSE, conf.level = NULL, B1 = 200L, B2 = 200L,
                    C.version = TRUE) {

  if (!is.null(y)) {
      if (length(x) != length(y)) stop("'x' and 'y' should have the same length.")
  } else {
    if (!is.matrix(x)) stop("'x' should be a matrix when 'y' is not supplied.")
    if (ncol(x) != 2) stop("'x' should be a matrix with two columns.")
    y <- x[, 2]
    x <- x[, 1]
  }
  Kmax <- as.integer(Kmax)
  Lmax <- as.integer(Lmax)
  K <- as.integer(K)
  L <- as.integer(L)
  B1 <- as.integer(B1)
  B2 <- as.integer(B2)
  if (K < 0) stop("'K' should be positive.")
  if (L < 0) stop("'L' should be positive.")
  if (!is.null(conf.level)) {  
    if (length(conf.level) > 1) stop("'conf.level' should be of length 1.")
    if ((B1 <= 0) | (B2 <= 0)) stop("'B1' and 'B2' should be positive if 'conf.level' is used.")
    eps <- 100 * .Machine$double.eps
    if (any(conf.level < -eps | conf.level > 1 + eps)) stop("'conf.level' outside [0,1]")
  }
  CIeta <- NA
  pval <- NA
  
  if (C.version) {
    res <- .C("hellcorC", as.double(cbind(x, y)), as.integer(2 * length(x)), statistic = as.double(0.0), as.integer(pval.comp),
              pvalue = as.double(0.0), Khat = as.integer(Kmax), Lhat = as.integer(Lmax), as.integer(K), as.integer(L), as.double(alpha),
              conf.level = as.double(0.0), as.integer(B1), as.integer(B2), CIleft = as.double(0.0), CIright = as.double(0.0), PACKAGE = "HellCor")
    if (!is.null(conf.level)) {
      res <- .C("hellcorC", as.double(cbind(x, y)), as.integer(2 * length(x)), statistic = as.double(0.0), as.integer(pval.comp),
                pvalue = as.double(0.0), Khat = as.integer(Kmax), Lhat = as.integer(Lmax), as.integer(K), as.integer(L), as.double(alpha),
                conf.level = as.double(conf.level), as.integer(B1), as.integer(B2), CIleft = as.double(0.0), CIright = as.double(0.0), PACKAGE = "HellCor")
      CIeta <- c(res$CIleft, res$CIright)
    }
  } else {
    n <- length(x)
    XX <- cbind(x, y)
    UU <- apply(XX, MARGIN = 2, FUN = rank) / (nrow(XX) + 1)
    TT <- apply(UU, MARGIN = 2, FUN = function(x) qbeta(x, alpha, alpha))
    weight  <- sqrt(dbeta(TT[, 1], alpha, alpha) * dbeta(TT[, 2], alpha, alpha))
    Dist <- dist(TT, diag = TRUE, upper = TRUE)
#    library("FNN")
    neighbors <- FNN::get.knn(TT, k = 2)
    R1 <- neighbors$nn.dist[, 1] # in row i: distance from point i to its closest neighbor
    R2 <- neighbors$nn.dist[, 2] # in row i: distance from point i to its second closest neighbor
    II <- neighbors$nn.index # in row i: index of closest and second closest neighbor to point i
    R1 <- R1 * weight
    R2 <- R2 * weight
#    library("orthopolynom")
    tmp <- orthopolynom::legendre.polynomials(max(Kmax, Lmax), normalized = TRUE) # Normalized on [0, 1]
    bfuncs <- lapply(tmp, function(y) {eval(parse(text = paste("function(x)", y)))}) # bfuncs[[k]] is b_{k - 1}
    bfuncs[[1]] <- function(x) rep(sqrt(2) / 2, length(x))
    cte1 <- 2 * sqrt(n - 1) / n
    betahat <- matrix(NA, nrow = Kmax + 1, ncol = Lmax + 1)
    for (k in 0:Kmax) { # Equ. (5.2)
      for (l in 0:Lmax) {
        betahat[k + 1, l + 1] <- cte1 * sum(R1 * sqrt(2) * bfuncs[[k + 1]](2 * UU[, 1] - 1) * sqrt(2) * bfuncs[[l + 1]](2 * UU[, 2] - 1))
      }
    }
    Aterm <- matrix(NA, nrow = Kmax + 1, ncol = Lmax + 1) # Denominator (5.3) without sqrt
    for (K in 0:Kmax) {
      for (L in 0:Lmax) {
        Aterm[K + 1 , L + 1] <- sum(betahat[1:(K + 1), 1:(L + 1)] ^ 2)
      }
    }
    cte2 <- 2 * sqrt(n - 2) / (n - 1)
    betahatminusi <- as.list(1:n) # Second equation after (C.1)
    for (i in 1:n) {
      rminusivec <- R1
      rminusivec[II[, 1] == i] <- R2[II[, 1] == i]
      rminusivec[i] <- 0
      betahatminusi[[i]] <- matrix(NA, nrow = Kmax + 1, ncol = Lmax + 1)
      for (k in 0:Kmax) {
        for (l in 0:Lmax) {
          betahatminusi[[i]][k + 1, l + 1]  <- cte2 * sum(rminusivec * sqrt(2) * bfuncs[[k + 1]](2 * UU[, 1] - 1) * sqrt(2) * bfuncs[[l + 1]](2 * UU[, 2] - 1)) 
        }
      }
    }
    
    Khat <- Lhat <- 0
    max <- -Inf
    for (K in 0:Kmax) {
      for (L in 0:Lmax) {
        term <- 0
        for (i in 1:n) {
          num  <- denom <- 0
          for (k in 0:K) {
            for (l in 0:L) {
              num <- num + betahatminusi[[i]][k + 1, l + 1] * sqrt(2) * bfuncs[[k + 1]](2 * UU[i, 1] - 1) * sqrt(2) * bfuncs[[l + 1]](2 * UU[i, 2] - 1)
              denom <- denom + (betahatminusi[[i]][k + 1, l + 1]) ^ 2
            }
          }
          term <- term + R1[i] * num / sqrt(denom)
        }
        term <- cte1 * term
        if (max < term) {
          max <- term
          Khat <- K
          Lhat <- L
        }
      }
    }
    if (Khat < 1) Khat  <-  1
    if (Lhat < 1) Lhat  <-  1
    
    etafunc <- function(H2) {
      B  <- 1 - H2
      return(2.0 * sqrt(-2.0 + B ^ 4 + sqrt(4.0 - 3.0 * B ^ 4)) / B ^ 2)
    }
    
    etahat  <-  etafunc(1.0 - betahat[1][1] / sqrt(Aterm[Khat + 1, Lhat + 1]))

    res <- vector("list")
    res$statistic <- etahat
    res$Khat <- Khat
    res$Lhat <- Lhat
      
    if (!is.null(conf.level)) {    
      XX <- cbind(x, y)
      n <- nrow(XX)
      RR <- apply(XX, 2, rank)
      etahatstar <- Zstar <- rep(NA, B1)
      etahatstarstar <- rep(NA, B2)
      for (b1 in 1:B1) {
        II <- sample(1:n, n, replace = TRUE)
        UUstar <- cbind(rbeta(n, RR[, 1][II], n + 1 - RR[, 1][II]), rbeta(n, RR[, 2][II], n + 1 - RR[, 2][II]))
        etacur <- HellCor(UUstar[, 1], UUstar[, 2], Kmax = Kmax, Lmax = Lmax, K = 0L, L = 0L, alpha = alpha, pval.comp = FALSE, conf.level = NULL)
        etahatstar[b1] <- etacur$Hcor
        RRstar <- apply(UUstar, 2, rank)
        for (b2 in 1:B2) {
          II <- sample(1:n, n, replace = TRUE)
          UUstarstar <- cbind(rbeta(n, RRstar[, 1][II], n + 1 - RRstar[, 1][II]), rbeta(n, RRstar[, 2][II], n + 1 - RRstar[, 2][II]))
          etahatstarstar[b2] <- HellCor(UUstarstar[, 1], UUstarstar[, 2], Kmax = Kmax, Lmax = Lmax,
                                        K = 0L, L = 0L, alpha = alpha, pval.comp = FALSE, conf.level = NULL)$Hcor
        }
        Vstarstar <- sd(etahatstarstar)
        Zstar[b1] <- (etacur$Hcor - res$statistic) / Vstarstar
      }
      Vstar <- sd(etahatstar)
      CIeta <- pmin(pmax(0, res$statistic - Vstar * quantile(Zstar, c(1 - (1 - conf.level) / 2, (1 - conf.level) / 2))), 1)
    }
  }

  if (pval.comp) {
    n <- length(x)
    twon <- 2 * n
    y <- runif(n)
    pval <- mean(replicate(10000, .C("hellcorC", as.double(cbind(runif(n), y)), as.integer(twon), statistic = as.double(0.0),
                                     as.integer(0), as.double(0.0), as.integer(4), as.integer(4), as.integer(0), as.integer(L),
                                     as.double(alpha), as.double(0.0), as.integer(B1), as.integer(B2),
                                     CIleft = as.double(0.0), CIright = as.double(0.0),
                                     PACKAGE = "HellCor")$statistic) > res$statistic)
  }

  
  if ((Kmax == 0L) & (Lmax == 0L)) {res$Khat <- NA; res$Lhat <- NA}
  return(list(Hcor = res$statistic, p.value = pval, conf.int = CIeta, Khat = res$Khat, Lhat = res$Lhat))
  
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
.datagenHeavisine <- function(n) {
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

.datagenBex <- function(n, depth = 2) { # depth >=1
  width <- 1 / 2 ^ depth
  which.be <- sample(1:2 ^ (2 * (depth - 1)), n, replace = TRUE)
  bex.center <- cbind(rep((1:2 ^ (depth - 1)) / 2 ^ (depth - 1),
                          2 ^ (depth - 1)),
                      rep((1:2 ^ (depth - 1)) / 2 ^ (depth - 1),
                          rep(2 ^ (depth - 1), 2 ^ (depth - 1)))) - 1 / 2 ^ depth
  x_inc <- runif(n) * 2 * width - width
  y_inc <- x_inc * sample(c(-1, 1), n, replace = TRUE)
  xy <- bex.center[which.be, ] + cbind(x_inc, y_inc)
  colnames(xy) <- c("x", "y")
  return(xy)
}

.datagenHilbertPeano <- function(n, depth = 1) { # depth >=1
  res <- .C("hilbertpeano", x = as.double(rep(0.0, 2 * n)), xlen = as.integer(2 * n), depth = as.integer(depth), setseed = 1L, PACKAGE = "HellCor")
  res <- matrix(res$x, nrow = n, ncol = 2)
  return(res)
}

.drawHilbertPeano <- function(depth = 1, drawlines = FALSE) { # depth >=1

  Double <- function(A) {
    #The matrix for H(n) is equal to Double(H(n-1)).
    m <- dim(A)[1]
    depth <- dim(A)[2]
    N <- m * depth
    B <- A + N
    C <- B + N
    D <- C + N
    E <- cbind(rbind(B, skew.transpose(A)), rbind(C, t(D)))
    return(E)
  }

  Rotate <- function(A) {
    #Rotates the matrix A clockwise.
    m <- dim(A)[1]
    depth <- dim(A)[2]
    N <- m * depth
    B <- matrix(0, m, depth)
    for (i in 1:m) for (j in 1:depth) B[j, depth + 1 - i] <- A[i, j]
    return(B)
  }

  skew.transpose <- function(A) {
    return(Rotate(Rotate(t(A))))
  }

  rowofx <- function(A, x) {
    #Returns the row index of the matrix A for entry equal to x.
    m <- dim(A)[1]
    depth <- dim(A)[2]
    for (i in 1:m) for (j in 1:depth) if (A[i, j] == x) return(i)
  }

  colofx <- function(A, x) {
    #Returns the column index of the matrix A for entry equal to x.
    m <- dim(A)[1]
    depth <- dim(A)[2]
    for (i in 1:m) for (j in 1:depth) if (A[i, j] == x) return(j)
  }

  Draw <- function(A, n) {
    #Draws a graphical representation of the matrix A.
    A <- Rotate(A)
    m <- dim(A)[1]
    depth <- dim(A)[2] - 1
    N <- m * (depth + 1)
    plot((colofx(A, 1) - 1) / depth, (rowofx(A, 1) - 1) / depth, pch = 19, cex = 0.5, ylim = c(0,1), xlim =c(0,1), ylab = character(1), xlab = character(1), axes = FALSE)
    for (i in 1:(N - 1)) lines(c((colofx(A, i) - 1) / depth, ((colofx(A, i + 1) - 1) / depth)), c((rowofx(A, i) - 1) / depth, ((rowofx(A, i + 1) - 1) / depth)), lwd = 1)
    points((colofx(A, N) - 1) / depth, (rowofx(A, N) - 1) / depth, pch = 19, cex = 0.5)
  }


  H <- function(depth) {
    if (depth == 0) return(matrix(c(2, 1, 3, 4), 2, 2))
    return(Double(H(depth - 1)))
  }
  
  H <- H(depth)
  A <- Rotate(H)
  m <- dim(A)[1]
  n <- dim(A)[2] - 1
  
  N <- m * (depth+1)
  pts <- runif(n, min = 1, max = N)
  left <- floor(pts)
  right <- ceiling(pts)
  alpha <- pts - left
  
  mat <- matrix(NA, ncol = 2, nrow = n)
  for (ind in 1:n) {
    i <- left[ind]
    x1 <- (rowofx(A, i) - 1) / depth
    x2 <- (rowofx(A, i + 1) - 1) / depth
    y1 <- (colofx(A, i) - 1) / depth
    y2 <- (colofx(A, i + 1) - 1) / depth
    mat[ind, ] <- c(x1 + alpha[ind] * (x2 - x1), y1 + alpha[ind] * (y2 - y1))
  }
  
  if (drawlines) Draw(H, n)
  
  return(mat)
}

# Generate samples from the 15 scenarios
.gendep <- function(n, alpha, beta, psi1, psi2, dimU2) {
  u1 <- runif(n)
  uu2 <- matrix(runif(n * dimU2), ncol = dimU2)
  w1 <- runif(n)
  ww2 <- matrix(runif(n * dimU2), ncol = dimU2)
  v1 <- runif(n)
  v2 <- runif(n)
  X1X2 <- cbind(psi1(u1, uu2), alpha * psi2(u1, uu2) + (1 - alpha) * psi2(w1, ww2)) + beta * cbind(qnorm(v1), qnorm(v2))
  return(X1X2)
}

.datagen15 <- function(n, alpha, beta, type) {
  switch(type,
         "W" = { # corrected
           dimU2 <- 0
           psi1 <- function(u1, uu2) (2 * u1 - 1)
           psi2 <- function(u1, uu2) 4 * ((2 * u1 - 1) ^ 2 - 0.5) ^ 2
         },
         "Diamond" = {
           dimU2 <- 1
           psi1 <- function(u1, uu2) 2 * u1
           psi2 <- function(u1, uu2) -1 + (2 * (uu2 > 0.5) - 1) * (1 - abs(2 * u1 - 1))
         },
         "Parabola" = { # corrected
           dimU2 <- 0
           psi1 <- function(u1, uu2) (2 * u1 - 1)
           psi2 <- function(u1, uu2) (2 * u1 - 1) ^ 2 / 2
         },
         "2Parabolas" = { # corrected
           dimU2 <- 1
           psi1 <- function(u1, uu2) (2 * u1 - 1)
           psi2 <- function(u1, uu2) (2 * u1 - 1) ^ 2 * (2 * (uu2 > 0.5) - 1)
         },
         "Circle" = { # corrected
           dimU2 <- 0
           psi1 <- function(u1, uu2) sin(2 * pi * u1 - pi)
           psi2 <- function(u1, uu2) cos(2 * pi * u1 - pi)
         },
         "4Clouds" = {
           dimU2 <- 0
           psi1 <- function(u1, uu2) (2 * (u1 > 0.5) - 1) + u1 / 100
           psi2 <- function(u1, uu2) (2 * (u1 > 0.5) - 1) + u1 / 100
         },
         "Cubic" = {
           dimU2 <- 0
           psi1 <- function(u1, uu2) (2 * u1 - 1) ^ 3
           psi2 <- function(u1, uu2) 4 * (2 * u1 - 1) ^ 3 + (2 * u1 - 1) ^ 2 - 4 * (2 * u1 - 1)
         },
         "Sine" = {
           dimU2 <- 0
           psi1 <- function(u1, uu2) 4 * u1
           psi2 <- function(u1, uu2) 4 * sin(4 * pi * u1)
           },
         "Wedge" = {
           dimU2 <- 1
           psi1 <- function(u1, uu2) u1
           psi2 <- function(u1, uu2) u1 * (2 * (uu2 > 0.5) - 1)
         },
         "Cross" = {
           dimU2 <- 1
           psi1 <- function(u1, uu2) u1
           psi2 <- function(u1, uu2) (2 * u1 - 1) * (2 * (uu2 > 0.5) - 1)
         },
         "Spiral" = {
           dimU2 <- 0
           psi1 <- function(u1, uu2) u1 *  6 * pi * cos(u1 * 6 * pi) / 3
           psi2 <- function(u1, uu2) u1 *  6 * pi * sin(u1 * 6 * pi) / 3
         },
         "4Circles" = {
           dimU2 <- 2
           psi1 <- function(u1, uu2) sqrt(2) * cos(u1 * 2 * pi) + (2 * (uu2[, 1] > 0.5) - 1)
           psi2 <- function(u1, uu2) sqrt(2) * sin(u1 * 2 * pi) + (2 * (uu2[, 2] > 0.5) - 1)
         },
         "Heavisine" = {
           dimU2 <- 0
           psi1 <- function(u1, uu2) 4 * u1
           psi2 <- function(u1, uu2) {y  <- (4 * sin(4 * pi * u1) - sign(u1 - 0.3) - sign(0.72 - u1)); return(y / sqrt(var(y)))}
         },
         "Doppler" = {
           dimU2 <- 0
           psi1 <- function(u1, uu2) 10 * u1
           psi2 <- function(u1, uu2) {eps <- 0.05; y <- sqrt(u1 * (1 - u1)) * sin((2 * pi * (1 + eps))/ (u1 + eps)); return(y / sqrt(var(y)))}
         },
         "5to15Clouds" = {
           dimU2 <- 0
           psi1 <- function(u1, uu2) {mu <- 1; ind1 <- (u1 < (1 / 5)); ind2  <- (u1 >= (1 / 5)) & (u1 < (2 / 5)); ind3  <- (u1 >= (2 / 5)) & (u1 < (3 / 5)); ind4  <- (u1 >= (3 / 5)) & (u1 < (4 / 5)); ind5 <- (u1 >= (4 / 5)); u1[ind1] <- u1[ind1] / 100 - mu; u1[ind2] <- u1[ind2] / 100 - mu; ; u1[ind3] <- u1[ind3] / 100; u1[ind4] <- u1[ind4] / 100 + mu; u1[ind5] <- u1[ind5] / 100 + mu; return(u1)}
           psi2 <- function(u1, uu2) {mu <- 1; ind1 <- (u1 < (1 / 5)); ind2  <- (u1 >= (1 / 5)) & (u1 < (2 / 5)); ind3  <- (u1 >= (2 / 5)) & (u1 < (3 / 5)); ind4  <- (u1 >= (3 / 5)) & (u1 < (4 / 5)); ind5 <- (u1 >= (4 / 5)); u1[ind1] <- u1[ind1] / 100 - mu; u1[ind2] <- u1[ind2] / 100 + mu; ; u1[ind3] <- u1[ind3] / 100; u1[ind4] <- u1[ind4] / 100 + mu; u1[ind5] <- u1[ind5] / 100 - mu; return(u1)}
         }
         )
  res <- .gendep(n, alpha, beta, psi1, psi2, dimU2)
  return(res)
}
