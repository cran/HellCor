#include <R.h>
#include "Rmath.h"

extern "C" {

  void hellcorC(double *x, int *xlen, double *statistic,
		int *pvalcomp, double *pvalue, int *Kmaxvalue, int *Lmaxvalue, int *Kset, int *Lset, double *alphavalue,
		double *conflevel, int *B1, int *B2, double *CIetaleft, double *CIetaright) {

    int n = xlen[0];
    
    int N = n / 2;
  
    if (N>3) {
      double LegendrePoly(double x, int k);
      double etafunc(double H2);
      double tmp = 0.0, tmp1, tmp2, max, cte1, cte2, etahat, num, denom, term;
      double *u1, *u2, *R1, **dist, *u1ties, *u2ties, *t1, *t2, *weight, *R2, **Rminusi, **Legtmp1, **Legtmp2, **betahat, **Aterm;
      double ***betahatminusi;
      int **II;
      int Kmax = Kmaxvalue[0], Lmax = Lmaxvalue[0], i, iprime, j, k, l, K, L, Khat = 0, Lhat = 0;
      double alpha = alphavalue[0];
      cte1 = 2.0 * sqrt((double)(N - 1)) / (double)N;
      cte2 = 2.0 * sqrt((double)(N - 2)) / (double)(N - 1);
      u1 = new double[N];
      u2 = new double[N];
      R1 = new double[N];
      dist = new double*[N];
      for (i = 0; i < N; i++) dist[i] = new double[N];
      u1ties = new double[N];
      u2ties = new double[N];
      t1 = new double[N];  
      t2 = new double[N];
      weight = new double[N];
      Legtmp1 = new double*[N];
      for (i = 0; i < N; i++) Legtmp1[i] = new double[Kmax + 1];
      Legtmp2 = new double*[N];
      for (i = 0; i < N; i++) Legtmp2[i] = new double[Lmax + 1];
      R2 = new double[N];
      Rminusi = new double*[N];
      for (i = 0; i < N; i++) Rminusi[i] = new double[N];
      II = new int*[N];
      for (i = 0; i < N; i++) II[i] = new int[2];
      betahat = new double*[Kmax + 1];
      for (k = 0; k <= Kmax; k++) betahat[k] = new double[Lmax + 1];
      Aterm = new double*[Kmax + 1];
      for (k = 0; k <= Kmax; k++) Aterm[k] = new double[Lmax + 1];
      betahatminusi = new double**[N];
      for (i = 0; i < N; i++) betahatminusi[i] = new double*[Kmax + 1];
      for (i = 0; i < N; i++) {for (k = 0; k <= Kmax; k++) betahatminusi[i][k] = new double[Lmax + 1];}

      // Variables for the computation of the confidence interval
      int *III;
      III = new int[N];
      double *etahatstar, *Zstar, *etahatstarstar;
      etahatstar = new double[B1[0]];
      Zstar = new double[B1[0]];
      etahatstarstar = new double[B2[0]];
      double R_unif_index(double);
      int b1, b2;
      double *UUstar;
      UUstar = new double[n];
      double rbeta(double aa, double bb);
      double *ustar1, *ustar2, *ustar1ties, *ustar2ties;
      ustar1 = new double[N];
      ustar2 = new double[N];
      ustar1ties = new double[N];
      ustar2ties = new double[N];
      double *RRstar;
      RRstar = new double[n];
      double *Hcor;
      Hcor = new double[1];
      int *pvalcomp2;
      pvalcomp2 = new int[1];
      pvalcomp2[0] = 0;
      int *Kmaxvalue2, *Lmaxvalue2, *Kset2, *Lset2, *B12, *B22;
      Kmaxvalue2 = new int[1];
      Lmaxvalue2 = new int[1];
      Kset2 = new int[1];
      Lset2 = new int[1];
      B12 = new int[1];
      B22 = new int[1];
      double *UUstarstar;
      UUstarstar = new double[n];
      double Vstarstar, sdetahatstarstar;
      double Vstar, sdetahatstar;
      void R_rsort (double* x, int n);

	double *index;
	index = new double[2];
	int *lo, *hi;
	lo = new int[2];
	hi = new int[2];
	double *qs;
	qs = new double[2];
	int *ii;
	ii = new int[2];
	double h;

      
      //  Now we define the "copula" observations (U_{i1}, U_{i,2}), i = 1, ..., n;
      //  see first paragraph of Section 5.2     
      for (i = 0; i < N; i++) {
	u1[i] = 0.0;
	u2[i] = 0.0;
	R1[i] = R_PosInf;
	for (j = 0; j < N; j++) {
	    dist[i][j] = 0.0;
	    if (x[j] <= x[i]) u1[i] = u1[i] + 1.0;
	    if (x[j + N] <= x[i + N]) u2[i] = u2[i] + 1.0;
	}
      }
      for (i = 0; i < N; i++) {
	k = 0;
	l = 0;
	for (j = 0; j < N; j++) {
	  if (u1[i] == u1[j]) k = k + 1;
	  if (u2[i] == u2[j]) l = l + 1;
	}
	if (k > 1) u1ties[i] = u1[i] - (double)(k - 1) * k / (double)(2 * k); else u1ties[i] = u1[i];
	if (l > 1) u2ties[i] = u2[i] - (double)(l - 1) * l / (double)(2 * l); else u2ties[i] = u2[i];
      }
      // (u1, u2) contains the empirical copula, i.e., UU <- apply(XX, 2, rank) / (n+1)
      // with ties.method = "average":
      for (i = 0; i < N; i++) {
	u1[i] = u1ties[i] / double(N + 1);
	u2[i] = u2ties[i] / double(N + 1);
      }

      // Now we apply the marginal transformations in Section B of the Appendix
      // and we compute the weights appearing in Equ (B.1)
      for (i = 0; i < N; i++) {
	t1[i] = Rf_qbeta(u1[i], alpha, alpha, 1, 0); // T_{i1}
	t2[i] = Rf_qbeta(u2[i], alpha, alpha, 1, 0); // T_{i2}
	// \sqrt{\xi_1}(T_{i1})\sqrt{\xi_2}(T_{i2})
	weight[i] = sqrt(Rf_dbeta(t1[i], alpha, alpha, 0) * Rf_dbeta(t2[i], alpha, alpha, 0));
      }

      // dist[i][j] = || T_j - T_i ||
      for (i = 0; i < (N - 1); i++) {
	for (j = i + 1; j < N; j++) {
	  dist[i][j] = sqrt(R_pow(t1[i] - t1[j], 2.0) + R_pow(t2[i] - t2[j], 2.0));
	}
      }

      // Compute, for i=1,...,N, the smallest R[i] and second smallest R2[i] of dist[i][j], for j != i
      for (i = 0; i < N; i++) {
	R1[i] = R_PosInf; // smallest distance from a point T_j (j != i) to T_i
	R2[i] = R_PosInf; //second smallest distance from a point T_j (j != i) to T_i
	II[i][0] = 1; // index of which point is the closest to the ith
	II[i][1] = 2; // index of which point is the second closest to the ith
	for (j = 0; j < N; j ++) {
	  if (j > i) tmp = dist[i][j]; else if (j < i) tmp = dist[j][i];
	  if (j != i) {
	    if (tmp < R1[i]) {
	      R2[i] = R1[i];
	      II[i][1] = II[i][0];
	      R1[i] = tmp;
	      II[i][0] = j + 1;
	    }  else if (tmp < R2[i]) {
	      R2[i] = tmp;
	      II[i][1] = j + 1;
	    }
	  }
	}
      }


      for (i = 0; i < N; i++) {
	R1[i] = R1[i] * weight[i];
	R2[i] = R2[i] * weight[i];
      }
	    
      // Apply the weights desribed in Equ (B.1)
      for (i = 0; i < N; i++) {
	for (iprime = 0; iprime < N; iprime++) {
	  Rminusi[i][iprime] = R1[iprime];
	}
	for (iprime = 0; iprime < N; iprime++) {
	  if (II[iprime][0] == i + 1) {
	    Rminusi[i][iprime] = R2[iprime];
	  }
	if (iprime == i) Rminusi[i][iprime] = 0.0;
	}
      }
    

	


      
      // We compute the b_k(\hat{U}_{ij}), i=1,...,n, j=1,2 in Equ (5.2)
      for (i = 0; i < N; i++) {
	for (k = 0; k <= Kmax; k++) {
	  Legtmp1[i][k] = sqrt(2.0) * LegendrePoly(2.0 * u1[i] - 1.0, k);
	  Legtmp2[i][k] = sqrt(2.0) * LegendrePoly(2.0 * u2[i] - 1.0, k);
	}
      }

      // We compute \hat{\beta}_{kl} in Equ (5.2)
      for (k = 0; k <= Kmax; k++) {
 	for (l = 0; l <= Lmax; l++) {
 	  tmp = 0.0;
 	  for (i = 0; i < N; i++) {
 	    tmp = tmp +  R1[i] * Legtmp1[i][k] * Legtmp2[i][l];
 	  }
	  betahat[k][l] =  cte1 * tmp;		
 	}
      }

      if ((Kmax > 0) || (Lmax > 0)) {
	
	// Aterm size is (Kmax + 1) x (Lmax + 1)
	// For K = 0, ..., Kmax+1 and L = 0, ..., Lmax+1
	//    Aterm[K][L] = \sum_{k=0}^K \sum_{l=0}^L \hat{\beta}_{kl}^2
	// Aterm[K][L] is thus the term in the denominator of Equ (5.3), without the 1/2 power
	for (L = 0; L <= Lmax; L++) {
	  tmp = 0.0;
	  for (K = 0; K <= Kmax; K++) {
	    tmp = tmp + R_pow(betahat[K][L], 2.0);
	    Aterm[K][L] = tmp; // For L fixed, Aterm[K][L] = \sum_{k=0}^K \hat{\beta}_{kL}^2
	  }
	}
	for (K = 0; K <= Kmax; K++) {
	  tmp = 0.0;
	  for (L = 0; L <= Lmax; L++) {
	    tmp = tmp + Aterm[K][L];
	    Aterm[K][L] = tmp; // For K fixed, Aterm[K][L] = \sum_{l=0}^L \sum_{k=0}^K \hat{\beta}_{kl}^2
	  }
	}

	
	// We compute \hat{\beta}_{kl}^{(-i)} in Appendix C.
	for (i = 0; i < N; i++) {
	  for (k = 0; k <= Kmax; k++) {
	    for (l = 0; l <= Lmax; l++) {
	      tmp = 0.0;
	      for (iprime = 0; iprime < N; iprime++) {
		if (iprime != i) {
		  //	Rprintf(" rminusvec: %g ", Rminusi[iprime][i]);
		  tmp = tmp +  Rminusi[i][iprime] * Legtmp1[iprime][k] * Legtmp2[iprime][l];
		}
	      }
	      betahatminusi[i][k][l] =  cte2 * tmp;
	    }
	  }
	}

	//	Rprintf(" %g ", betahatminusi[0][0][0]);


	
	// We compute Khat and Lhat. No need to premultiply by a constant here since it won't change the minimization results
	Khat = 0;
	Lhat = 0;
	max = R_NegInf;
	for (K = 0; K <= Kmax; K++) {	
	  for (L = 0; L <= Lmax; L++) {
	    term = 0.0;
	    for (i = 0; i < N; i++) {
	      num = 0.0;
	      denom = 0.0;
	      for (k = 0; k <= K; k++) {	
		for (l = 0; l <= L; l++) {	
		  num = num + betahatminusi[i][k][l] * Legtmp1[i][k] * Legtmp2[i][l];
		  denom = denom + R_pow(betahatminusi[i][k][l], 2.0);
		}
	      }
	      term = term + R1[i] * num  / sqrt(denom);	    
	    }
	    term = cte1 * term;
	    if (max < term) {
	      max = term;
	      Khat = K;
	      Lhat = L;
	    }
	  }
	}
	if (Khat < 1) Khat = 1;
	if (Lhat < 1) Lhat = 1;
	Kmaxvalue[0] = Khat;
	Lmaxvalue[0] = Lhat;
	
      
	// H ^ 2 = 1 - B (see Equ (4.2))
	// Just above Equ (5.2), we see that B = \beta_{00}
	// The denominator is from Equ (5.3)
	// We use the transformation etafunc in Equ (4.3) to get Equ (5.4)
	etahat = etafunc(1.0 - betahat[0][0] / sqrt(Aterm[Khat][Lhat]));

      } else if ((Kset[0] > 0) || (Lset[0] > 0)) {

	// Aterm[K][L] is thus the term in the denominator of Equ (5.3), without the 1/2 power
	tmp2 = 0.0;
	for (k = 0; k <= Kset[0]; k++) {
	  for (l = 0; l <= Lset[0]; l++) {
	    tmp1 = 0.0;
	    for (i = 0; i < N; i++) {
	      tmp1 = tmp1 +  R1[i] * sqrt(2.0) * LegendrePoly(2.0 * u1[i] - 1.0, k) * sqrt(2.0) * LegendrePoly(2.0 * u2[i] - 1.0, l);
	    }
	    tmp2 = tmp2 + R_pow(cte1 * tmp1, 2.0);
	  }
	}
	
	etahat = etafunc(1 - betahat[0][0] / sqrt(tmp2));
	
      } else {
	
	etahat = etafunc(1 - betahat[0][0]);
	
      }
      
      statistic[0] = etahat; // Here is the test statistic value


      if ((B1[0] > 0) && (B2[0] > 0) && (conflevel[0] > 0.0)) {

	GetRNGstate();

	Vstar = 0.0;
	for (b1 = 1; b1 <= B1[0]; b1++) {
      
	  for (i = 0; i < N; i++) III[i] =  (int)R_unif_index(N) + 1; // sample(1:N, N, replace = TRUE)

	  
	  //      UUstar <- cbind(rbeta(n, RR[, 1][II], n + 1 - RR[, 1][II]), rbeta(n, RR[, 2][II], n + 1 - RR[, 2][II]))
	  for (i = 0; i < N; i++) UUstar[i] = rbeta((double)(N + 1) * u1[III[i] - 1], (double)(N + 1) * (1.0 - u1[III[i] - 1]));
	  for (i = 0; i < N; i++) UUstar[N + i] = rbeta((double)(N + 1) * u2[III[i] - 1], (double)(N + 1) * (1.0 - u2[III[i] - 1]));
	  
	  Hcor[0] = 0.0;
	  Kmaxvalue2[0] = Kmax;
	  Lmaxvalue2[0] = Lmax;
	  Kset2[0] = 0;
	  Lset2[0] = 0;
	  B12[0] = 0;
	  B22[0] = 0;
	  hellcorC(UUstar, xlen, Hcor, pvalcomp2, pvalue, Kmaxvalue2, Lmaxvalue2, Kset2, Lset2,
		   alphavalue, conflevel, B12, B22, CIetaleft, CIetaright);
	  
	  etahatstar[b1 - 1] = Hcor[0];

	  for (i = 0; i < N; i++) {
	    ustar1[i] = 0.0;
	    ustar2[i] = 0.0;
	    for (j = 0; j < N; j++) {
	      if (UUstar[j] <= UUstar[i]) ustar1[i] = ustar1[i] + 1.0;
	      if (UUstar[j + N] <= UUstar[i + N]) ustar2[i] = ustar2[i] + 1.0;
	    }
	  }
	  for (i = 0; i < N; i++) {
	    k = 0;
	    l = 0;
	    for (j = 0; j < N; j++) {
	      if (ustar1[i] == ustar1[j]) k = k + 1;
	      if (ustar2[i] == ustar2[j]) l = l + 1;
	    }
	    if (k > 1) ustar1ties[i] = ustar1[i] - (double)(k - 1) * k / (double)(2 * k); else ustar1ties[i] = ustar1[i];
	    if (l > 1) ustar2ties[i] = ustar2[i] - (double)(l - 1) * l / (double)(2 * l); else ustar2ties[i] = ustar2[i];
	  }
	  // RRstar <- apply(UUstar, 2, rank)
	  // with ties.method = "average":
	  for (i = 0; i < N; i++) {
	    RRstar[i] = ustar1ties[i];
	    RRstar[N + i] = ustar2ties[i];
	  }
	  
	  Vstarstar = 0.0;
	  for (b2 = 1; b2 <= B2[0]; b2++) {
	    //        II <- sample(1:n, n, replace = TRUE)
	    for (i = 0; i < N; i++) III[i] =  (int)R_unif_index(N) + 1; // sample(1:N, N, replace = TRUE)

	    //        UUstarstar <- cbind(rbeta(n, RRstar[, 1][II], n + 1 - RRstar[, 1][II]), rbeta(n, RRstar[, 2][II], n + 1 - RRstar[, 2][II]))
	    for (i = 0; i < N; i++) UUstarstar[i] = rbeta(RRstar[III[i] - 1], (double)(N + 1) - RRstar[III[i] - 1]);
	    for (i = 0; i < N; i++) UUstarstar[N + i] = rbeta(RRstar[N + III[i] - 1], (double)(N + 1) - RRstar[N + III[i] - 1]);

	    Hcor[0] = 0.0;
	    Kmaxvalue2[0] = Kmax;
	    Lmaxvalue2[0] = Lmax;
	    Kset2[0] = 0;
	    Lset2[0] = 0;
	    B12[0] = 0;
	    B22[0] = 0;
	    hellcorC(UUstarstar, xlen, Hcor, pvalcomp2, pvalue, Kmaxvalue2, Lmaxvalue2, Kset2, Lset2, alphavalue, conflevel, B12, B22, CIetaleft, CIetaright);
	    
	    etahatstarstar[b2 - 1] = Hcor[0];
	    Vstarstar = Vstarstar + etahatstarstar[b2 - 1];
	  }
	  
	  //      Vstarstar <- sd(etahatstarstar)
	  Vstarstar = Vstarstar / B2[0];
	  sdetahatstarstar = 0.0;
	  for (b2 = 1; b2 <= B2[0]; b2++) {
	    sdetahatstarstar = sdetahatstarstar + R_pow(etahatstarstar[b2 - 1] - Vstarstar, 2.0);
	  }
	  Vstarstar = sqrt(sdetahatstarstar / (double)(B2[0] - 1));

	  // Zstar[b1] <- (etacur$Hcor - res$statistic) / Vstarstar
	  Zstar[b1 - 1] = (etahatstar[b1 - 1] - statistic[0]) / Vstarstar;
	  
	  Vstar = Vstar + etahatstar[b1 - 1];
	}

   
	//    Vstar <- sd(etahatstar)
	Vstar = Vstar / B1[0];
	sdetahatstar = 0.0;
	for (b1 = 1; b1 <= B1[0]; b1++) {
	  sdetahatstar = sdetahatstar + R_pow(etahatstar[b1 - 1] - Vstar, 2.0);
	}
	Vstar = sqrt(sdetahatstar / (double)(B1[0] - 1));

	//     CIeta <- pmin(pmax(0, res$statistic - Vstar * quantile(Zstar, c(1 - (1 - conf.level) / 2, (1 - conf.level) / 2))), 1)
	R_rsort (Zstar, B1[0]);

	// R function quantile(x, probs)    
	// 	n <- length(x) 
	//	np <- length(probs)    here np = 2
	//      index <- 1 + (n - 1) * probs
	index[0] = 1.0 + (double)(B1[0] - 1) * (1.0 - (1.0 - conflevel[0]) / 2.0);
	index[1] = 1.0 + (double)(B1[0] - 1) * (1.0 - conflevel[0]) / 2.0;
	//      lo <- floor(index)
	//      hi <- ceiling(index)
	lo[0] = (int)std::floor(index[0]);
	lo[1] = (int)std::floor(index[1]);
	hi[0] = (int)std::ceil(index[0]);
	hi[1] = (int)std::ceil(index[1]);


	//    qs <- x[lo]
	qs[0] = Zstar[lo[0] - 1];
	qs[1] = Zstar[lo[1] - 1];

	// i <- which(index > lo & x[hi] != qs)
	if ((index[0] > lo[0]) && (Zstar[hi[0] -1] != qs[0])) ii[0] = 1; else ii[0] = 0;
	if ((index[1] > lo[1]) && (Zstar[hi[1] -1] != qs[1])) ii[1] = 2; else ii[1] = 0;

	//	h <- (index - lo)[i]
	//      qs[i] <- (1 - h) * qs[i] + h * x[hi[i]]
	for (i = 0; i <= 1; i++) {
	  if (ii[i] != 0) {
	    h = index[ii[i] - 1] - (double)lo[ii[i] - 1];
	    qs[ii[i] - 1] = (1.0 - h) * qs[ii[i] - 1] + h * Zstar[hi[ii[i] - 1] - 1];
	  }
	}


	//	CIetaleft[0] = statistic[0] - Vstar * Zstar[(int)(B1[0] * (1.0 - (1.0 - conflevel[0]) / 2.0)) - 1];
	//	CIetaright[0] = statistic[0] - Vstar * Zstar[(int)(B1[0] * (1.0 - conflevel[0]) / 2.0) - 1];

	
	CIetaleft[0] = statistic[0] - Vstar * qs[0];
	CIetaright[0] = statistic[0] - Vstar * qs[1];


	if (CIetaleft[0] < 0) CIetaleft[0] = 0.0;
	if (CIetaright[0] < 0) CIetaright[0] = 0.0;
	if (CIetaleft[0] > 1) CIetaleft[0] = 1.0;
	if (CIetaright[0] > 1) CIetaright[0] = 1.0;
	

	PutRNGstate();
	
      }


      if (pvalcomp[0] == 0) {
      
	pvalue[0] = 0.0;

      } else {

	// TODO
	
      }

    
      
 // We free the unused array of pointers
      delete[] u1;
      delete[] u2;
      delete[] R1;
      for (i = 0; i < N; i++) delete[] dist[i];
      delete[] dist;
      delete[] u1ties;
      delete[] u2ties;
      delete[] t1;
      delete[] t2;
      delete[] weight;
      delete[] R2;
      for (i = 0; i < N; i++) delete[] Rminusi[i];
      delete[] Rminusi;
      for (i = 0; i < N; i++) delete[] II[i];
      delete[] II;
      for (i = 0; i < N; i++) delete[] Legtmp1[i];
      delete[] Legtmp1;
      for (i = 0; i < N; i++) delete[] Legtmp2[i];
      delete[] Legtmp2;
      for (k = 0; k <= Kmax; k++) delete[] betahat[k];
      delete[] betahat;
      for (k = 0; k <= Kmax; k++) delete[] Aterm[k];
      delete[] Aterm;
      for (i = 0; i < N; i++) {for (k = 0; k <= Kmax; k++) delete[] betahatminusi[i][k];}
      for (i = 0; i < N; i++) delete[] betahatminusi[i];
      delete[] betahatminusi;     

      delete[] III;
      delete[] etahatstar;
      delete[] Zstar;
      delete[] etahatstarstar;
      delete[] UUstar;
      delete[] ustar1;
      delete[] ustar2;
      delete[] ustar1ties;
      delete[] ustar2ties;
      delete[] RRstar;
      delete[] Hcor;
      delete[] pvalcomp2;
      delete[] Kmaxvalue2;
      delete[] Lmaxvalue2;
      delete[] Kset2;
      delete[] Lset2;
      delete[] B12;
      delete[] B22;
      delete[] UUstarstar;

      delete[] index;
      delete[] lo;
      delete[] hi;
      delete[] qs;
      delete[] ii;

    }

// We return
    return;
          
  }

  
  double etafunc(double H2) { // This is h^{-1}(B) in Equ (4.3) but note that above (4.2) we have h^{-1}(H^2) with H^2 = 1 - B
  double B = 1.0 - H2;
  double B2 = B * B;
  double B4 = B2 * B2;
  return 2.0 * sqrt(-2.0 + B4 + sqrt(4.0 - 3.0 * B4)) / B2;
}

  double LegendrePoly(double x, int k) {
    /*
The polynomials M_n below is equal to \sqrt{(2n+1)/2} * P_n where P_n is defined 
at https://en.wikipedia.org/wiki/Legendre_polynomials#Shifted_Legendre_polynomials
The P_n are orthogonal on [-1,1] with P_0=1, P_1=x
We have 
\int_{-1}^1 P_M(x)P_n(x)dx = (2/(2n+1))\delta_{mn}
(n+1)P_{n+1}(x) = (2n+1)xP_n(x)-nP_{n-1}(x)
The orthonormal Legendre polynomials on [0,1] are defined as \sqrt{2} * P_n(2x-1)
     */
    if (k == 0) {
      return 0.707106781186547;
    } else if (k == 1) {
      return 1.22474487139159 * x;
    } else if (k == 2) {
      return -0.790569415042095 + 2.37170824512628 * R_pow(x, 2.0);
    } else if (k == 3) {
      return -2.80624304008046*x + 4.67707173346743 * R_pow(x, 3.0);
    } else if (k == 4) {
      return 0.795495128834866 - 7.95495128834866 * R_pow(x, 2.0) + 9.28077650307344 * R_pow(x, 4.0);
    } else if (k == 5) {
      return 4.39726477483446*x - 20.5205689492275 * R_pow(x, 3.0) + 18.4685120543048 * R_pow(x, 5.0);
    } else if (k == 6) {
      return -0.796721798998873 + 16.7311577789763 * R_pow(x, 2.0) - 50.193473336929 * R_pow(x, 4.0) +  
	36.8085471137479 * R_pow(x, 6.0);
    } else if (k == 7) {
      return -5.99071547271275*x + 53.9164392544148 * R_pow(x, 3.0) - 118.616166359713 * R_pow(x, 5.0) +  
	73.4290553655363 * R_pow(x, 7.0);
    } else if (k == 8) {
      return 0.797200454373381 - 28.6992163574417 * R_pow(x, 2.0) + 157.845689965929 * R_pow(x, 4.0) -  
	273.599195940944 * R_pow(x, 6.0) + 146.570997825506 * R_pow(x, 8.0);
    } else if (k == 9) {
      return 7.58511879271573 * x - 111.248408959831 * R_pow(x, 3.0) + 433.86879494334 * R_pow(x, 5.0) -  
	619.812564204771 * R_pow(x, 7.0) + 292.689266430031 * R_pow(x, 9.0);
    } else if (k == 10) {
      return -0.797434890624405 + 43.8589189843422 * R_pow(x, 2.0) - 380.110631197633 * R_pow(x, 4.0) +  
	1140.3318935929 * R_pow(x, 6.0) - 1384.68872793423 * R_pow(x, 8.0) + 584.646351794454 * R_pow(x, 10.0);
    } else if (k == 11) {
      return -9.17998960606603 * x + 198.899774798097 * R_pow(x, 3.0) - 1193.39864878858 * R_pow(x, 5.0) +  
	2898.25386134371 * R_pow(x, 7.0) - 3059.26796475169 * R_pow(x, 9.0) + 1168.0841319961 * R_pow(x, 11.0);
    } else if (k == 12) {
      return 0.797566730732873 - 62.2102049971641 * R_pow(x, 2.0) + 777.627562464552 * R_pow(x, 4.0) -  
	3525.2449498393 * R_pow(x, 6.0) + 7176.39150503 * R_pow(x, 8.0) - 6697.96540469467 * R_pow(x, 10.0) +  
	2334.13945921178 * R_pow(x, 12.0);
    } else if (k == 13) {
      return 10.7751235804364 * x - 323.253707413091 * R_pow(x, 3.0) + 2747.65651301127 * R_pow(x, 5.0) -  
	9943.89976137412 * R_pow(x, 7.0) + 17401.8245824047 * R_pow(x, 9.0) - 14554.2532871021 * R_pow(x, 11.0) +  
	4664.82477150709 * R_pow(x, 13.0);
    } else if (k == 14) {
      return -0.797648110941312 + 83.7530516488378 * R_pow(x, 2.0) - 1423.80187803024 * R_pow(x, 4.0) +  
	9017.41189419154 * R_pow(x, 6.0) - 27052.2356825746 * R_pow(x, 8.0) + 41480.0947132811 * R_pow(x, 10.0) -  
	31424.3141767281 * R_pow(x, 12.0) + 9323.69761287537 * R_pow(x, 14.0);
    } else if (k == 15) {
      return -12.37042008527 * x + 490.693330049044 * R_pow(x, 3.0) - 5593.9039625591 * R_pow(x, 5.0) +  
	27969.5198127955 * R_pow(x, 7.0) - 71477.6617438108 * R_pow(x, 9.0) + 97469.5387415601 * R_pow(x, 11.0) -  
	67478.9114364647 * R_pow(x, 13.0) + 18637.0326824522 * R_pow(x, 15.0);
    } else if (k == 16) {
      return 0.79770183004505 - 108.487448886127 * R_pow(x, 2.0) + 2404.80511697581 * R_pow(x, 4.0) -  
	20200.3629825968 * R_pow(x, 6.0) + 82965.7765356655 * R_pow(x, 8.0) - 184368.392301479 * R_pow(x, 10.0) +  
	226270.299642724 * R_pow(x, 12.0) - 144216.234937121 * R_pow(x, 14.0) + 37255.8606920895 * R_pow(x, 16.0);
    } else if (k == 17) {
      return 13.9658239139855 * x - 707.601744975264 * R_pow(x, 3.0) + 10401.7456511364 * R_pow(x, 5.0) -  
	68354.3285646105 * R_pow(x, 7.0) + 237341.41862712 * R_pow(x, 9.0) - 466052.240213253 * R_pow(x, 11.0) +  
	519827.498699398 * R_pow(x, 13.0) - 306945.761136787 * R_pow(x, 15.0) + 74479.4861581911 * R_pow(x, 17.0);
    } else if (k == 18) {
      return -0.797739132849908 + 136.413391717334 * R_pow(x, 2.0) - 3819.57496808536 * R_pow(x, 4.0) +  
	40996.7713241162 * R_pow(x, 6.0) - 219625.560664908 * R_pow(x, 8.0) + 658876.681994724 * R_pow(x, 10.0) -  
	1158025.68350588 * R_pow(x, 12.0) + 1183476.79742909 * R_pow(x, 14.0) - 650912.238585997 * R_pow(x, 16.0) +  
	148901.492486993 * R_pow(x, 18.0);
    } else if (k == 19) {
      return -15.5613022863318 * x + 980.362044038905 * R_pow(x, 3.0) - 18038.6616103158 * R_pow(x, 5.0) +  
	150322.180085965 * R_pow(x, 7.0) - 676449.810386844 * R_pow(x, 9.0) + 1783367.68192895 * R_pow(x, 11.0) -  
	2835097.34050244 * R_pow(x, 13.0) + 2673091.77818801 * R_pow(x, 15.0) - 1375856.06230265 * R_pow(x, 17.0) +  
	297699.849738001 * R_pow(x, 19.0);
    } else if (k == 20) {
      return 0.797766083041056 - 167.530877438622 * R_pow(x, 2.0) + 5779.81527163245 * R_pow(x, 4.0) -  
	77064.203621766 * R_pow(x, 6.0) + 520183.37444692 * R_pow(x, 8.0) - 2011375.71452809 * R_pow(x, 10.0) +  
	4723685.39017961 * R_pow(x, 12.0) - 6851939.2472935 * R_pow(x, 14.0) + 5995446.84138182 * R_pow(x, 16.0) -  
	2899758.60302127 * R_pow(x, 18.0) + 595213.607988577 * R_pow(x, 20.0);
    } else {
      return sqrt((double)(2 * k  + 1)) * ( sqrt((double)(2 * k - 1)) * x * LegendrePoly(x, k - 1) - (k - 1) * LegendrePoly(x, k - 2) / sqrt((double)(2 * k - 3))) / (double)k;
    }
  }
  

}
