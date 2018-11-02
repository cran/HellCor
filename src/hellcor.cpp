#include <R.h>
#include "Rmath.h"

extern "C" {

  void hellcorC(double *x, int *xlen, double *statistic, int *pvalcomp, double *pvalue, int *KLmaxvalue, double *alphavalue) {

    int n = xlen[0];
    
    int N = n / 2;
  
    if (N>3) {

      double tmp, tmp1, tmp2, tmp3, min, cte, cte2, etahat, *toto, *titi, *tata, *tutu;
      double *u1, *u2, *u1ties, *u2ties, *s1, *s2, **Legtmp1, **Legtmp2, *weight, *R, *R2, **Rstar, **dist, **betahat, **Aterm, **Aij, **Cterm;
      int KLmax = KLmaxvalue[0], k, l, i, j, i1, i2, **II, Khat = 0, Lhat = 0;
      double alpha = alphavalue[0];
      u1 = new double[N];
      u2 = new double[N];
      u1ties = new double[N];
      u2ties = new double[N];
      s1 = new double[N];
      s2 = new double[N];
      Legtmp1 = new double*[N];
      for (i = 0; i < N; i++) Legtmp1[i] = new double[KLmax];
      Legtmp2 = new double*[N];
      for (i = 0; i < N; i++) Legtmp2[i] = new double[KLmax];
      weight = new double[N];
      R = new double[N];
      R2 = new double[N];
      Rstar = new double*[N];
      for (i = 0; i < N; i++) Rstar[i] = new double[N];
      dist = new double*[N];
      for (i = 0; i < N; i++) dist[i] = new double[N];
      betahat = new double*[KLmax];
      for (i = 0; i < KLmax; i++) betahat[i] = new double[KLmax];
      Aterm = new double*[KLmax];
      for (i = 0; i < KLmax; i++) Aterm[i] = new double[KLmax];
      Aij = new double*[KLmax];
      for (i = 0; i < KLmax; i++) Aij[i] = new double[KLmax];
      Cterm = new double*[KLmax];
      for (i = 0; i < KLmax; i++) Cterm[i] = new double[KLmax];
      II = new int*[N];
      for (i = 0; i < N; i++) II[i] = new int[2];
      toto = new double[KLmax];
      titi = new double[KLmax];
      tata = new double[KLmax];
      tutu = new double[KLmax];

      double LegendrePoly(double x, int k);
      double etafunc(double H2);
      
      for (i = 0; i < N; i++) {
	u1[i] = 0.0;
	u2[i] = 0.0;
	R[i] = R_PosInf;
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

      for (i = 0; i < N; i++) {
	u1[i] = u1ties[i] / double(N + 1);
	u2[i] = u2ties[i] / double(N + 1);
      }

      // (u1, u2) now contains the empirical copula, i.e., UU <- apply(XX, 2, rank) / (n+1)
      // with ties.method = "average"

      for (i = 0; i < N; i++) {
	s1[i] = Rf_qbeta(u1[i], alpha, alpha, 1, 0); 
	s2[i] = Rf_qbeta(u2[i], alpha, alpha, 1, 0);
	weight[i] = sqrt(Rf_dbeta(s1[i], alpha, alpha, 0) * Rf_dbeta(s2[i], alpha, alpha, 0));
      }
      
      for (i = 0; i < (N - 1); i++) {
	for (j = i + 1; j < N; j++) {
	  dist[i][j] = sqrt(R_pow(s1[i] - s1[j], 2.0) + R_pow(s2[i] - s2[j], 2.0));
	}
      }

      for (i = 0; i < N; i++) {
	R[i] = R_PosInf; // min1 = smallest
	R2[i] = R_PosInf; // min2 = second smallest (can be equal to min1)
	II[i][0] = 1;
	II[i][1] = 2;
	for (j = 0; j < N; j ++) {
	  if (j > i) tmp = dist[i][j]; else if (j < i) tmp = dist[j][i];
	  if (j != i) {
	    if (tmp < R[i]) {
	      R2[i] = R[i];
	      II[i][1] = II[i][0];
	      R[i] = tmp;
	      II[i][0] = j + 1;
	    }
	    if (tmp < R2[i]) {
	      R2[i] = tmp;
	      II[i][1] = j + 1;
	    }
	  }
	}
      }

      cte = 2.0 * sqrt((double)(N - 1)) / (double)N;
      cte2 = 4.0 * sqrt((double)(N - 2)) / (double)(N - 1);
      for (i = 0; i < N; i++) {
	R[i] = R[i] * weight[i];
	R2[i] = R2[i] * weight[i];
	for (j = 0; j < N; j++) {
	  Rstar[i][j] = cte2 * R[i];
	}
	//	if (II[i][0] < 1 && II[i][0] > N) Rprintf(" %g", II[i][0]);  // CHECK!!
	Rstar[i][II[i][0] - 1] = cte2 * R2[i];
	Rstar[i][i] = 0.0;
	R[i] = cte * R[i];
	//	Rprintf(" %g", R[i]);
      }

      for (i = 0; i < N; i++) {
	for (k = 0; k < KLmax; k++) {
	  Legtmp1[i][k] = sqrt(2.0) * LegendrePoly(2.0 * u1[i] - 1.0, k);
	  Legtmp2[i][k] = sqrt(2.0) * LegendrePoly(2.0 * u2[i] - 1.0, k);
	}
      }

      for (k = 0; k < KLmax; k++) {
 	for (l = 0; l < KLmax; l++) {
 	  tmp = 0.0;
 	  for (i = 0; i < N; i++) {
 	    tmp = tmp +  R[i] * Legtmp1[i][k] * Legtmp2[i][l];
 	  }
	  betahat[k][l] =  tmp;
 	}
      }

      // Aterm is KLmax x KLmax
      for (l = 0; l < KLmax; l++) {
	tmp = 0.0;
	for (k = 0; k < KLmax; k++) {
	  tmp = tmp + R_pow(betahat[k][l], 2.0);
	  Aterm[k][l] = tmp;
	}
      }
      for (k = 0; k < KLmax; k++) {
	tmp = 0.0;
	for (l = 0; l < KLmax; l++) {
	  tmp = tmp + Aterm[k][l];
	  Aterm[k][l] = tmp;
	  Aij[k][l] = 0.0;
	}
      }

      for (i2 = 0; i2 < N; i2++) {
      	for (l = 0; l < KLmax; l++) {
      	  toto[l] = Legtmp1[i2][l];
      	  titi[l] = Legtmp2[i2][l];
      	}
      	tmp1 = R[i2];
      	for (i1 = 0; i1 < N; i1++) {
      	  tmp2 = R[i1] * (tmp1 - Rstar[i2][i1]);
      	  for (l = 0; l < KLmax; l++) {
      	    tata[l] = tmp2 * toto[l] * Legtmp1[i1][l];
      	    tutu[l] = titi[l] * Legtmp2[i1][l];
      	  }
      	  for (k = 0; k < KLmax; k++) {
      	    tmp3 = tata[k];
      	    for (l = 0; l < KLmax; l++) {
      	      Aij[k][l] =  Aij[k][l] + tmp3 * tutu[l] ;
      	    }
      	  }
      	}
      }

      // Cterm is KLmax x KLmax
      for (l = 0; l < KLmax; l++) {
	tmp = 0.0;
	for (k = 0; k < KLmax; k++) {
	  tmp = tmp + Aij[k][l];
	  Cterm[k][l] = tmp;
	}
      }
      for (k = 0; k < KLmax; k++) {
	tmp = 0.0;
	for (l = 0; l < KLmax; l++) {
	  tmp = tmp + Cterm[k][l];
	  Cterm[k][l] = tmp;
	}
      }

      min = Cterm[0][0] / sqrt(Aterm[0][0]);
      for (k = 0; k < KLmax; k++) {
	for (l = 0; l < KLmax; l++) {
	  tmp = Cterm[k][l] / sqrt(Aterm[k][l]);
	  if (tmp < min) {
	    min = tmp;
	    Khat = k;
	    Lhat = l;
	  }
	}
      }
      if (Khat < 1) Khat = 1;
      if (Lhat < 1) Lhat = 1;

      etahat = etafunc(1.0 - betahat[0][0] / sqrt(Aterm[Khat][Lhat]));
      
      statistic[0] = etahat; // Here is the test statistic value

      if (pvalcomp[0] == 0) {
      
	pvalue[0] = 0.0;

      } else {

	// TODO
	
      }
      
 // We free the unused array of pointers
      for (i = 0; i < N; i++) delete[] Legtmp1[i];
      delete[] Legtmp1;
      for (i = 0; i < N; i++) delete[] Legtmp2[i];
      delete[] Legtmp2;
      delete[] u1;
      delete[] u2;
      delete[] u1ties;
      delete[] u2ties;
      delete[] s1;
      delete[] s2;
      delete[] weight;
      delete[] R;
      delete[] R2;
      for (i = 0; i < N; i++) delete[] dist[i];
      delete[] dist;
      for (i = 0; i < N; i++) delete[] Rstar[i];
      delete[] Rstar;
      for (i = 0; i < KLmax; i++) delete[] betahat[i];
      delete[] betahat;
      for (i = 0; i < KLmax; i++) delete[] Aterm[i];
      delete[] Aterm;
      for (i = 0; i < KLmax; i++) delete[] Aij[i];
      delete[] Aij;
      for (i = 0; i < KLmax; i++) delete[] Cterm[i];
      delete[] Cterm;
      for (i = 0; i < N; i++) delete[] II[i];
      delete[] II;
      delete[] toto;
      delete[] titi;
      delete[] tata;
      delete[] tutu;
	}

// We return
    return;
          
  }

  
double etafunc(double H2) {
  double B = 1.0 - H2;
  double B2 = B * B;
  double B4 = B2 * B2;
  return 2.0 * sqrt(-2.0 + B4 + sqrt(4.0 - 3.0 * B4)) / B2;
}

  double LegendrePoly(double x, int k) {
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
      return 0.0;
    }
  }
  

}
