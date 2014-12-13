#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <math.h>
#include <time.h>

int rand_lim1(int limit) {
    int divisor = RAND_MAX/(limit+1);
    int ret;
    while(1) {
      ret = rand() / divisor;
      if(ret <= limit) return ret;
    }
}

SEXP avSGD2(SEXP X, SEXP Y, SEXP dim, SEXP params) {
  
  srand(time(NULL));
  int n, m;
  dim = coerceVector(dim, INTSXP);
  n = INTEGER(dim)[0];
  m = INTEGER(dim)[1];
  //Rprintf("n = %d, m = %d \n", n, m);
  
  double *x, *y;
  x = REAL(X);
  y = REAL(Y);
  
  double C, eta, lambda_1, lambda_2;
  params = coerceVector(params, REALSXP);
  C = REAL(params)[0];
  eta = REAL(params)[1];
  lambda_1 = REAL(params)[2];
  lambda_2 = REAL(params)[3];
  
  Rprintf("C = %f \n", C);
  Rprintf("eta = %f \n", eta);
  Rprintf("lambda_1= %f \n", lambda_1);
  Rprintf("lambda_2 = %f \n", lambda_2);
  Rprintf("length of X: %d \n", length(X));
  
  //Rprintf("X , y \n");
  //for(int i=0; i < n; i++) {
  //  Rprintf("%f \t %f \t %f \n", x[i], x[i + n], y[i]);
  //}
  
  SEXP W = PROTECT(allocVector(REALSXP, m));
  SEXP WBAR = PROTECT(allocVector(REALSXP, m));
  SEXP GRADIENT = PROTECT(allocVector(REALSXP, m));
  
  double *w = REAL(W);
  double *w_bar = REAL(WBAR);
  double *gradient = REAL(GRADIENT);
  
  for(int t=0; t <= n * 5; t++) {
    int i = rand_lim1(n - 1);
    int j;
    do { j = rand_lim1(n - 1); } while(i == j);
    
    double idot = 0;
    for(int p=0; p < m; p++) {
        idot += x[i + p * n] * w[p];
    }
    
    double jdot = 0;
      for(int p=0; p < m; p++) {
        jdot += x[j + p * n] * w[p];
    }
    for(int k=0; k < m; k++) {
      gradient[k] = 4 * lambda_1 * idot * x[i + k * n];
      gradient[k] -= 4 * lambda_1 * jdot * y[i] * y[j] * x[i + k * n];
      gradient[k] += w[k];
      gradient[k] -= lambda_2 * y[i] * x[i + k * n];
      gradient[k] -= C * n * y[i] * x[i + k * n] * ((y[i] * idot < 1) ? 1 : 0);
    }
    for(int k=0; k < m; k++) {
      w[k] -= eta * gradient[k];
      w_bar[k] += (1 / t) * (w[k] - w_bar[k]); 
    }
  
  }
  
    UNPROTECT(3);
  return(WBAR);
 
}