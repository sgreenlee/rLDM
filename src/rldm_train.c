#include <stdlib.h>
#include <R.h>
#include <math.h>
#include <time.h>


void coordDescent(int *n, double *A, double *GY, double *alpha, double *beta,
                    double *beta_old, double *del_beta, double *h, double *cost)
{
  double const EPSILON = 1e-6;
  double delta_beta;
  int counter = 0;
  do
  {
    for (int i=0; i<*n; i++)
    {
      del_beta[i] = 0; // update gradient of beta
      for (int k=0; k<*n; k++)
      {
        del_beta[i] += GY[i * *n + k] * alpha[k];
      }
      del_beta[i] -= 1;
      beta_old[i] = beta[i];
      beta[i] -= del_beta[i]/h[i];
      beta[i] = beta[i] < *cost ? (beta[i] > 0 ? beta[i] : 0) : *cost;
      for (int k=0; k<*n; k++)
      {
        alpha[k] += (beta[i] - beta_old[i]) * A[i * *n + k];
      }
    }
    
    delta_beta = 0;
    for (int i=0; i<*n; i++)
    {
      delta_beta += pow((beta[i] - beta_old[i]), 2);
    }
    delta_beta = sqrt(delta_beta);
    
    counter++;
    if (counter > 5000) break;
  } while (delta_beta > EPSILON);
}


void makeGramMatrix(double *X, int *n, int* m, double*G, 
                      int *kernel_type, double *tuning_params) 
{  
  switch(*kernel_type)
  {
    case 0: // linear kernel //
    {
      Rprintf("Kernel type: linear");
      for (int j=0;j<*n;j++)
      {
        for (int i=0;i<*n;i++) 
        {
          //  set G[i,j] //
          double tmp = 0;
          for (int k=0;k<*m;k++)
          {
            tmp += X[i + k * *n] * X[j + k * *n];
          }
          G[j * *n + i] =  tmp;
          }
        }
      break;
    }
    
    case 1: // polynomial kernel //
    {
      Rprintf("Kernel type: polynomial");
      double degree = tuning_params[3];
      double c = tuning_params[2];
      for (int j=0;j<*n;j++)
      {
        for (int i=0;i<*n;i++) 
        {
          //  set G[i,j] //
          double tmp = 0;
          for (int k=0;k<*m;k++)
          {
            tmp += X[i + k * *n] * X[j + k * *n];
          }
          G[j * *n + i] =  pow(tmp + c, degree);
        }
      }
      break;
    }
    
    case 2: // radial-basis-function kernel //
    {
      double rbf_gamma = tuning_params[4];
      for (int j=0;j<*n;j++)
      {
        for (int i=0;i<*n;i++) 
        {
          //  set G[i,j] //
          if (i==j) G[j * *n + i] = 1;
          else if (j > i) G[j * *n + i] = G[i * *n + j];
          else {
            double tmp = 0;
            for (int k=0;k<*m;k++) // ||Xi - Xj||_2
            {
              tmp += pow(X[i + k * *n] - X[j + k * *n], 2);
            }
            G[j * *n + i] =  exp(-rbf_gamma * tmp);
          }
        }
      }
      break;
    }
    
    case 3: // sigmoid kernel //
    {
      Rprintf("Kernel type: sigmoid kernel");
      double a = tuning_params[1];
      double r = tuning_params[0];
      for (int j=0;j<*n;j++)
      {
        for (int i=0;i<*n;i++) 
        {
          double tmp = 0;
          for (int k=0;k<*m;k++) // ||Xi - Xj||_2
          {
            tmp += X[i + k * *n] * X[j + k * *n];
          }
          G[j * *n + i] =  tanh(a * tmp + r);
        }
      }
      break;
    } 
  }
}

void ldmPredict(double *alpha, double *rbf_gamma, double *newdata, int *dim, 
                  double *modelMatrix, double *pred)
{
  int n = dim[0]; // # of new instances to predict
  int m = dim[1]; // # of instances in training set
  int p = dim[2]; // # of features
  
  
  double tmp;
  for (int i=0;i<n;i++)
  {
    pred[i] = 0;
    for (int j=0;j<m;j++) 
    {
      tmp = 0;
      for (int k=0;k<p;k++)
      {
        tmp += pow(newdata[i + k * p] - modelMatrix[j + k * p], 2);
      }
      tmp = alpha[j] * exp(-*rbf_gamma * tmp);
      pred[i] += tmp;
    }
  }
}

int rand_lim(int limit) {
    int divisor = RAND_MAX/(limit+1);
    int ret;
    while(1) {
      ret = rand() / divisor;
      if(ret <= limit) return ret;
    }
}

void avSGD(int *n, int *m, double *X, double *y, double *C, double *eta, double *w,
            double *w_bar, double *gradient, double *lambda_1, double *lambda_2) {
              
  Rprintf("Address of n: %p \n", n);
  Rprintf("Address of m: %p \n", m);
  Rprintf("Address of X: %p \n", X);
  Rprintf("Address of y: %p \n", y);
  Rprintf("Address of C: %p \n", C);
  Rprintf("Address of eta: %p \n", eta);
  Rprintf("Address of w: %p \n", w);
  Rprintf("Address of w_bar: %p \n", w_bar);
  Rprintf("Address of gradient: %p \n", gradient);
  Rprintf("Address of lambda_1: %p \n", lambda_1);
  Rprintf("Address of lambda_2: %p \n", lambda_2);
              
  srand(time(NULL));
  
  for(int t=0; t <= *n * 5; t++) {
    int i = rand_lim(*n - 1);
    int j;
    do { j = rand_lim(*n - 1); } while(i == j);
    
    double idot = 0;
    for(int p=0; p < *m; p++) {
        idot += X[i + p * *n] * w[p];
    }
    int C = (y[i] * idot < 1) ? 1 : 0;
    double jdot = 0;
      for(int p=0; p < *m; p++) {
        jdot += X[j + p * *n] * w[p];
    }
    for(int k=0; k < *m; k++) {
      gradient[k] = 4 * *lambda_1 * idot * X[i + k * *n];
      gradient[k] -= 4 * *lambda_1 * jdot * y[i] * y[j] * X[i + k * *n];
      gradient[k] += w[k];
      gradient[k] -= *lambda_2 * y[i] * X[i + k * *n];
      gradient[k] -= C * *n * y[i] * X[i + k * *n];
    }
    for(int k=0; k < *m; i++) {
      w[k] -= *eta * gradient[k];
      w_bar[k] += (1 / t) * (w[k] - w_bar[k]); 
    }
  
  }
}
