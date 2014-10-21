#include <stdlib.h>
#include <R.h>
#include <math.h>


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
      for (int j=0; j<*n; j++)
      {
        for (int i=0; i<*n; i++)
        {
          for (int k=0; k<*m; k++)
          {
            G[j * *n + i] += X[i * *n + k] * X[j * *n + k];
          }
        }
      }
    }
    
    case 1: // polynomial kernel //
    {
      
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
    }
    
    case 3: // sigmoid kernel //
    {
    }
    
      
  }
}
