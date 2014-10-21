#include <stdlib.h>
#include <R.h>
#include <math.h>
#include <R_ext/Lapack.h>

void ldm_train( double *x, 
                double *y, 
                int *N, 
                int *M,
                int *kernel_type, 
                double *tuning_params,
                double *G,
                double *Gy,
                double *GY,
                double *Q,
                int *IPIV,
                double *WORK,
                double *A,
                double *h,
                double *alpha,
                double *beta,
                double *beta_old,
                double *del_beta,
                double *result)
                
{ 
  Rprintf("Line 10");







  double cost = tuning_params[0];
  double lambda_1 = tuning_params[1];
  double lambda_2 = tuning_params[2];
  int n = *N;
  int m = *M;
  
  Rprintf("Line 24\n");
  
  // make gram matrix //
  
  Rprintf("Line 30\n");
  switch(*kernel_type)
  {
    case 0: // linear kernel //
    {
      for (int j=0; j<n; j++)
      {
        for (int i=0; i<n; i++)
        {
          for (int k=0; k<n; k++)
          {
            G[j * n + i] += x[i * n + k] * x[j * n + k];
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
      for (int j=0;j<n;j++)
      {
        for (int i=0;i<n;i++) 
        {
          //  set G[i,j] //
          if (i==j) G[j * n + i] = 1;
          else if (j > i) G[j * n + i] = G[i * n + j];
          else {
            double tmp = 0;
            for (int k=0;k<m;k++) // ||Xi - Xj||_2
            {
              tmp += pow(x[i + k * n] - x[j + k * n], 2);
            }
            G[j * n + i] =  exp(-rbf_gamma * tmp);
          }
        }
      }
    }
    
    case 3: // sigmoid kernel //
    {
    }
    
      
  }
  Rprintf("Line 80\n");
  // set Gy //
  for (int i=0;i<n;i++)
  {
    Gy[i] = 0;
    for(int k=0;k<n;k++)
    {
      Gy[i] += G[k*n + i] * y[k];
    }
  }
  
  // GY //
  Rprintf("Line 93\n");
  for (int j=0;j<n;j++)
  {
    for(int i=0;i<n;i++)
    {
      GY[j * n + i] = G[j * n + i] * y[j];
    }
  }
  Rprintf("Line 102\n");
  // set Q //
  for (int j=0;j<n;j++)
  {
    for (int i=0;i<n;i++)
    {
      if (j > i) Q[j * n + i] = Q[i * n + j]; // symmetry of Q //
      else
      {
        Q[j * n + i] = 0;
        for (int k=0;k<n;k++)
        {
          Q[j * n + i] +=  G[i * n + k] * G[j * n + k];
        }
        Q[j * n + i] /= n;
        Q[j * n + i] -= (Gy[i] * Gy[j])/(n * n);
        Q[j * n + i] *= 4 * lambda_1;
        Q[j * n + i] += G[j * n + i];
      }
    }
  }
  Rprintf("Line 124\n");
  // Q -> Q^-1
  int INFO;
  int LWORK = n * n;
  F77_CALL(dgetrf)(&n, &n, Q, &n, IPIV, &INFO); // LU factorization of Q //
  F77_CALL(dgetri)(&n, Q, &n, IPIV, WORK, &LWORK, &INFO); // Inverse of Q //
  Rprintf("Line 132\n");
  // set A
  for (int j=0; j<n; j++)
  {
    for (int i=0; i<n; i++)
    {
      A[j * n + i] = 0;
      for (int k=0; k<n; k++)
      {
        A[j * n + i] += Q[k * n + i] * GY[j * n + k];
      }
    }
  }
  Rprintf("Line 146\n");
  // set h
  for (int i=0; i<n; i++)
  {
    h[i] = 0;
    for (int k=0; k<n; k++)
      {
        h[i * n + i] += A[i * n + k] * GY[i * n + k];
      }
  }
  Rprintf("Line 157\n");
  for (int i=0; i<n; i++)
  {
    alpha[i] = 0;
    for (int k=0; k<n; k++)
    {
      alpha[i] += Q[k * n + i] * Gy[k];
    }
    alpha[i] *= (lambda_2 / n);
  }
  Rprintf("Line 173\n");
  double delta_beta = 0;
  Rprintf("Line 177\n");
  
  int counter = 0;
  do
  {
    for (int i=0; i<n; i++)
    {
      del_beta[i] = 0; // update gradient of beta
      for (int k=0; k<n; k++)
      {
        del_beta[i] += GY[i * n + k] * alpha[k];
      }
      del_beta[i] -= 1;
      beta_old[i] = beta[i];
      beta[i] -= del_beta[i]/h[i];
      if (beta[i] <= 0) beta[i] = 0;
      else if (beta[i] >= cost) beta[i] = cost;
      for (int k=0; k<n; k++)
      {
        alpha[k] += (beta[i] - beta_old[i]) * A[i*n +k];
      }
    }
    
    delta_beta = 0;
    for (int i=0; i<n; i++)
    {
      delta_beta += pow((beta[i] - beta_old[i]), 2);
    }
    delta_beta = sqrt(delta_beta);
    
    counter++;
    if (counter % 50 == 0) {
      Rprintf("Training iterations: %d \n", counter);
      Rprintf("DELTA BETA: %d \n", delta_beta);
    }
    if (counter > 1000) break;
  } while (delta_beta > 1.0e-3);
  Rprintf("Line 203");
  
  // final result
  
  for (int i=0; i<n; i++)
  {
    result[i] = 0;
    for (int k=0; k<n; k++)
    {
      result[i] += A[k * n + i] * (lambda_2 / n + beta[k]);
    }
  }
  Rprintf("Line 215\n");
}

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
      if (beta[i] <= 0) beta[i] = 0;
      else if (beta[i] >= *cost) beta[i] = *cost;
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
    if (counter % 50 == 0) {
      Rprintf("Training iterations: %d \n", counter);
      Rprintf("DELTA BETA: %f \n", delta_beta);
    }
    if (counter > 1000) break;
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
