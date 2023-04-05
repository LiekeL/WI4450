
#include "pcg_solver.hpp"
#include "operations.hpp"

#include <cmath>
#include <stdexcept>

#include <iostream>
#include <iomanip>
#include <omp.h>

void pcg_solver(stencil3d const* op, int n, double* x, double const* b,
        double tol, int maxIter,
        double* resNorm, int* numIter,
        int verbose)
{
  if (op->nx * op->ny * op->nz != n)
  {
    throw std::runtime_error("mismatch between stencil and vector dimension passed to cg_solver");
  }
  
  int num_procs = omp_get_max_threads(); 

  double *p = new double[n];
  double *q = new double[n];
  double *r = new double[n];

  double alpha, beta, rho=1.0, rho_old=0.0;

  // r = op * x
  apply_stencil3d(op, x, r);
  
  // r = b - r;
  axpby(n, 1.0, b, -1.0, r);


  // p = q = 0
  init(n, p, 0.0);
  init(n, q, 0.0);

  // start CG iteration

  #pragma omp parallel for
  for (int proc=0; proc<num_procs; proc++){
    #pragma omp task
       {
          #pragma omp critical
          {
            int iter = -1;
            while (true)
            {
              iter++;

              // rho = <r, r>
              rho = dot(n,r,r);

              if (verbose)
              {
                std::cout << std::setw(4) << iter << "\t" << std::setw(8) << std::setprecision(4) << sqrt(rho) << std::endl;
                // std::cout << std::setprecision(4) << rho << ","; // if you want to plot it in python
              }

              // check for convergence or failure
              
              if ((std::sqrt(rho) < tol) || (iter > maxIter))
              {
                break;
              }

              if (rho_old==0.0)
              {
                alpha = 0.0;
              }
              else
              {
                alpha = rho / rho_old;
              }
              // p = r + alpha * p
              axpby(n, 1.0, r, alpha, p);

              // q = op * p
              apply_stencil3d(op, p, q);

              // beta = <p,q>
              beta = dot(n, p, q);

              alpha = rho / beta;

              // x = x + alpha * p
              axpby(n, alpha, p, 1.0, x);

              // r = r - alpha * q
              axpby(n, -alpha, q, 1.0, r);

              std::swap(rho_old, rho);
            }// end of while-loop
            
            *numIter = iter;
          }
       }
  }
  #pragma omp taskwait

  // return number of iterations and achieved residual
  *resNorm = rho;

  // clean up
  delete [] p;
  delete [] q;
  delete [] r;

  return;
}
