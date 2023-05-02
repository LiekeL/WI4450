
#include "time_integration.hpp"
#include "operations.hpp"

#include <cmath>
#include <stdexcept>
#include <algorithm>

#include <iostream>
#include <iomanip>
#include <omp.h>

void time_integration(stencil3d const* op, int n, double* x, double const* b,
        double tol, double delta_t, int maxIter,
        double* resNorm, int* numIter,
        int verbose)
{
  if (op->nx * op->ny * op->nz != n)
  {
    throw std::runtime_error("mismatch between stencil and vector dimension passed to cg_solver");
  }
  
  
  // make the rhs from b and zero's
  double *rhs = new double[n*maxIter];
  double *rhs_0 = new double[n*(maxIter-1)];
  init(n*(maxIter-1), rhs_0, 0.0);
  std::copy(b, b + n, rhs);
  std::copy(rhs_0, rhs_0 + n*(maxIter-1), rhs + n);
  delete [] rhs_0;

  double *r = new double[n*maxIter];
  double *Ax = new double[n*maxIter];

  double *x_j1 = new double[n];
  double *x_j = new double[n];
  double *Ax_j1 = new double[n];

  double rho=1.0;

  // start at x = rhs
  axpby(n*maxIter, 1.0, rhs, 0.0, x);

  for (int iter=0; iter<maxIter; iter++){

    // calculate Ax (TODO: could be more efficient)
    for (int j=0; j<maxIter; j++){

      // copy x[j*n to (j+1)*n] into x_j
      for(int l=j*n; l<(j+1)*n; l++){
          x_j[l-j*n] = x[l];
        }

      if (j==0){
        init(n, x_j1, 0);
      } else {
        // copy x[(j-1)*n to j*n] into x_j1
        for(int l=(j-1)*n; l<j*n; l++){
          x_j1[l-(j-1)*n] = x[l];
        }

        // calculate (-I+delta_t*A)*x_{j-1} = -(x_j1)+delta_t*A(x_j1)
        // Ax_j1 = op*x_j1
        apply_stencil3d(op, x_j1, Ax_j1);
        // x_j1 = - x_j1 + delta_t*Ax_j1
        axpby(n, delta_t, Ax_j1, -1.0, x_j1);
      }
      
      // x_j = x_j + x_j1 (x_j1 = -I + delta*L*x_{j-1})
      axpby(n, 1.0, x_j1, 1.0, x_j);

      // copy x_j into Ax
      for(int l=j*n; l<(j+1)*n; l++){
          Ax[l] = x_j[l-j*n];
      }
    }

    // r = rhs - Ax
    axpby(n*maxIter, 1.0, rhs, 0.0, r); // first copy rhs to r
    axpby(n*maxIter, -1.0, Ax, 1.0, r);


    // calculate norm of r and print it
    rho = dot(n*maxIter,r,r);

    if (verbose)
    {
      std::cout << std::setw(4) << iter << "\t" << std::setw(8) << std::setprecision(4) << rho << std::endl;
      // std::cout << std::setprecision(4) << rho << ","; // if you want to plot it in python
    }

    // x = x + r
    axpby(n*maxIter, 1.0, r, 1.0, x);

  }

  delete [] rhs;
  delete [] r;
  delete [] Ax;
  delete [] x_j1;
  delete [] x_j;
  delete [] Ax_j1;

  *resNorm = rho;

  return;
}















// int num_procs = omp_get_max_threads(); 
  
//   maxIter = num_procs * maxIter/num_procs;
//   int no_blocks = num_procs;
//   printf("no_blocks = %i\n", no_blocks);

//   double *x_p = new double[n*maxIter/no_blocks];

//   // initialize rhs and x
//   double *x_0 = new double[n*maxIter];
//   double *Ax = new double[n*(maxIter/no_blocks)];
//   double *zeros = new double[n*(maxIter/no_blocks-1)];

//   // use that x_0 = rhs and use copy to make the rhs from b and zero's
//   double *x_0_0 = new double[n*(maxIter-1)];
//   init(n*(maxIter-1), x_0_0, 0.0);
//   std::copy(b, b + n, x_0);
//   std::copy(x_0_0, x_0_0 + n*(maxIter-1), x_0 + n);

//   // x = x_0 (x = 0*x + 1*x_0)
//   axpby(n*maxIter, 1.0, x_0, 0, x);

//   for(int iter;iter<maxIter;iter++){
//     for(int proc; proc<no_blocks;proc++){

//       // initialize part of solution vector x
//       init(n, x_p, 0.0);

//       // communicate necessary part of x_i to new process TODO
//       double *x_i = new double[n];
//       #pragma omp parallel for
//       for (int i=n*(proc*maxIter/no_blocks-1); i<n*proc*maxIter/no_blocks; i++){
//         x_i[i-n*(proc*maxIter/no_blocks-1)] = x[i];
//       } 

//       // calculate (A-M)*x_i
//       double *new_x_i = new double[n];
//       // A*x_i
//       apply_stencil3d(op, x_i, new_x_i);
//       // new_x_i = delta_t*new_x_i - I*x_i
//       axpby(n, -1, x_i, delta_t, new_x_i);

//       // make (A-M)*x_i complete
//       init(n*(maxIter/no_blocks-1), zeros, 0.0);
//       std::copy(new_x_i, new_x_i + n, Ax);
//       std::copy(zeros, zeros + n*(maxIter/no_blocks-1), Ax + n);

//       // x_p = b - (A-M)*x_i
//       #pragma omp parallel for
//       for (int i=proc*(n*maxIter/no_blocks); i<(proc+1)*(n*maxIter/no_blocks); i++){
//         x_p[i-proc*n*maxIter/no_blocks] = x_0[i]-Ax[i];
//       }

//       // solve M * x_(i+1) = x_p using forward substitution
      
//     }
//   }
//   // clean up
//   delete [] x_0;
//   delete [] x_0_0;
//   delete [] x_p;
