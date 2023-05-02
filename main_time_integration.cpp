#include "operations.hpp"
#include "time_integration.hpp"
#include "timer.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <limits>

#include <cmath>

// Main program that solves the time integration in 3D
// on a unit cube. The grid size (nx,ny,nz) can be 
// passed to the executable like this:
//
// ./main_time_integration <nx> <ny> <nz>
//
// or simply ./main_time_integration <nx> for ny=nz=nx.
// If no arguments are given, the default nx=ny=nz=128 is used.
//
// Boundary conditions and forcing term f(x,y,z) are
// hard-coded in this file. See README.md for details
// on the PDE and boundary conditions.

stencil3d laplace3d_stencil(int nx, int ny, int nz)
{
  if (nx<=2 || ny<=2 || nz<=2) throw std::runtime_error("need at least two grid points in each direction to implement boundary conditions.");
  stencil3d L;
  L.nx=nx; L.ny=ny; L.nz=nz;
  double dx=1.0/(nx-1), dy=1.0/(ny-1), dz=1.0/(nz-1);
  L.value_c = 2.0/(dx*dx) + 2.0/(dy*dy) + 2.0/(dz*dz);
  L.value_n = -1.0/(dy*dy);
  L.value_e = -1.0/(dx*dx);
  L.value_s = -1.0/(dy*dy);
  L.value_w = -1.0/(dx*dx);
  L.value_t = -1.0/(dz*dz);
  L.value_b = -1.0/(dz*dz);
  return L;
}

int main(int argc, char* argv[])
{
  {Timer t("main");
  int nx, ny, nz;

  if      (argc==1) {nx=128;           ny=128;           nz=128;}
  else if (argc==2) {nx=atoi(argv[1]); ny=nx;            nz=nx;}
  else if (argc==4) {nx=atoi(argv[1]); ny=atoi(argv[2]); nz=atoi(argv[3]);}
  else {std::cerr << "Invalid number of arguments (should be 0, 1 or 3)"<<std::endl; exit(-1);}
  if (ny<0) ny=nx;
  if (nz<0) nz=nx;

  // total number of unknowns
  int n=nx*ny*nz;

  double dx=1.0/(nx-1), dy=1.0/(ny-1), dz=1.0/(nz-1);

  // Laplace operator
  stencil3d L = laplace3d_stencil(nx,ny,nz);

  // initial value: initial value for the time integration method included in the rhs
  double *b = new double[n];
  init(n, b, 1.0);

  // solve the linear system of equations using parallel forward euler
  int numIter=0, maxIter=40;
  double resNorm=10e6, tol=std::sqrt(std::numeric_limits<double>::epsilon());
  double delta_t = 10e-6;

  // solution vector: start with a 0 vector
  double *x = new double[n*maxIter];
  init(n, x, 0.0);

  try {
  Timer t("time_integration");
  time_integration(&L, n, x, b, tol, delta_t, maxIter, &resNorm, &numIter, 1.0);
  std::cout << std::setw(4) << maxIter << "\t" << std::setw(8) << std::setprecision(4) << resNorm << std::endl;
  } catch(std::exception e)
  {
    std::cerr << "Caught an exception in time_integation: " << e.what() << std::endl;
    exit(-1);
  }
  delete [] x;
  }
  Timer::summarize();

  return 0;
}
