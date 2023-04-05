#include "operations.hpp"
#include "cg_solver.hpp"
#include "timer.hpp"

#include <iostream>
#include <cmath>
#include <limits>

#include <cmath>

// Forcing term
double f(double x, double y, double z)
{
  return z*sin(2*M_PI*x)*std::sin(M_PI*y) + 8*z*z*z;
}

// boundary condition at z=0
double g_0(double x, double y)
{
  return x*(1.0-x)*y*(1-y);
}

stencil3d laplace3d_stencil(int nx, int ny, int nz)
{
  if (nx<=2 || ny<=2 || nz<=2) throw std::runtime_error("need at least two grid points in each direction to implement boundary conditions.");
  stencil3d L;
  L.nx=nx; L.ny=ny; L.nz=nz;
  double dx=1.0/(nx-1), dy=1.0/(ny-1), dz=1.0/(nz-1);
  L.value_c = 2.0/(dx*dx) + 2.0/(dy*dy) + 2.0/(dz*dz);
  L.value_n =  -1.0/(dy*dy);
  L.value_e = -1.0/(dx*dx);
  L.value_s =  -1.0/(dy*dy);
  L.value_w = -1.0/(dx*dx);
  L.value_t = -1.0/(dz*dz);
  L.value_b = -1.0/(dz*dz);
  return L;
}

int main(int argc, char* argv[]){
    
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

    double *x = new double[n];
    double *r = new double[n];

    int n_runs = 50;

    for (int i=0; i<n_runs; i++){
    
    // intialize x = i with init
    {Timer t("init");
    init(n, x, double(i));
    }

    // calculate r = L*x with apply_stencil3d
    {Timer t("stencil");
    apply_stencil3d(&L, x, r);
    }

    // dot product of r and x with dot
    {Timer t("dot");
    dot(n, x, r);
    }

    // add 2x+3y with axpby 
    {Timer t("axpby");
    axpby(n, 2.0, x, -3.0, r);
    }
    }

    delete [] x;
    delete [] r;
    }
    Timer::summarize();
    return 0;
}