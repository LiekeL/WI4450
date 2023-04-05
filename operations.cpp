#include "operations.hpp"
#include <omp.h>

void init(int n, double* x, double value)
{
  #pragma omp parallel for
  for (int i=0; i<n; i++){
    x[i]=value;
  }

  return;
}

double dot(int n, double const* x, double const* y)
{
  double dot_product = 0.0;
  #pragma omp parallel for reduction(+: dot_product)
  for (int i=0; i<n; i++){
    dot_product += x[i]*y[i];
  }

  return dot_product;
}

void axpby(int n, double a, double const* x, double b, double* y)
{
  #pragma omp parallel for
  for (int i=0; i<n; i++){
    y[i] = a*x[i] + b*y[i];
  }

  return;
}

//! apply a 7-point stencil to a vector
void apply_stencil3d(stencil3d const* S,
        double const* u, double* v)
{
  init((S->nx)*(S->ny)*(S->nz), v, 0.0);

  #pragma omp parallel for collapse(3)
  for (int k=0; k<S->nz; k++){
    for (int j=0; j<S->ny; j++){
        for (int i=0; i<S->nx; i++){
        
        v[S->index_c(i,j,k)] = S->value_c * u[S->index_c(i,j,k)];

        if (i != 0)       {v[S->index_c(i,j,k)] += S->value_w * u[S->index_w(i,j,k)];}
        if (i != S->nx-1) {v[S->index_c(i,j,k)] += S->value_e * u[S->index_e(i,j,k)];}
        if (j != 0)       {v[S->index_c(i,j,k)] += S->value_s * u[S->index_s(i,j,k)];}
        if (j != S->ny-1) {v[S->index_c(i,j,k)] += S->value_n * u[S->index_n(i,j,k)];}
        if (k != 0)       {v[S->index_c(i,j,k)] += S->value_b * u[S->index_b(i,j,k)];}
        if (k != S->nz-1) {v[S->index_c(i,j,k)] += S->value_t * u[S->index_t(i,j,k)];} 
        
      }
    }
  }

  return;
}

