#include "gtest_mpi.hpp"

#include "operations.hpp"

#include <iostream>

// note: you may add any number of tests to verify
// your code behaves correctly, but do not change
// the existing tests.

TEST(stencil, bounds_check)
{
  stencil3d S;
  S.nx=5;
  S.ny=3;
  S.nz=2;
  EXPECT_THROW(S.index_c(-1,0,0), std::runtime_error);
  EXPECT_THROW(S.index_c(S.nx,0,0), std::runtime_error);
  EXPECT_THROW(S.index_c(0,-1,0), std::runtime_error);
  EXPECT_THROW(S.index_c(0,S.ny,0), std::runtime_error);
  EXPECT_THROW(S.index_c(0,0,-1), std::runtime_error);
  EXPECT_THROW(S.index_c(0,0,S.nz), std::runtime_error);
}

TEST(stencil, index_order_kji)
{
  stencil3d S;
  S.nx=50;
  S.ny=33;
  S.nz=21;

  int i=10, j=15, k=9;

  EXPECT_EQ(S.index_c(i,j,k), S.index_c(i-1,j,k)+1);
  EXPECT_EQ(S.index_c(i,j,k), S.index_c(i,j-1,k)+S.nx);
  EXPECT_EQ(S.index_c(i,j,k), S.index_c(i,j,k-1)+S.nx*S.ny);
}

TEST(operations, init)
{
  const int n=15;
  double x[n];
  for (int i=0; i<n; i++) x[i]=double(i+1);

  double val=42.0;
  init(n, x, val);

  double err=0.0;
  for (int i=0; i<n; i++) err = std::max(err, std::abs(x[i]-val));

  // note: EXPECT_NEAR uses a tolerance relative to the size of the target,
  // near 0 this is very small, so we use an absolute test instead by 
  // comparing to 1 instead of 0.
  EXPECT_NEAR(1.0+err, 1.0, std::numeric_limits<double>::epsilon());
}


TEST(operations, dot) {
  const int n=150;
  double x[n], y[n];

  for (int i=0; i<n; i++)
  {
    x[i] = double(i+1);
    y[i] = 1.0/double(i+1);
  }

  double res = dot(n, x, y);
  EXPECT_NEAR(res, (double)n, n*std::numeric_limits<double>::epsilon());
}

TEST(operations, axpby) {
  const int n=150;
  double x[n], y[n];

  for (int i=0; i<n; i++)
  {
    x[i] = double(i+1);
    y[i] = double(n-i-1);
  }

  double a = 42.0;
  double b = a;
  axpby(n, a, x, b, y);

  double err=0.0;

  for (int i=0; i<n; i++) err = std::max(err, std::abs(y[i]-a*n));

  EXPECT_NEAR(1.0+err, 1.0, std::numeric_limits<double>::epsilon());
}

TEST(operations,stencil3d_symmetric)
{
//  const int nx=3, ny=4, nz=5;
  const int nx=2, ny=2, nz=2;
  const int n=nx*ny*nz;
  double* e=new double[n];
  for (int i=0; i<n; i++) e[i]=0.0;
  double* A=new double[n*n];

  stencil3d S;

  S.nx=nx; S.ny=ny; S.nz=nz;
  S.value_c = 8;
  S.value_n = 2;
  S.value_e = 4;
  S.value_s = 2;
  S.value_w = 4;
  S.value_b = 1;
  S.value_t = 1;

  for (int i=0; i<n; i++)
  {
    e[i]=1.0;
    if (i>0) e[i-1]=0.0;
    apply_stencil3d(&S, e, A+i*n);
  }

  int wrong_entries=0;
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
    {
      if (A[i*n+j]!=A[j*n+i]) wrong_entries++;
    }
  EXPECT_EQ(0, wrong_entries);

  if (wrong_entries)
  {
    std::cout << "Your matrix (computed on a 2x2x2 grid by apply_stencil(I)) is ..."<<std::endl;
    for (int j=0; j<n; j++)
    {
      for (int i=0; i<n; i++)
      {
        std::cout << A[i*n+j] << " ";
      }
      std::cout << std::endl;
    }
  }
  delete [] e;
  delete [] A;
}

TEST(operations,stencil3d_diagdom)
{
  const int nx=3, ny=4, nz=5;
  const int n=nx*ny*nz;
  double* e=new double[n];
  for (int i=0; i<n; i++) e[i]=0.0;
  double* A=new double[n*n];

  stencil3d S;

  S.nx=nx; S.ny=ny; S.nz=nz;
  S.value_c = 9;
  S.value_n = 1;
  S.value_e = 2;
  S.value_s = 1;
  S.value_w = 2;
  S.value_b = 1;
  S.value_t = 1;

  for (int i=0; i<n; i++)
  {
    e[i]=1.0;
    if (i>0) e[i-1]=0.0;
    apply_stencil3d(&S, e, A+i*n);
  }

  int wrong_rows=0;
  for (int i=0; i<n; i++){
      // if (A[i*n+j]!=A[j*n+i]) wrong_entries++;
      double row_sum = 0.0;
      for (int j=0; j<n; j++){
        if (j != i) row_sum += A[i*n+j];
      }
      if (A[i*n+i] < row_sum) wrong_rows++;
    }
  EXPECT_EQ(0, wrong_rows);

  if (wrong_rows)
  {
    std::cout << "Your matrix (computed on a 3x4x5 grid by apply_stencil(I)) is ..."<<std::endl;
    for (int j=0; j<n; j++)
    {
      for (int i=0; i<n; i++)
      {
        std::cout << A[i*n+j] << " ";
      }
      std::cout << std::endl;
    }
  }
  delete [] e;
  delete [] A;
}
