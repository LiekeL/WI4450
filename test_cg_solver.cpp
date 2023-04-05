#include "gtest_mpi.hpp"

#include "operations.hpp"
#include "cg_solver.hpp"

#include <iostream>
#include <cmath>

TEST(cg_solver, poisson){
    const int nx=10, ny=8, nz=6;    // define dimensions
    const int n=nx*ny*nz;           // define length vector
    double* rhs=new double[n];      // define right hand side vector
    double* u_known=new double[n];  // define solution vector
    double* u=new double[n];        // define solution vector

    for (int i=0;i<n;i++){      
        u_known[i] = i%3;   // solution is known
        u[i] = 1.0;         // starting vector for the conjugate gradient method
    }

    stencil3d S;            // create a symmetric operation

    S.nx=nx; S.ny=ny; S.nz=nz;
    S.value_c = 6;
    S.value_n = -1;
    S.value_e = -1;
    S.value_s = -1;
    S.value_w = -1;
    S.value_b = -1;
    S.value_t = -1;

    apply_stencil3d(&S,u_known,rhs);    // compute the right hand side with the known u_known

    double cg_tol = 10e-20;
    double tol = 10e-6;
    int maxIter = 1000;
    int numIter = 0;
    double resNorm = 10e6;

    cg_solver(&S, n, u, rhs, cg_tol, maxIter, &resNorm, &numIter, 0);

    double err=0.0;

    for (int i=0; i<n; i++) err = std::max(err, std::abs(u[i]-u_known[i]));

    EXPECT_NEAR(1.0+err, 1.0, tol);
    
    delete [] rhs;
    delete [] u_known;
    delete [] u;
}

TEST(cg_solver, identity){
    const int nx=3, ny=8, nz=6;    // define dimensions
    const int n=nx*ny*nz;           // define length vector
    double* u_known=new double[n];  // define solution vector
    double* u=new double[n];        // define solution vector

    for (int i=0;i<n;i++){      
        u_known[i] = i%3;   // solution is known
        u[i] = 1.0;         // starting vector for the conjugate gradient method
    }

    stencil3d S;            // create the identity stencil

    S.nx=nx; S.ny=ny; S.nz=nz;
    S.value_c = 1;
    S.value_n = 0;
    S.value_e = 0;
    S.value_s = 0;
    S.value_w = 0;
    S.value_b = 0;
    S.value_t = 0;

    double cg_tol = 10e-20;
    double tol = 10e-6;
    int maxIter = 1;
    int numIter = 0;
    double resNorm = 10e6;

    cg_solver(&S, n, u, u_known, cg_tol, maxIter, &resNorm, &numIter, 0); // solve I*u=u_known

    double err=0.0;

    for (int i=0; i<n; i++) err = std::max(err, std::abs(u[i]-u_known[i]));

    EXPECT_NEAR(1.0+err, 1.0, tol);
    
    delete [] u_known;
    delete [] u;
}