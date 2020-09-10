/*************************************************************************
 * IPM_MC_PK (primal-dual predictor-corrector interior-point method       *
 * for basic SDP relaxation for Max-Cut) solves:                         *
 *                                                                       *
 * max tr(LX), subject to diag(X) = e, X psd                             *   
 * min e'y, subject to Z = Diag(y) - L psd, y unconstrained              *
 *                                                                       *
 * input:  L    ... Laplacian matrix of the graph (1/4 of L for max-cut) *
 *         n    ... size of the problem                                  *
 *         print... print level                                          *
 * output: phi  ... optimal value of SDP (value of the dual problem)     *
 *         X    ... optimal primal matrix                                *
 *         y    ... optimal dual vector                                  *   
 *         Z    ... optimal dual matrix                                  *
 *************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "admm.h"

/* NOTE: C uses row-major, but blas and lapack routines (written in Fortran)
 * use column-major --> be careful when multiplying non-symmetric matrices) */

void ipm_mc_pk(double *L, int n, double *X, double *y, double *Z, double *phi, int print) {

    /* variables for blas and lapack routines */
    int inc = 1;
    char up = 'U';              // for lapack take upper triangular part of matrix 
    char side = 'L';            // in matrix product AB, the left matrix is symmetric 
    int info;                   // test whether lapack function succeded
    double alpha, beta;         // scalars in linear combination (lapack)

    /* other variables */
    int i, j, k;                // loop iter
    double *p, *p2, *p3;        // pointers in loops
    double psi;                 // value of primal problem
    double mu;                  // ZX = mu * I (parametrized optimality condition)
    double alpha_p, alpha_d;    // step lengths

    double *b;                  // vector of ones
    double *dX, *dy, *dZ;       
    double *Zi;                 // inv(Z)
    double *M;                  // M * dy = rhs
    double *dy1, *dX1;
    double *tmp, *tmp2;         // need non-symm matrices when computing for instance Zi*diag(dy)*X       

    /*************************************************
     * initial positive definite matrices X, Z and y *
     * primal,dual cost, gap                         *
     *************************************************/
    int nn = n*n;

    /* set y to zero vector */
    for (int i = 0; i < n; ++i)
        y[i] = 0.0;

    /* y = 1.1 + sum(abs(L))' */
    for (i = 1, p = L, p2 = y; i <= nn; ++i, ++p) {
        *p2 += fabs(*p);
        if ((i % n) == 0) {
            *p2 += 1.1;
            ++p2;
        }
    }

    /* vector b of all ones */
    alloc_vector(b, n, double);
    for (i = 0; i < n; ++i)
        b[i] = 1.0;
    
    /* X = eye(n) = Diag(b) */
    Diag(X, b, n);

    /* Z = Diag(y) - L */
    Diag(Z, y, n);                           /* Z = Diag(y) */
    alpha = -1.0;
    daxpy_(&nn,&alpha,L,&inc,Z,&inc);        /* Z = Z - L */


    /* phi = ones(n,1)'*y */                 /* initial dual value */
    *phi = ddot_(&n,b,&inc,y,&inc);

    /* psi = L(:)'*X(:) */                   /* initial primal value */
    psi = ddot_(&nn,L,&inc,X,&inc);

    /* mu = Z(:)'*X(:) / (2*n); */           /* initial complementarity */ 
    mu = ddot_(&nn,Z,&inc,X,&inc) / (2.0 * n);


    /* print output */
    if (print) {
        puts("iter     log10(gap)     primal        dual");
        puts("*******************************************"); 
    }

    /*  allocate space */
    alloc_matrix(dX, n, double);
    alloc_vector(dy, n, double);
    alloc_matrix(dZ, n, double);

    alloc_vector(dy1, n, double);
    alloc_matrix(dX1, n, double);
    alloc_matrix(tmp, n, double);
    alloc_matrix(tmp2, n, double);

    alloc_matrix(Zi, n, double);    
    alloc_matrix(M, n, double);   



    /*************
     * main loop *
     *************/
    for (i = 1; *phi - psi > 1e-4; ++i) {/* while duality gap too large */
        
        /******** compute inverse of Z ********/
        dcopy_(&nn,Z,&inc,Zi,&inc);         /* copy Z to Zi */   

        dpotrf_(&up,&n,Zi,&n,&info);        /* computes Cholesky factorization */
        if (info != 0) {
            fprintf(stderr, "%s: Problem with Cholesky factorization \
                (line: %d).\n", __func__, __LINE__);
            exit(EXIT_FAILURE);
        }

        dpotri_(&up,&n,Zi,&n,&info);        /* computes Zi = inv(Z) */
        if (info != 0) {
            fprintf(stderr, "%s: Problem with computation of inverse matrix \
                (line: %d).\n", __func__, __LINE__);
            exit(EXIT_FAILURE);
        }

        /* NOTE: only upper triangular part of Zi is ok and can be used! */
        /* copy strictly upper triangular part of Zi to strictly lower triangular part */
        for (j = 0; j < n; ++j) {       
            for (k = 0; k < j; ++k)
                Zi[j+n*k] = Zi[k+n*j];
        }


        /********   predictor step (mu = 0) solves:       ********
         ********   Z * X + Diag(dy1) * X + Z * dX1 = 0   ********/

        /* M = Zi .* X */
        for (j = 0, p = M; j < nn; ++j, ++p)
            *p = Zi[j] * X[j];

        /* copy -b to dy1 */
        alpha = -1.0;
        dcopy_(&n,b,&inc,dy1,&inc);  // copy
        dscal_(&n,&alpha,dy1,&inc);  // scale     

        /* solve the system: dy1 = (Zi .* X) \ (-e) */
        dposv_(&up, &n, &inc, M, &n, dy1, &n, &info);
        /* NOTE: M and dy1 are changed on exit! */

        if (info != 0) {
            fprintf(stderr, "%s: predictor step: problem in solving linear system \
                (line: %d).\n", __func__, __LINE__);
            exit(EXIT_FAILURE);
        }

        /* dX1 = -Zi*diag(dy1)*X - X */
        // NOTE: be carefull to multiply in column-major ordering!

        /* 1. step: tmp = -diag(dy1)*X */
        /* = multiply j-th row of X by -dy1[j]        (FORTRAN) */
        /* = multiply j-th column of X by -dy[j]      (C)       */
        p3 = tmp; 
        for (j = 0, p = X, p2 = dy1; j < nn; ++j, ++p, ++p3) {
            *p3 = -(p2[j%n]) * (*p);
        }

        /* 2. step: Zi * tmp */
        alpha = 1.0;
        beta = 0.0;
        dsymm_(&side,&up,&n,&n,&alpha,Zi,&n,tmp,&n,&beta,dX1,&n);
     
        /* dX1 = dX1 - X */
        alpha = -1.0;
        daxpy_(&nn,&alpha,X,&inc,dX1,&inc);

        /* symmetrize: dX1 = (dX1 + dX1')/2 */
        for (j = 0; j < n; ++j) {       
            for (k = 0; k <= j; ++k)    
                dX1[k+n*j] = dX1[j+n*k] = 0.5 * (dX1[j+n*k] + dX1[k+n*j]);
        }

  
        /*************** corrector step solves:  ******************/  
        /******** diag(dy2)*X + Z*dX2 - mu*I + diag(dy1)*dX1 = 0 **/

        /* dy2 = M \ (mu*diag(Zi) - (Zi .* dX1)*dy1) */

        /* dy = diag(Zi) */
        diag(Zi,dy,n);     

        /* tmp = Zi .* dX1 */
        for (j = 0, p = tmp; j < nn; ++j, ++p)
            *p = Zi[j] * dX1[j];

        /* dy = -tmp*dy1 + mu*dy */
        alpha = -1.0;
        dsymv_(&up,&n,&alpha,tmp,&n,dy1,&inc,&mu,dy,&inc);

        // NOTE: M was changed during dposv
        /* M = Zi .* X */
        for (j = 0, p = M; j < nn; ++j, ++p)
            *p = Zi[j] * X[j];

        /* dy2 = M \ dy */
        /* dy2 is saved in dy! */
        dposv_(&up, &n, &inc, M, &n, dy, &n, &info);
        if (info != 0) {
            fprintf(stderr, "%s: corrector step: problem in solving linear system \
                (line: %d).\n", __func__, __LINE__);
            exit(EXIT_FAILURE);
        }

        /* dX2 = mu*Zi - Zi*( diag(dy2) * X + diag(dy1) * dX1) */

        /* 1. step: tmp = -diag(dy2)*X */
        /* = multiply j-th row of X by -dy2[j]       (FORTRAN) */
        /* = multiply j-th column of X by -dy2[j]    (C)       */    
        /* NOTE: dy2 = dy */
        p3 = tmp; 
        for (j = 0, p = X, p2 = dy; j < nn; ++j, ++p, ++p3) {
            *p3 = -(p2[j%n]) * (*p);
        }

        /* 2. step: tmp2 = -diag(dy1)*dX1 */
        p3 = tmp2; 
        for (j = 0, p = dX1, p2 = dy1; j < nn; ++j, ++p, ++p3) {
            *p3 = -(p2[j%n]) * (*p);
        }
        
        /* tmp = tmp + tmp2 */
        alpha = 1.0;
        daxpy_(&nn,&alpha,tmp2,&inc,tmp,&inc);

        /* dX2 = Zi * tmp*/
        /* NOTE: dX2 is stored in dX */
        alpha = 1.0;
        beta = 0.0;
        dsymm_(&side,&up,&n,&n,&alpha,Zi,&n,tmp,&n,&beta,dX,&n);
        
        /* dX = dX2 = mu*Zi + dX2 */
        daxpy_(&nn,&mu,Zi,&inc,dX,&inc);

      
        /**** final steps ****/
        alpha = 1.0;
        daxpy_(&n,&alpha,dy1,&inc,dy,&inc);     /* dy = dy1 + dy */
        daxpy_(&nn,&alpha,dX1,&inc,dX,&inc);    /* dX = dX1 + dX */
        
        /* symmetrize: dX = (dX + dX')/2 */
        for (j = 0; j < n; ++j) {       
            for (k = 0; k <= j; ++k)    
                dX[k+n*j] = dX[j+n*k] = 0.5 * (dX[j+n*k] + dX[k+n*j]);
        }

        /* dZ = Diag(dy) */
        Diag(dZ,dy,n);
                

        /*********** find step lengths alpha_p and alpha_d ***********/
        
        /* line search on primal: X  = X + alpha_p * dX  psd matrix */
        alpha_p = 1.0;
        info = 1;
        while (info != 0) {
            dcopy_(&nn,X,&inc,tmp,&inc);            /* tmp = X */
            daxpy_(&nn,&alpha_p,dX,&inc,tmp,&inc);  /* tmp = alpha_p * dX + tmp */
            dpotrf_(&up,&n,tmp,&n,&info);

            if (info != 0)
                alpha_p *= 0.8;
        }
        
        if (alpha_p < 1.0)                          /* stay away from boundary */
            alpha_p *= 0.95;


        /* line search on dual */
        alpha_d = 1.0;
        info = 1;
        while (info != 0) {
            dcopy_(&nn,Z,&inc,tmp,&inc);            /* tmp = Z */
            daxpy_(&nn,&alpha_d,dZ,&inc,tmp,&inc);  /* tmp = alpha_d * dZ + tmp */
            dpotrf_(&up,&n,tmp,&n,&info);

            if (info != 0)
                alpha_d *= 0.8;
        }
        
        if (alpha_d < 1.0)                          /* stay away from boundary */
            alpha_d *= 0.95;


        /******** update ********/
        daxpy_(&nn,&alpha_p,dX,&inc,X,&inc);        /* X = alpha_p * dX + X */    
        daxpy_(&n,&alpha_d,dy,&inc,y,&inc);         /* y = alpha_d * dy + y */
        daxpy_(&nn,&alpha_d,dZ,&inc,Z,&inc);        /* Z = alpha_d * dZ + Z */

        /* mu = Z(:)'*X(:) / (2*n); */            
        mu = ddot_(&nn,Z,&inc,X,&inc) / (2.0 * n);   

        /* speed up for long steps */
        if (alpha_p + alpha_d > 1.6)
            mu *= 0.5;
        if (alpha_p + alpha_d > 1.9)
            mu *= 0.2;

        /***** objective values *****/

        /* phi = ones(n,1)'*y */                 
        *phi = ddot_(&n,b,&inc,y,&inc);

        /* psi = L(:)'*X(:) */                   
        psi = ddot_(&nn,L,&inc,X,&inc);


        /* print output */
        if (print)
            printf("%3d %11.2f %14.5f %14.5f \n",i,log10(*phi-psi),psi,*phi);

    } // end of main loop

    if (print)
        puts("*******************************************");

    /************* free memory *************/
    free(b);    
    free(dX);
    free(dy);
    free(dZ);
    free(Zi);
    free(M);
    free(dy1);
    free(dX1);
    free(tmp);
    free(tmp2);
}
