#include <math.h>

#include "admm.h"

extern Triangle_Inequality *Cuts;           // vector of triangle inequality constraints
extern Pentagonal_Inequality *Pent_Cuts;    // vector of pentagonal inequality constraints
extern Heptagonal_Inequality *Hepta_Cuts;   // vector of heptagonal inequality constraints

/* primal and dual matrices */
extern double *X;  
extern double *Q;

/* for eigenvalues and eigenvectors */
extern int M; 
extern double *W; 
extern double *Z; 

extern int *ISUPPZ; 
extern double *WORK; 
extern int *IWORK; 
extern int sizeWORK, sizeIWORK;

// NOTE: matrices are symmetric --> only need to use upper triangle
// NOTE: be careful when calling fortran routines (blas/lapack)
//       because it stores data as column major --> use L for lower triangular part!


/***************** diag *********************/
/* 
 * operator diag: save diagonal of matrix X in vector y
 */
void diag(const double *X, double *y, int n) {

    for (int i = 0; i < n; ++i) {
        y[i] = X[i + i * n];
    }
}


/***************** Diag *********************/
/* 
 * operator Diag: make diagonal matrix X from vector y
 */
void Diag(double *X, const double *y, int n) {

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j)
                X[i + i * n] = y[i];
            else
                X[j + i * n] = 0.0;
        }
    }
}


/***************** op_B *********************/
/*
 * computes y = y + B(X), where operator B
 * corresponds to triangle inequalities
 */
void op_B(const Problem *P, double *y, double *X) {

    int N = P->n;
    int type, ii, jj, kk, ll, mm, nn, oo;

    /* copy upper triangle to lower */
    for (int j = 0; j < N; ++j) {       
        for (int k = 0; k < j; ++k)    
            X[k + N * j] = X[j + N * k];
    }

    
    for (int ineq = 0; ineq < P->NIneq; ++ineq) {

        type = Cuts[ineq].type;
        ii   = Cuts[ineq].i;
        jj   = Cuts[ineq].j;
        kk   = Cuts[ineq].k;

        switch (type) {
            // triangle inequalities are always with i,j,k for which it holds
            // k < j < i !!
            // --> use upper triangle of matrix X
            case 1:
                y[ineq] += (-X[ii + jj * N] - X[ii + kk * N] - X[jj + kk * N]);
                break;
            case 2:
                y[ineq] += (-X[ii + jj * N] + X[ii + kk * N] + X[jj + kk * N]);
                break;
            case 3:
                y[ineq] += (X[ii + jj * N] - X[ii + kk * N] + X[jj + kk * N]);
                break;
            default: //type == 4
                y[ineq] += (X[ii + jj * N] + X[ii + kk * N] - X[jj + kk * N]);
        }
    }    

    /* pentagonal inequalities */
    for (int ineq = 0; ineq < P->NPentIneq; ++ineq) {

        type = Pent_Cuts[ineq].type;
        ii   = Pent_Cuts[ineq].permutation[0];
        jj   = Pent_Cuts[ineq].permutation[1];
        kk   = Pent_Cuts[ineq].permutation[2];
        ll   = Pent_Cuts[ineq].permutation[3];
        mm   = Pent_Cuts[ineq].permutation[4];

        switch (type) {

            case 1: // H1 = ee^T, where e is vector of all ones 
                y[P->NIneq + ineq] += 0.5 * (-X[ii + jj * N] - X[ii + kk * N] - X[ii + ll * N] - X[ii + mm * N] - X[jj + kk * N] - X[jj + ll * N] - X[jj + mm * N] - X[kk + ll * N] - X[kk + mm * N] - X[ll + mm * N] );
                break;
            case 2: // e[0] = -1; H2 = ee^T
                y[P->NIneq + ineq] += 0.5 * (X[ii + jj * N] + X[ii + kk * N] + X[ii + ll * N] + X[ii + mm * N] - X[jj + kk * N] - X[jj + ll * N] - X[jj + mm * N] - X[kk + ll * N] - X[kk + mm * N] - X[ll + mm * N] );
                break;
            default: // type = 3; // e[0] = -1, e[1] = -1; H3 = ee^T
                y[P->NIneq + ineq] += 0.5 * (-X[ii + jj * N] + X[ii + kk * N] + X[ii + ll * N] + X[ii + mm * N] + X[jj + kk * N] + X[jj + ll * N] + X[jj + mm * N] - X[kk + ll * N] - X[kk + mm * N] - X[ll + mm * N] );                    
        }

    }

    /* heptagonal inequalities */
    for (int ineq = 0; ineq < P->NHeptaIneq; ++ineq) {

        type = Hepta_Cuts[ineq].type;
        ii   = Hepta_Cuts[ineq].permutation[0];
        jj   = Hepta_Cuts[ineq].permutation[1];
        kk   = Hepta_Cuts[ineq].permutation[2];
        ll   = Hepta_Cuts[ineq].permutation[3];
        mm   = Hepta_Cuts[ineq].permutation[4]; 
        nn   = Hepta_Cuts[ineq].permutation[5];
        oo   = Hepta_Cuts[ineq].permutation[6];   

        switch (type) {

            case 1: // H1 = ee^T, where e is vector of all ones 
                y[P->NIneq + P->NPentIneq + ineq] += 1.0/3.0 * (-X[ii + jj * N] - X[ii + kk * N] - X[ii + ll * N] - X[ii + mm * N] - X[ii + nn * N] - X[ii + oo * N] - X[jj + kk * N] - X[jj + ll * N] - X[jj + mm * N] - X [jj + nn * N] - X[jj + oo * N] - X[kk + ll * N] - X[kk + mm * N] - X[kk + nn * N] - X[kk + oo * N] - X[ll + mm * N] - X[ll + nn * N] - X[ll + oo * N] - X[mm + nn * N] - X[mm + oo * N] - X[nn + oo * N]);
                break;
            case 2: // e[0] = -1; H2 = ee^T
                y[P->NIneq + P->NPentIneq + ineq] += 1.0/3.0 * (X[ii + jj * N] + X[ii + kk * N] + X[ii + ll * N] + X[ii + mm * N] + X[ii + nn * N] + X[ii + oo * N] - X[jj + kk * N] - X[jj + ll * N] - X[jj + mm * N] - X[jj + nn * N] - X[jj + oo * N] - X[kk + ll * N] - X[kk + mm * N] - X[kk + nn * N] - X[kk + oo * N] - X[ll + mm * N] - X[ll + nn * N] - X[ll + oo * N] - X[mm + nn * N] - X[mm + oo * N] - X[nn + oo * N]);
                break;    
            case 3: // H1 = ee^T, where e is vector of all ones 
                y[P->NIneq + P->NPentIneq + ineq] += 1.0/3.0 * (-X[ii + jj * N] + X[ii + kk * N] + X[ii + ll * N] + X[ii + mm * N] + X[ii + nn * N] + X[ii + oo * N] + X[jj + kk * N] + X[jj + ll * N] + X[jj + mm * N] + X[jj + nn * N] + X[jj + oo * N] - X[kk + ll * N] - X[kk + mm * N] - X[kk + nn * N] - X[kk + oo * N] - X[ll + mm * N] - X[ll + nn * N] - X[ll + oo * N] - X[mm + nn * N] - X[mm + oo * N] - X[nn + oo * N]);
                break;   
            default: // e[0] = -1, e[1] = -1, e[2] = -1; H4 = ee^T
                y[P->NIneq + P->NPentIneq + ineq] += 1.0/3.0 * (-X[ii + jj * N] - X[ii + kk * N] + X[ii + ll * N] + X[ii + mm * N] + X[ii + nn * N] + X[ii + oo * N] - X[jj + kk * N] + X[jj + ll * N] + X[jj + mm * N] + X [jj + nn * N] + X[jj + oo * N] + X[kk + ll * N] + X[kk + mm * N] + X[kk + nn * N] + X[kk + oo * N] - X[ll + mm * N] - X[ll + nn * N] - X[ll + oo * N] - X[mm + nn * N] - X[mm + oo * N] - X[nn + oo * N]);   
        }  
    }
}


/***************** op_Diag *********************/
/*
 * computes X = X - Diag(y)
 */
void op_Diag(const Problem *P, double *X, const double *yy) {

    int N = P->n;

    for (int i = 0; i < N; ++i) {
        X[i + i * N] -= yy[i];
    }
}

/***************** op_Bt *********************/
/*
 * computes X = X - Bt(t), where operator B
 * corresponds to triangle inequalities
 */
void op_Bt(const Problem *P, double *X, const double *tt) {

    int N = P->n;
    int type, ii, jj, kk, ll, mm, nn, oo;
    double temp;

    /***** triangle inequalitiess *****/ 
    for (int ineq = 0; ineq < P->NIneq; ++ineq) {

        type = Cuts[ineq].type;
        ii   = Cuts[ineq].i;
        jj   = Cuts[ineq].j;
        kk   = Cuts[ineq].k;

        temp = 0.5 * tt[ineq];

        switch (type) {
            
            case 1:
                X[ii + jj * N] += temp; 
                X[ii + kk * N] += temp; 
                X[jj + kk * N] += temp; 
                X[jj + ii * N] += temp;
                X[kk + ii * N] += temp;
                X[kk + jj * N] += temp;
                break;
            case 2:
                X[ii + jj * N] += temp; 
                X[ii + kk * N] -= temp; 
                X[jj + kk * N] -= temp; 
                X[jj + ii * N] += temp;
                X[kk + ii * N] -= temp;
                X[kk + jj * N] -= temp;
                break;
            case 3:
                X[ii + jj * N] -= temp; 
                X[ii + kk * N] += temp; 
                X[jj + kk * N] -= temp; 
                X[jj + ii * N] -= temp;
                X[kk + ii * N] += temp;
                X[kk + jj * N] -= temp;
                break;
            default: //type == 4
                X[ii + jj * N] -= temp; 
                X[ii + kk * N] -= temp; 
                X[jj + kk * N] += temp; 
                X[jj + ii * N] -= temp;
                X[kk + ii * N] -= temp;
                X[kk + jj * N] += temp;
        }
    }

    /***** pentagonal inequalities *****/
    for (int ineq = 0; ineq < P->NPentIneq; ++ineq) {

        type = Pent_Cuts[ineq].type;
        ii   = Pent_Cuts[ineq].permutation[0];
        jj   = Pent_Cuts[ineq].permutation[1];
        kk   = Pent_Cuts[ineq].permutation[2];
        ll   = Pent_Cuts[ineq].permutation[3];
        mm   = Pent_Cuts[ineq].permutation[4];

        temp = 0.25 * tt[P->NIneq + ineq]; // 0.5 due to symmetry and 0.5 due to pentagonal inequality -0.5*(...) <= 1 

        switch (type) {
            case 1:
                X[ii + jj * N] += temp; 
                X[ii + kk * N] += temp; 
                X[ii + ll * N] += temp;
                X[ii + mm * N] += temp;
                X[jj + kk * N] += temp;
                X[jj + ll * N] += temp;
                X[jj + mm * N] += temp;
                X[kk + ll * N] += temp;
                X[kk + mm * N] += temp;
                X[ll + mm * N] += temp;

                X[jj + ii * N] += temp; 
                X[kk + ii * N] += temp; 
                X[ll + ii * N] += temp;
                X[mm + ii * N] += temp;
                X[kk + jj * N] += temp;
                X[ll + jj * N] += temp;
                X[mm + jj * N] += temp;
                X[ll + kk * N] += temp;
                X[mm + kk * N] += temp;
                X[mm + ll * N] += temp;
                break;

            case 2:
                X[ii + jj * N] -= temp; 
                X[ii + kk * N] -= temp; 
                X[ii + ll * N] -= temp;
                X[ii + mm * N] -= temp;
                X[jj + kk * N] += temp;
                X[jj + ll * N] += temp;
                X[jj + mm * N] += temp;
                X[kk + ll * N] += temp;
                X[kk + mm * N] += temp;
                X[ll + mm * N] += temp;

                X[jj + ii * N] -= temp; 
                X[kk + ii * N] -= temp; 
                X[ll + ii * N] -= temp;
                X[mm + ii * N] -= temp;
                X[kk + jj * N] += temp;
                X[ll + jj * N] += temp;
                X[mm + jj * N] += temp;
                X[ll + kk * N] += temp;
                X[mm + kk * N] += temp;
                X[mm + ll * N] += temp;
                break;    

            case 3:
                X[ii + jj * N] += temp; 
                X[ii + kk * N] -= temp; 
                X[ii + ll * N] -= temp;
                X[ii + mm * N] -= temp;
                X[jj + kk * N] -= temp;
                X[jj + ll * N] -= temp;
                X[jj + mm * N] -= temp;
                X[kk + ll * N] += temp;
                X[kk + mm * N] += temp;
                X[ll + mm * N] += temp;

                X[jj + ii * N] += temp; 
                X[kk + ii * N] -= temp; 
                X[ll + ii * N] -= temp;
                X[mm + ii * N] -= temp;
                X[kk + jj * N] -= temp;
                X[ll + jj * N] -= temp;
                X[mm + jj * N] -= temp;
                X[ll + kk * N] += temp;
                X[mm + kk * N] += temp;
                X[mm + ll * N] += temp;
                break;      
        }

    }  

    /***** heptagonal inequalities *****/
    for (int ineq = 0; ineq < P->NHeptaIneq; ++ineq) {

        type = Hepta_Cuts[ineq].type;
        ii   = Hepta_Cuts[ineq].permutation[0];
        jj   = Hepta_Cuts[ineq].permutation[1];
        kk   = Hepta_Cuts[ineq].permutation[2];
        ll   = Hepta_Cuts[ineq].permutation[3];
        mm   = Hepta_Cuts[ineq].permutation[4]; 
        nn   = Hepta_Cuts[ineq].permutation[5];
        oo   = Hepta_Cuts[ineq].permutation[6];

        temp = 1.0/6.0 * tt[P->NIneq + P->NPentIneq + ineq]; // 0.5 due to symmetry and 1/3 due to heptagonal inequality -1/3*(...) <= 1 

        switch (type) {

            case 1:
                X[ii + jj * N] += temp; 
                X[ii + kk * N] += temp; 
                X[ii + ll * N] += temp;
                X[ii + mm * N] += temp;
                X[ii + nn * N] += temp;
                X[ii + oo * N] += temp;
                X[jj + kk * N] += temp;
                X[jj + ll * N] += temp;
                X[jj + mm * N] += temp;
                X[jj + nn * N] += temp;
                X[jj + oo * N] += temp;
                X[kk + ll * N] += temp;
                X[kk + mm * N] += temp;
                X[kk + nn * N] += temp;
                X[kk + oo * N] += temp;
                X[ll + mm * N] += temp;
                X[ll + nn * N] += temp;
                X[ll + oo * N] += temp;
                X[mm + nn * N] += temp;
                X[mm + oo * N] += temp;
                X[nn + oo * N] += temp;

                X[jj + ii * N] += temp; 
                X[kk + ii * N] += temp; 
                X[ll + ii * N] += temp;
                X[mm + ii * N] += temp;
                X[nn + ii * N] += temp;
                X[oo + ii * N] += temp;
                X[kk + jj * N] += temp;
                X[ll + jj * N] += temp;
                X[mm + jj * N] += temp;
                X[nn + jj * N] += temp;
                X[oo + jj * N] += temp;
                X[ll + kk * N] += temp;
                X[mm + kk * N] += temp;
                X[nn + kk * N] += temp;
                X[oo + kk * N] += temp;
                X[mm + ll * N] += temp;
                X[nn + ll * N] += temp;
                X[oo + ll * N] += temp;
                X[nn + mm * N] += temp;
                X[oo + mm * N] += temp;
                X[oo + nn * N] += temp;
                break;

            case 2:
                X[ii + jj * N] -= temp; 
                X[ii + kk * N] -= temp; 
                X[ii + ll * N] -= temp;
                X[ii + mm * N] -= temp;
                X[ii + nn * N] -= temp;
                X[ii + oo * N] -= temp;
                X[jj + kk * N] += temp;
                X[jj + ll * N] += temp;
                X[jj + mm * N] += temp;
                X[jj + nn * N] += temp;
                X[jj + oo * N] += temp;
                X[kk + ll * N] += temp;
                X[kk + mm * N] += temp;
                X[kk + nn * N] += temp;
                X[kk + oo * N] += temp;
                X[ll + mm * N] += temp;
                X[ll + nn * N] += temp;
                X[ll + oo * N] += temp;
                X[mm + nn * N] += temp;
                X[mm + oo * N] += temp;
                X[nn + oo * N] += temp;

                X[jj + ii * N] -= temp; 
                X[kk + ii * N] -= temp; 
                X[ll + ii * N] -= temp;
                X[mm + ii * N] -= temp;
                X[nn + ii * N] -= temp;
                X[oo + ii * N] -= temp;
                X[kk + jj * N] += temp;
                X[ll + jj * N] += temp;
                X[mm + jj * N] += temp;
                X[nn + jj * N] += temp;
                X[oo + jj * N] += temp;
                X[ll + kk * N] += temp;
                X[mm + kk * N] += temp;
                X[nn + kk * N] += temp;
                X[oo + kk * N] += temp;
                X[mm + ll * N] += temp;
                X[nn + ll * N] += temp;
                X[oo + ll * N] += temp;
                X[nn + mm * N] += temp;
                X[oo + mm * N] += temp;
                X[oo + nn * N] += temp;
                break;  

            case 3:
                X[ii + jj * N] += temp; 
                X[ii + kk * N] -= temp; 
                X[ii + ll * N] -= temp;
                X[ii + mm * N] -= temp;
                X[ii + nn * N] -= temp;
                X[ii + oo * N] -= temp;
                X[jj + kk * N] -= temp;
                X[jj + ll * N] -= temp;
                X[jj + mm * N] -= temp;
                X[jj + nn * N] -= temp;
                X[jj + oo * N] -= temp;
                X[kk + ll * N] += temp;
                X[kk + mm * N] += temp;
                X[kk + nn * N] += temp;
                X[kk + oo * N] += temp;
                X[ll + mm * N] += temp;
                X[ll + nn * N] += temp;
                X[ll + oo * N] += temp;
                X[mm + nn * N] += temp;
                X[mm + oo * N] += temp;
                X[nn + oo * N] += temp;

                X[jj + ii * N] += temp; 
                X[kk + ii * N] -= temp; 
                X[ll + ii * N] -= temp;
                X[mm + ii * N] -= temp;
                X[nn + ii * N] -= temp;
                X[oo + ii * N] -= temp;
                X[kk + jj * N] -= temp;
                X[ll + jj * N] -= temp;
                X[mm + jj * N] -= temp;
                X[nn + jj * N] -= temp;
                X[oo + jj * N] -= temp;
                X[ll + kk * N] += temp;
                X[mm + kk * N] += temp;
                X[nn + kk * N] += temp;
                X[oo + kk * N] += temp;
                X[mm + ll * N] += temp;
                X[nn + ll * N] += temp;
                X[oo + ll * N] += temp;
                X[nn + mm * N] += temp;
                X[oo + mm * N] += temp;
                X[oo + nn * N] += temp;
                break;      

            case 4:
                X[ii + jj * N] += temp; 
                X[ii + kk * N] += temp; 
                X[ii + ll * N] -= temp;
                X[ii + mm * N] -= temp;
                X[ii + nn * N] -= temp;
                X[ii + oo * N] -= temp;
                X[jj + kk * N] += temp;
                X[jj + ll * N] -= temp;
                X[jj + mm * N] -= temp;
                X[jj + nn * N] -= temp;
                X[jj + oo * N] -= temp;
                X[kk + ll * N] -= temp;
                X[kk + mm * N] -= temp;
                X[kk + nn * N] -= temp;
                X[kk + oo * N] -= temp;
                X[ll + mm * N] += temp;
                X[ll + nn * N] += temp;
                X[ll + oo * N] += temp;
                X[mm + nn * N] += temp;
                X[mm + oo * N] += temp;
                X[nn + oo * N] += temp;

                X[jj + ii * N] += temp; 
                X[kk + ii * N] += temp; 
                X[ll + ii * N] -= temp;
                X[mm + ii * N] -= temp;
                X[nn + ii * N] -= temp;
                X[oo + ii * N] -= temp;
                X[kk + jj * N] += temp;
                X[ll + jj * N] -= temp;
                X[mm + jj * N] -= temp;
                X[nn + jj * N] -= temp;
                X[oo + jj * N] -= temp;
                X[ll + kk * N] -= temp;
                X[mm + kk * N] -= temp;
                X[nn + kk * N] -= temp;
                X[oo + kk * N] -= temp;
                X[mm + ll * N] += temp;
                X[nn + ll * N] += temp;
                X[oo + ll * N] += temp;
                X[nn + mm * N] += temp;
                X[oo + mm * N] += temp;
                X[oo + nn * N] += temp;
                break;     
        }        
    } 

}


/********** project matrix onto PSD cone **********/
/* input: X = L + X/sigma - Diag(y) - Bt(y)
 *
 * computes X_+ and Q = X - X_+
 */
void project_SDP(const Problem *P) {

    int N = P->n;
    int nn = N * N;

    /* Copy -X into Q (used later for Q = X - X_+) */
    int INCX = 1;
    int INCQ = 1;
    double ALPHA = -1.0;

    dcopy_(&nn, X, &INCX, Q, &INCQ);
    dscal_(&nn, &ALPHA, Q, &INCQ);

    // upper and lower bound on positive eigenvalues
    double VU;        
    double VL = 1e-8;
    
    char UPLO = 'L';
    int LDX = N;

    /* compute operator norm for X to get an upper bound for eigenvalues
     * |X|_2 <= sqrt{ |X|_1 |X|_inf }  (Holder's inequality)
     *        = |X|_1 = |X|_inf        (since X is symmetric)
     *
     * Frobenius norm is also an upper bound on the spectral radius of X:
     *      |X|_2 <= |X|_F
     * 
     * --> use the one that is smaller
     */
    char NORM = 'I';
    double norm_inf = dlansy_(&NORM, &UPLO, &N, X, &LDX, WORK);

    NORM = 'F';
    double norm_fro = dlansy_(&NORM, &UPLO, &N, X, &LDX, WORK);

    // VU = min(norm_inf, norm_fro)
    VU = (norm_inf < norm_fro) ? norm_inf : norm_fro;

    // Ensure that VL <= VU.
    if (VU < VL) {
        VU = 2.0 * VL;
    }

    /* Compute the positive eigenvalues and associated eigenvectors of X.
     *
     * The M columns of Z will contain the orthonormal eigenvectors of the 
     * matrix X corresponding to the positive eigenvalues, the i-th column 
     * of Z holding the eigenvector associated with W[i].
     */
    char JOBZ = 'V';
    char RANGE = 'V';
    int IL = 0;
    int IU = 0;
    double ABSTOL = 1e-8;
    int LDZ = N;
    int LWORK = sizeWORK;
    int LIWORK = sizeIWORK;
    int INFO;
    dsyevr_(&JOBZ, &RANGE, &UPLO, &N, X, &LDX, &VL, &VU, &IL, &IU, &ABSTOL, 
            &M, W, Z, &LDZ, ISUPPZ, WORK, &LWORK, IWORK, &LIWORK, &INFO);

    // Check if the decomposition failed (i.e., if INFO != 0)
    if (INFO != 0) {
        VL = 0.0;
        ABSTOL = 0.0;
        dsyevr_(&JOBZ, &RANGE, &UPLO, &N, X, &LDX, &VL, &VU, &IL, &IU, &ABSTOL, 
                &M, W, Z, &LDZ, ISUPPZ, WORK, &LWORK, IWORK, &LIWORK, &INFO);

        if (INFO != 0) {
            fprintf(stderr, "Error: eigenvalue computation failed\n");
            exit(1);
        }
    }

    // Compute Z = Z*Diag(W)^{1/2}
    int INCZ = 1;
    double temp;
    for (int j = 0; j < M; ++j) {

        // Scale j-th column of Z by sqrt(W[j])
        temp = sqrt(W[j]);
        dscal_(&N, &temp, Z + j*N, &INCZ);
    }

    /* Compute X_+ = ALPHA*Z*Z^T + BETA*X = Z*Z^T
     * Only writes to upper-triangular part of X.
     * When M = 0, we will obtain X = 0.
     */
    char TRANS = 'N';
    ALPHA = 1.0;
    double BETA = 0.0;

    dsyrk_(&UPLO, &TRANS, &N, &M, &ALPHA, Z, &LDZ, &BETA, X, &LDX);

    /* Compute Q = -X_- = X_+ - X */
    ALPHA = 1.0;
    daxpy_(&nn, &ALPHA, X, &INCX, Q, &INCQ);

}
