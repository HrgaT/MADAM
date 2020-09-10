#include <math.h>
#include "admm.h"
#include "cholmod.h"                       // sparse Cholesky solver 

extern Parameters params;
extern double TIME;
extern double f;                           // opt. value of SDP 
extern int M;                              // number of eigenvalues
extern double *Z;                          // contains the eigenvectors

/* primal and dual variables corresponding to cutting planes */
extern double *u;
extern double *t;
extern double *s;

/* dsyevr variables */
extern int *ISUPPZ;         
extern double *WORK;        
extern int *IWORK;          
extern int sizeWORK, sizeIWORK;

/* vector of triangle inequality constraints */
extern Triangle_Inequality *Cuts;   
extern Pentagonal_Inequality *Pent_Cuts;
extern Heptagonal_Inequality *Hepta_Cuts; 

/* macro to handle the errors if Blas/Lapack routine failed */
#define ERROR_INFO(cond,message)\
        if ((cond)) {\
            fprintf(stderr, "\nError: "#message"\n");\
            exit(1);\
        }

/******************** ADMM solver *********************/
/* ADMM method for solving Max-Cut SDP relaxation strengthened
 * with triangle inequality constraints.
 *
 * Additional variable s is introduced in the dual SDP to
 * avoid solving QP for t. Now one linear system for t is
 * solved, projection onto PSD cone to obtain Q and X,
 * projection onto nonnegative orthant to obtain s + additional dual update to get u.
 *
 * variables:
 *   DUAL: y (free), t (free), Q (psd), s (nonnegative)
 * PRIMAL: X (dual variable for equality constraint L-Diag(y)-BT(t)+Z = 0)
 *         u slack variable for triangle inequalities B(X) <= e
 ******************************************************/

void ADMM_solver(const Problem *PP, double *X, double *y, double *Q, double *sigma, double tol, int *nbit, int print) {

    /* constraints: 
     * n equality constraints : diag(X) = e 
     * PP->NIneq + PP->NPentIneq + PP->NHeptaIneq inequality constraints: B(X) <= e */

     // size of problem
    int n = PP->n;   
    int nn = n*n;  

    // number of cutting planes
    int m = PP->NIneq + PP->NPentIneq + PP->NHeptaIneq; 
    
    int inc = 1;        // for BLAS/LAPACK
    double alpha;
    char NORM = 'F';    // Frobenious norm
    char UPLO = 'L';    // only upper triangular part is referenced!

    double primal;
    double dual;
    double gap;
    double err_p;       // primal error
    double err_d;       // dual error
    double err;         // overall error
    double digits_d;
    double digits_p;
    double secs; 
    

    /**** Initialize t, s, and u ****/
    for (int ineq = 0; ineq < PP->NIneq; ++ineq) {
        s[ineq] = Cuts[ineq].y;
        t[ineq] = Cuts[ineq].y;
        u[ineq] = 0.0;
    }

    for (int ineq = 0; ineq < PP->NPentIneq; ++ineq) {
        s[ineq + PP->NIneq] = Pent_Cuts[ineq].y;
        t[ineq + PP->NIneq] = Pent_Cuts[ineq].y;
        u[ineq + PP->NIneq] = 0.0;
    }

    for (int ineq = 0; ineq < PP->NHeptaIneq; ++ineq) {
        s[ineq + PP->NIneq + PP->NPentIneq] = Hepta_Cuts[ineq].y;
        t[ineq + PP->NIneq + PP->NPentIneq] = Hepta_Cuts[ineq].y;
        u[ineq + PP->NIneq + PP->NPentIneq] = 0.0;
    }

    /******************** using CHOLMOD ********************/
    cholmod_dense *b_cholmod = NULL;
    cholmod_dense *t_cholmod = NULL;
    cholmod_factor *L_cholmod;
    cholmod_sparse *B_cholmod;
    cholmod_triplet *triplet;
    double beta[2] = {1, 0};

    /* ---------------------------------------------------------------------- */
    /* start CHOLMOD and set parameters */
    /* ---------------------------------------------------------------------- */
    cholmod_common Common;
    cholmod_start(&Common);

    /* ---------------------------------------------------------------------- */
    /* create matrix B */
    /* ---------------------------------------------------------------------- */

    // allocate space for triple: m rows, n^2 columns, 6m + ... nonzero elements, 0 = unsymmetric
    // nonzero elements (3 for triag + 10 for penta + 21 for hepta) x 2 for symmetrix matrices (lower part) 
    triplet = cholmod_allocate_triplet(m, n*n, 6*PP->NIneq + 20*PP->NPentIneq + 42*PP->NHeptaIneq, 0, CHOLMOD_REAL, &Common); // all pointers already initialized
    
    // construct triplet using Cuts array
    int ii, jj, kk, ll, mm, nnn, oo, type, l;

    /* triangles */
    for (int ineq = 0; ineq < PP->NIneq; ++ineq) {

        type = Cuts[ineq].type;
        ii   = Cuts[ineq].i;
        jj   = Cuts[ineq].j;
        kk   = Cuts[ineq].k;

        /*** rows ***/
        for (l = 6*ineq; l < 6*(ineq+1); ++l){
            ((int*)triplet->i)[l] = ineq;
        } 

        /*** columns ***/
        ((int*)triplet->j)[6*ineq]   = jj * n + ii;
        ((int*)triplet->j)[6*ineq+1] = ii * n + jj;
        ((int*)triplet->j)[6*ineq+2] = ii * n + kk;
        ((int*)triplet->j)[6*ineq+3] = kk * n + ii;
        ((int*)triplet->j)[6*ineq+4] = kk * n + jj;
        ((int*)triplet->j)[6*ineq+5] = jj * n + kk;

        /*** values ***/
        if (type == 1){
            for (l = 6*ineq; l < 6*(ineq+1); ++l){
                ((double*)triplet->x)[l] = -0.5;
            } 
        }
        else if (type == 2){
            ((double*)triplet->x)[6*ineq]   = -0.5;
            ((double*)triplet->x)[6*ineq+1] = -0.5;
            ((double*)triplet->x)[6*ineq+2] = 0.5;
            ((double*)triplet->x)[6*ineq+3] = 0.5;
            ((double*)triplet->x)[6*ineq+4] = 0.5;
            ((double*)triplet->x)[6*ineq+5] = 0.5;
        }
        else if (type == 3){
            ((double*)triplet->x)[6*ineq]   = 0.5;
            ((double*)triplet->x)[6*ineq+1] = 0.5;
            ((double*)triplet->x)[6*ineq+2] = -0.5;
            ((double*)triplet->x)[6*ineq+3] = -0.5;
            ((double*)triplet->x)[6*ineq+4] = 0.5;
            ((double*)triplet->x)[6*ineq+5] = 0.5;
        }
        else if (type == 4){
            ((double*)triplet->x)[6*ineq]   = 0.5;
            ((double*)triplet->x)[6*ineq+1] = 0.5;
            ((double*)triplet->x)[6*ineq+2] = 0.5;
            ((double*)triplet->x)[6*ineq+3] = 0.5;
            ((double*)triplet->x)[6*ineq+4] = -0.5;
            ((double*)triplet->x)[6*ineq+5] = -0.5;
        }
    }

    /* pentagonals */
    for (int ineq = 0; ineq < PP->NPentIneq; ++ineq) {

        type = Pent_Cuts[ineq].type;
        ii   = Pent_Cuts[ineq].permutation[0];
        jj   = Pent_Cuts[ineq].permutation[1];
        kk   = Pent_Cuts[ineq].permutation[2];
        ll   = Pent_Cuts[ineq].permutation[3];
        mm   = Pent_Cuts[ineq].permutation[4];

        /*** rows ***/
        for (l = 6*PP->NIneq + 20*ineq; l < 6*PP->NIneq + 20*(ineq+1); ++l){
            ((int*)triplet->i)[l] = PP->NIneq + ineq;
        } 

        /*** columns ***/
        ((int*)triplet->j)[6*PP->NIneq + 20*ineq]   = jj * n + ii;
        ((int*)triplet->j)[6*PP->NIneq + 20*ineq+1] = ii * n + jj;
        ((int*)triplet->j)[6*PP->NIneq + 20*ineq+2] = ii * n + kk;
        ((int*)triplet->j)[6*PP->NIneq + 20*ineq+3] = kk * n + ii;
        ((int*)triplet->j)[6*PP->NIneq + 20*ineq+4] = ii * n + ll;
        ((int*)triplet->j)[6*PP->NIneq + 20*ineq+5] = ll * n + ii;
        ((int*)triplet->j)[6*PP->NIneq + 20*ineq+6] = ii * n + mm;
        ((int*)triplet->j)[6*PP->NIneq + 20*ineq+7] = mm * n + ii;
        ((int*)triplet->j)[6*PP->NIneq + 20*ineq+8] = jj * n + kk;
        ((int*)triplet->j)[6*PP->NIneq + 20*ineq+9] = kk * n + jj;
        ((int*)triplet->j)[6*PP->NIneq + 20*ineq+10] = jj * n + ll;
        ((int*)triplet->j)[6*PP->NIneq + 20*ineq+11] = ll * n + jj;
        ((int*)triplet->j)[6*PP->NIneq + 20*ineq+12] = jj * n + mm;
        ((int*)triplet->j)[6*PP->NIneq + 20*ineq+13] = mm * n + jj;
        ((int*)triplet->j)[6*PP->NIneq + 20*ineq+14] = kk * n + ll;
        ((int*)triplet->j)[6*PP->NIneq + 20*ineq+15] = ll * n + kk;
        ((int*)triplet->j)[6*PP->NIneq + 20*ineq+16] = kk * n + mm;
        ((int*)triplet->j)[6*PP->NIneq + 20*ineq+17] = mm * n + kk;
        ((int*)triplet->j)[6*PP->NIneq + 20*ineq+18] = ll * n + mm;
        ((int*)triplet->j)[6*PP->NIneq + 20*ineq+19] = mm * n + ll;

        /*** values ***/
        if (type == 1){
            for (l = 6*PP->NIneq + 20*ineq; l < 6*PP->NIneq + 20*(ineq+1); ++l){
                ((double*)triplet->x)[l] = -0.25;
            } 
        }
        else if (type == 2){

            ((double*)triplet->x)[6*PP->NIneq + 20*ineq] = 0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 1] = 0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 2] = 0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 3] = 0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 4] = 0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 5] = 0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 6] = 0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 7] = 0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 8] = -0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 9] = -0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 10] = -0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 11] = -0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 12] = -0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 13] = -0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 14] = -0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 15] = -0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 16] = -0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 17] = -0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 18] = -0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 19] = -0.25;

        }
        else if (type == 3){
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq] = -0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 1] = -0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 2] = 0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 3] = 0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 4] = 0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 5] = 0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 6] = 0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 7] = 0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 8] = 0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 9] = 0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 10] = 0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 11] = 0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 12] = 0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 13] = 0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 14] = -0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 15] = -0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 16] = -0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 17] = -0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 18] = -0.25;
            ((double*)triplet->x)[6*PP->NIneq + 20*ineq + 19] = -0.25;
        }
    
    }

    /* heptagonal */
    for (int ineq = 0; ineq < PP->NHeptaIneq; ++ineq) {

        type = Hepta_Cuts[ineq].type;
        ii   = Hepta_Cuts[ineq].permutation[0];
        jj   = Hepta_Cuts[ineq].permutation[1];
        kk   = Hepta_Cuts[ineq].permutation[2];
        ll   = Hepta_Cuts[ineq].permutation[3];
        mm   = Hepta_Cuts[ineq].permutation[4]; 
        nnn   = Hepta_Cuts[ineq].permutation[5];
        oo   = Hepta_Cuts[ineq].permutation[6];   

        /*** rows ***/
        for (l = 6*PP->NIneq + 20*PP->NPentIneq + 42*ineq; l < 6*PP->NIneq + 20*PP->NPentIneq + 42*(ineq+1); ++l){
            ((int*)triplet->i)[l] = PP->NIneq + PP->NPentIneq + ineq;
        }

        /*** columns ***/
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq]       = jj * n + ii;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 1]   = ii * n + jj;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 2]   = kk * n + ii;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 3]   = ii * n + kk;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 4]   = ll * n + ii;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 5]   = ii * n + ll;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 6]   = mm * n + ii;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 7]   = ii * n + mm;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 8]   = nnn * n + ii;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 9]   = ii * n + nnn;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 10]   = oo * n + ii;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 11]   = ii * n + oo;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 12]   = jj * n + kk;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 13]   = kk * n + jj;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 14]   = jj * n + ll;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 15]   = ll * n + jj;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 16]   = jj * n + mm;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 17]   = mm * n + jj;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 18]   = jj * n + nnn;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 19]   = nnn * n + jj;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 20]   = jj * n + oo;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 21]   = oo * n + jj;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 22]   = kk * n + ll;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 23]   = ll * n + kk;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 24]   = kk * n + mm;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 25]   = mm * n + kk;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 26]   = kk * n + nnn;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 27]   = nnn * n + kk;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 28]   = kk * n + oo;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 29]   = oo * n + kk;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 30]   = ll * n + mm;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 31]   = mm * n + ll;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 32]   = ll * n + nnn;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 33]   = nnn * n + ll;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 34]   = ll * n + oo;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 35]   = oo * n + ll;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 36]   = mm * n + nnn;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 37]   = nnn * n + mm;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 38]   = mm * n + oo;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 39]   = oo * n + mm;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 40]   = nnn * n + oo;
        ((int*)triplet->j)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 41]   = oo * n + nnn;
        

        /*** values ***/
        if (type == 1){
            for (l = 6*PP->NIneq + 20*PP->NPentIneq + 42*ineq; l < 6*PP->NIneq + 20*PP->NPentIneq + 42*(ineq+1); ++l){
                ((double*)triplet->x)[l] = -1.0/6.0;
            } 
        }
        else if (type == 2){

            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq]       = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 1]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 2]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 3]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 4]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 5]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 6]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 7]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 8]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 9]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 10]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 11]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 12]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 13]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 14]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 15]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 16]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 17]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 18]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 19]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 20]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 21]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 22]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 23]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 24]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 25]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 26]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 27]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 28]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 29]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 30]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 31]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 32]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 33]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 34]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 35]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 36]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 37]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 38]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 39]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 40]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 41]   = -1.0/6.0;

        }    
        else if (type == 3){

            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq]       = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 1]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 2]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 3]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 4]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 5]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 6]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 7]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 8]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 9]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 10]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 11]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 12]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 13]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 14]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 15]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 16]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 17]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 18]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 19]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 20]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 21]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 22]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 23]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 24]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 25]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 26]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 27]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 28]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 29]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 30]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 31]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 32]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 33]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 34]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 35]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 36]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 37]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 38]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 39]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 40]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 41]   = -1.0/6.0;

        }
        else if (type == 4){

            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq]       = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 1]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 2]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 3]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 4]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 5]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 6]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 7]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 8]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 9]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 10]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 11]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 12]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 13]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 14]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 15]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 16]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 17]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 18]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 19]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 20]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 21]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 22]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 23]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 24]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 25]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 26]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 27]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 28]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 29]   = 1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 30]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 31]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 32]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 33]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 34]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 35]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 36]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 37]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 38]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 39]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 40]   = -1.0/6.0;
            ((double*)triplet->x)[6*PP->NIneq + 20*PP->NPentIneq + 42*ineq + 41]   = -1.0/6.0;

        }

    }


    triplet->nnz = 6*PP->NIneq + 20*PP->NPentIneq + 42*PP->NHeptaIneq;

    // triplet to sparse and allocate space for rhs of linear system (BBt + I)t = b_cholmod
    B_cholmod = cholmod_triplet_to_sparse(triplet, 6*PP->NIneq + 20*PP->NPentIneq + 42*PP->NHeptaIneq, &Common);
    b_cholmod = cholmod_zeros(m, 1, B_cholmod->xtype, &Common);

    /* ---------------------------------------------------------------------- */
    /* analyze and factorize */
    /* ---------------------------------------------------------------------- */
    L_cholmod = cholmod_analyze(B_cholmod, &Common);
    cholmod_factorize_p(B_cholmod, beta, NULL, 0, L_cholmod, &Common); // factorize BBt + I
    /*******************************************************/


    // counts number of ADMM iterations
    *nbit = 0; 
    
    /* print output */
    if (print) {
        puts("======================================================================================");
        printf("%4s  %6s  %8s  %8s  %5s  %9s  %9s  %12s\n", 
            "iter", "time", "primal", "dual", "gap", "log10(rp)", "log10(rd)", "log10(sigma)");
        puts("======================================================================================"); 
    }

    /* some temporary variables */
    double *tmp;
    alloc_matrix(tmp, n, double);

    double *rhs;
    alloc_vector(rhs, m, double);

    double *res;
    alloc_vector(res, n, double);

    double *v;
    alloc_vector(v, m, double);


    // NOTE: we use the fact that operators diag and B are orthogonal
    // --> diag(Bt(t)) = 0 and B(Diag(y)) = 0

    /***** main loop *****/
    while ( (*nbit) <= params.ADMM_itermax ) {

        // increase iteration count
        (*nbit)++;

        /* compute y = diag(L + Q + X/sigma) - e/sigma */
        for (int i = 0; i < n; ++i) {
            y[i] = PP->L[i + i * n] + Q[i + i * n] + (X[i + i * n] - 1.0) / (*sigma);
        }

        /* compute tmp = L + Q + X/sigma */
        for (int i = 0; i < nn; ++i) {
            tmp[i] = PP->L[i] + Q[i] + X[i]/(*sigma);
        }

        /*** compute t ***/
        // rhs = B(tmp) + s + (u-b)/sigma --> copy into t
        dcopy_(&m,s,&inc,t,&inc);                 // copy s to t
        op_B(PP, t, tmp);                         // t += B(tmp)

        for (int i = 0; i < m; ++i) {
            t[i] += (u[i] - 1) / (*sigma);        // t += (u-b)/sigma
            ((double*)b_cholmod->x)[i] = t[i];    // copy t into cholmod object
        }

        /*** CHOLMOD: solve linear system (BBt + I)t = rhs ***/
        t_cholmod = cholmod_solve(CHOLMOD_A, L_cholmod, b_cholmod, &Common);

        // copy cholmod solution to t
        for (int i = 0; i < m; ++i) {
            t[i] = ((double*)t_cholmod->x)[i];
        }

        // free memory for t_cholmod --> otherwise memory leak
        cholmod_free_dense(&t_cholmod, &Common);
        
        /**** compute s = (t-u/sigma)_+ = v_+ ****/
        for (int i = 0; i < m; ++i) {
            v[i] = t[i] - u[i]/(*sigma);
            s[i] = ( ( v[i] > 0 ) ? v[i] : 0 );
        }

        /**** compute W = L + X/sigma - Diag(y) - Bt(t) ****/
        // NOTE: save it in X!
        for (int i = 0; i < nn; ++i) {
            X[i] = PP->L[i] + X[i]/(*sigma);
        }
        op_Diag(PP, X, y);
        op_Bt(PP, X, t);          
        
        /* 
         * projection onto PSD cone:
         * computes X_+ and Q = -(X_-)
         */ 
        project_SDP(PP);  


        /**** update primal variables ****/

        // Scale: X = sigma * X
        dscal_(&nn, sigma, X, &inc);

        // Update Z so that X = Z*Z'
        double sca = sqrt(*sigma);
        int nM = n * M;
        dscal_(&nM, &sca, Z, &inc);  // Z = Z * sqrt(sigma)

        // update: u = u + sigma * (s-t) 
        //           = -sigma * v_- 
        for (int i = 0; i < m; ++i) {
            u[i] = ( ( v[i] < 0 ) ? -(*sigma) * v[i] : 0 );
        }


        /* 
         * compute inner (primal) error || Diag(X) - e || + || max(B(X)-e,0) ||)
         * and divide it by 1 + sqrt(n) 
         */

        // || Diag(X) - e ||
        diag(X, res, n);
        for (int i = 0; i < n; ++i) {
            res[i] -= 1.0;
        }

        err_p = dnrm2_(&n, res, &inc);

        // || max(B(X)-b,0) ||
        for (int i = 0; i < m; ++i) {
            rhs[i] = -1.0;         
        }

        // rhs += B(X)
        op_B(PP, rhs, X);                

        for (int i = 0; i < m; ++i) {
            rhs[i] = ( (rhs[i] > 0) ? rhs[i] : 0);
        }

        err_p += dnrm2_(&m, rhs, &inc);
        err_p /= (1 + sqrt(n));

        /* 
         * compute outer (dual) error || Z - Diag(y) - Bt(t) + L || + || s - t ||
         * and divide it by 1 + norm(L, 'fro')
         */
        dcopy_(&m,s,&inc,rhs,&inc);                 // copy s to rhs
        alpha = -1.0;
        daxpy_(&m,&alpha,t,&inc,rhs,&inc);          // rhs -= t

        err_d = dnrm2_(&m,rhs,&inc);

        // tmp = L + Z - Diag(y) - Bt(t)
        for (int i = 0; i < nn; ++i) {
            tmp[i] = PP->L[i] + Q[i];
        }
        op_Diag(PP, tmp, y);
        op_Bt(PP, tmp, t);

        err_d += dlansy_(&NORM, &UPLO, &n, tmp, &n, WORK);
        err_d /= (1 + dlansy_(&NORM, &UPLO, &n, PP->L, &n, WORK));


        /**** compute primal and dual objective values ****/
        // NOTE: do not use ddot_(&nn,PP->L,&inc,X,&inc)
        // --> lower  triangular part of X contains garbage
        // and should not be referenced!
        primal = 0.0;
        for (int i = 0; i < n; ++i) {
            for (int j = i; j < n; ++j) {
                if (i==j) {
                    primal += PP->L[i + i*n] * X[i + i*n];
                }
                else {
                    primal += 2 * PP->L[j + i*n] * X[j + i*n];
                }
            }
        }

        dual = 0.0;
        for (int i = 0; i < n; ++i) {
            dual += y[i];
        }
        for (int ineq = 0; ineq < m; ++ineq) {
            dual += t[ineq]; 
        }

        /**** duality gap ****/
        gap = fabs(dual - primal) / (1 + fabs(dual) + fabs(primal));

        /*** error = max(err_p, err_d, gap) ***/
        err = ( err_p > err_d ? (err_p > gap ? err_p : gap) : (err_d > gap ? err_d : gap) ); 

        /**** time ****/
        secs = MPI_Wtime() - TIME;

        /**** number of digits ****/
        digits_p = log10(err_p);
        digits_d = log10(err_d);

        /**** print output (every 10 iterations) ****/
        if (print) {
            if ( (*nbit) % 10 == 0) {
                printf("%4d  %6.2f  %8.2f  %8.2f  %5.2f  %9.2f  %9.2f  %12.2f\n", 
                        *nbit, secs, primal, dual, gap, digits_p, digits_d, log10(*sigma));
            }
        }

        /**** adjust sigma ****/
        if (digits_d - digits_p < -0.5)
            *sigma /= params.scaleSigma;
        else if (digits_d - digits_p > 0.5)
            *sigma *= params.scaleSigma;


        /* stopping condition */
        if ( err < tol ) {
            if (print) {
                printf("%4d  %6.2f  %8.2f  %8.2f  %5.2f  %9.2f  %9.2f  %12.2f\n", 
                        *nbit, secs, primal, dual, gap, digits_p, digits_d, log10(*sigma));
            }
            break;
        }

    }

    /* ---------------------------------------------------------------------- */
    /* free matrices and finish CHOLMOD */
    /* ---------------------------------------------------------------------- */
    cholmod_free_dense(&b_cholmod, &Common);
    cholmod_free_dense(&t_cholmod, &Common);
    cholmod_free_factor(&L_cholmod, &Common);
    cholmod_free_sparse (&B_cholmod, &Common);
    cholmod_free_triplet(&triplet, &Common);
    cholmod_finish (&Common);
    /* ---------------------------------------------------------------------- */


    /******************** POST-PROCESSING to obtain safe upper bound ********************/
    /* STEPS:
     * Z_new = L-diag(y)-mat(B'*s);
     * y = y + max(eig(Z_new)) * e;
     * f = sum(y) + sum(s);
     */
    double *Z_new;
    alloc_matrix(Z_new, n, double);

    dcopy_(&nn,PP->L,&inc,Z_new,&inc);     // copy L into Z_new
    op_Diag(PP, Z_new, y);
    op_Bt(PP, Z_new, s);

    /**** compute largest eigenvalue ****/
    char JOBZ = 'N';
    char RANGE = 'I';   
    int LDA = n;
    double VL, VU;      // not referenced       
    int IL = n;
    int IU = n;
    double ABSTOL = 1e-8;
    int NUM_EIG;
    double lambda_max;
    int LDZ = n;
    int LWORK = sizeWORK;
    int LIWORK = sizeIWORK;
    int INFO;
    
    dsyevr_(&JOBZ, &RANGE, &UPLO, &n, Z_new, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, 
        &NUM_EIG, &lambda_max, Z, &LDZ, ISUPPZ, WORK, &LWORK, IWORK, &LIWORK, &INFO);

    if (INFO != 0) {
        fprintf(stderr, "Error: eigenvalue computation failed\n");
        exit(1);
    }

    /*** safe bound ***/
    free(Z_new);

    f = n * lambda_max;
    for (int i = 0; i < n; ++i) {
        f += y[i];
    }
    for (int ineq = 0; ineq < m; ++ineq) {
        f += s[ineq]; 
    }    
    /*---------------------------------------*/

    free(tmp);
    free(rhs);
    free(res);  
    free(v);

}


