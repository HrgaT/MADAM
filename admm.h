#ifndef ADMM_H
#define ADMM_H

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

/************************************************************************************************************/

// MESSAGES for MPI
typedef enum Message{              
    SEND_FREEWORKERS,           // send ranks of free workers
    IDLE,                       // worker is free, his local queue of subproblems is empty
    NEW_VALUE                   // better lower bound found
} Message;

// TAGS in MPI messages
typedef enum Tags{
    OVER,                       // info to finish
    MESSAGE,                    // type of message
    FREEWORKER,                 // when receiving/sending rank of free worker
    NUM_FREE_WORKERS,
    PROBLEM,
    LOWER_BOUND,                // new lower bound
    SOLUTION                    // solution vector
} Tags;

/************************************************************************************************************/


#define BIG_NUMBER 1e+9

/* Maximum number of cutting planes (triangle, pentagonal and heptagonal inequalities) allowed to add */
#define MaxTriIneqAdded 50000
#define MaxPentIneqAdded 50000
#define MaxHeptaIneqAdded 50000

/* Branching strategies */
#define LEAST_FRACTIONAL  0
#define MOST_FRACTIONAL   1

/* macros for allocating vectors and matrices */
#define alloc_vector(var, size, type)\
    var = (type *) calloc((size) , sizeof(type));\
    if(var == NULL){\
        fprintf(stderr,\
                "\nError: Memory allocation problem for variable "#var" in %s line %d\n",\
                __FILE__,__LINE__);\
        exit(1);\
    }

#define alloc(var, type) alloc_vector(var, 1, type)
#define alloc_matrix(var, size, type) alloc_vector(var, (size)*(size), type)


// Parameters and their default values
#ifndef PARAM_FIELDS
#define PARAM_FIELDS \
    P(double,   sigma0,              "%lf",              1.0) \
    P(double,   scaleSigma,          "%lf",            1.001) \
    P(int,      ADMM_itermax,        "%d",              2000) \
    P(double,   tol0,                "%lf",             5e-2) \
    P(double,   scaleTol,            "%lf",              0.6) \
    P(double,   minTol,              "%lf",             1e-4) \
    P(int,      triag_iter,          "%d",                 5) \
    P(int,      pent_iter,           "%d",                 5) \
    P(int,      hept_iter,           "%d",                 5) \
    P(int,      max_outer_iter,      "%d",                20) \
    P(int,      extra_iter,          "%d",                10) \
    P(double,   violated_TriIneq,    "%lf",             1e-3) \
    P(int,      TriIneq,             "%d",              5000) \
    P(int,      adjust_TriIneq,      "%d",                 1) \
    P(int,      PentIneq,            "%d",              5000) \
    P(int,      HeptaIneq,           "%d",              5000) \
    P(int,      Pent_Trials,         "%d",                60) \
    P(int,      Hepta_Trials,        "%d",                50) \
    P(int,      include_Pent,        "%d",                 1) \
    P(int,      include_Hepta,       "%d",                 1) \
    P(int,      root,                "%d",                 0) \
    P(int,      use_diff,            "%d",                 1) \
    P(int,      branchingStrategy,   "%d",   MOST_FRACTIONAL)  
#endif


typedef struct Parameters { 
#define P(type, name, format, def_value) type name;
    PARAM_FIELDS
#undef P
} Parameters;


/* Structure for storing triangle inequalities */
typedef struct Triangle_Inequality {
    int i;
    int j;
    int k;
    int type;       // type: 1-4
    double value;   // cut violation 
    double y;       // corresponding dual multiplier
} Triangle_Inequality;

/* Structure for storing pentagonal inequalities */
typedef struct Pentagonal_Inequality {
    int type;           // type: 1-3 (based on H1 = ee^T, ...)
    int permutation[5];
    double value;       // cut violation 
    double y;           // corresponding dual multiplier
} Pentagonal_Inequality;


/* Structure for storing heptagonal inequalities */
typedef struct Heptagonal_Inequality {
    int type;           // type: 1-4 (based on H1 = ee^T, ...)
    int permutation[7];
    double value;       // cut violation 
    double y;           // corresponding dual multiplier
} Heptagonal_Inequality;


/* The main problem and any subproblems are stored using the following structure. */
typedef struct Problem {
    double *L;          // Objective matrix 
    int n;              // size of L
    int NIneq;          // number of triangle inequalities
    int NPentIneq;      // number of pentagonal inequalities
    int NHeptaIneq;     // number of heptagonal inequalities
} Problem;


/* Maximum number of variables */
#define NMAX 1024

/* Solution of the problem */
typedef struct BabSolution {
    /*
     * Vector X: Binary vector that stores the solution of the branch-and-bound algorithm
     */
    int X[NMAX];
} BabSolution;


/*
 * Node of the branch-and-bound tree.
 * Structure that represent a node of the branch-and-bound tree and stores all the 
 * useful information.
 */
typedef struct BabNode {
    int xfixed[NMAX];       // 0-1 vector specifying which nodes are fixed
    BabSolution sol;        // 0-1 solution vector
    double fracsol[NMAX];   // fractional vector obtained from primal matrix X (last column except last element)
                            // from bounding routine. Used for determining the next branching variable.
    int level;              // level (depth) of the node in B&B tree     
    double upper_bound;     // upper bound on solution value of max-cut, i.e. MC <= upper_bound.
                            // Used for determining the next node in priority queue.  
} BabNode;


/* heap (data structure) declaration */
typedef struct Heap {   
    int size;               /* maximum number of elements in heap */
    int used;               /* current number of elements in heap */
    BabNode** data;         /* array of BabNodes                  */  
} Heap;



/****** BLAS  ******/

// level 1 blas 
extern void dscal_(int *n, double *alpha, double *X, int *inc);
extern void dcopy_(int *n, double *X, int *incx, double *Y, int *incy);
extern double dnrm2_(int *n, double *x, int *incx);
extern void daxpy_(int *n, double *alpha, double *X, int *incx, double *Y, int *incy);
extern double ddot_(int *n, double *X, int *incx, double *Y, int *incy);

// level 2 blas
extern void dsymv_(char *uplo, int *n, double *alpha, double *A, int *lda, double *x, int *incx, double *beta, double *y, int *incy);
extern void dsyr_(char *uplo, int *n, double *alpha, double *x, int *incx, double *A, int *lda);

// level 3 blas
extern void dsymm_(char *side, char *uplo, int *m, int *n, double *alpha, double *A, int *lda, double *B, int *ldb, double *beta, double *C, int *ldc);
extern void dsyrk_(char *UPLO, char *TRANS, int *N, int *K, double *ALPHA, double *A, int *LDA, double *BETA, double *C, int *LDC);


/****** LAPACK  ******/

// computes matrix norm
extern double dlansy_(char *NORM, char *UPLO, int *N, double *A, int *LDA, double *WORK);

// computes Cholesky factorization of positive definite matrix
extern void dpotrf_(char *uplo, int *n, double *X, int *lda, int *info);

// computes the inverse of a real symmetric positive definite
// matrix  using the Cholesky factorization  
extern void dpotri_(char *uplo, int *n, double *X, int *lda, int *info);

// computes solution to a real system of linear equations with positive definite matrix  
extern void dposv_(char *uplo, int *n, int *nrhs, double *A, int *lda, double *B, int *ldb, int *info);

// computes selected eigenvalues/eigenvectors of symmetric matrix
extern void dsyevr_(
        char *JOBZ, char *RANGE, char *UPLO, int *N, double *A, int *LDA, 
        double *VL, double *VU, int *IL, int *IU, double *ABSTOL, int *M, 
        double *W, double *Z, int *LDZ, int *ISUPPZ, double *WORK, int *LWORK, 
        int *IWORK, int *LIWORK, int *INFO);



/**** Declarations of functions per file ****/

/* admm.c */
void ADMM_solver(const Problem *PP, double *X, double *y, double *Q, double *sigma, double tol, int *nbit, int print);

/* allocate_free.c */
void allocMemory(void);
void freeMemory(void);

/* bab_functions.c */
void initializeBabSolution(void);
int Init_PQ(void);
int Bab_Init(int argc, char **argv, int rank);
double evaluateSolution(int *sol);
int updateSolution(int *x);
void master_Bab_Main(Message message, int source, int *busyWorkers, int numbWorkers, int *numbFreeWorkers, MPI_Datatype BabSolutiontype);
void worker_Bab_Main(MPI_Datatype BabSolutiontype, MPI_Datatype BabNodetype, int rank);
double time_CPU(void);
void printSolution(FILE *file);
void printFinalOutput(FILE *file, int num_nodes);
void Bab_End(void);
int getBranchingVariable(BabNode *node);
int countFixedVariables(BabNode *node);

/* bounding.c */
double SDPbound(BabNode *node, Problem *SP, Problem *PP, int rank);

/* cutting_planec.c */
double evaluateTriangleInequality(double *XX, int N, int type, int ii, int jj, int kk);
double getViolated_TriangleInequalities(double *X, int N, Triangle_Inequality *List, int *ListSize);
double updateTriangleInequalities(Problem *PP, double *y, int *NumAdded, int *NumSubtracted);
double getViolated_PentagonalInequalities(double *X, int N, Pentagonal_Inequality *Pent_List, int *ListSize);
double updatePentagonalInequalities(Problem *PP, double *y, int *NumAdded, int *NumSubtracted, int triag);
double getViolated_HeptagonalInequalities(double *X, int N, Heptagonal_Inequality *Hepta_List, int *ListSize);
double updateHeptagonalInequalities(Problem *PP, double *y, int *NumAdded, int *NumSubtracted, int hept_index);

/* evaluate.c */
double Evaluate(BabNode *node, Problem *SP, Problem *PP, int rank);
void createSubproblem(BabNode *node, Problem *SP, Problem *PP);
double getFixedValue(BabNode *node, Problem *SP);

/* heap.c */
double Bab_LBGet(void);                                 // returns global lower bound
int Bab_numEvalNodes(void);                             // returns number of evaluated nodes
void Bab_incEvalNodes(void);                            // increment the number of evaluated nodes
int isPQEmpty(void);                                    // checks if queue is empty
int Bab_LBUpd(double new_lb, BabSolution *bs);          // checks and updates lower bound if better found, returns 1 if success
BabNode* newNode(BabNode *parentNode);                  // create child node from parent
BabNode* Bab_PQPop(void);                               // take and remove the node with the highest priority
void Bab_PQInsert(BabNode *node);                       // insert node into priority queue based on intbound and level 
void Bab_LBInit(double lowerBound, BabSolution *bs);    // initialize global lower bound and solution vector
Heap* Init_Heap(int size);                              // allocates space for heap (array of BabNode*)

/* heuristic.c */
double runHeuristic(Problem *P0, Problem *P, BabNode *node, int *x);
double GW_heuristic(Problem *P0, Problem *P, BabNode *node, int *x, int num);
double mc_1opt(int *x, Problem *P0);
int update_best(int *xbest, int *xnew, double *best, Problem *P0);

/* ipm_mc_pk.c */
void ipm_mc_pk(double *L, int n, double *X, double *y, double *Z, double *phi, int print);

/* operators.c */
void diag(const double *X, double *y, int n);
void Diag(double *X, const double *y, int n);
void op_B(const Problem *P, double *y, double *X);
void op_Bt(const Problem *P, double *X, const double *tt);
void op_Diag(const Problem *P, double *X, const double *yy);
void project_SDP(const Problem *P);

/* process_input.c */
void print_symmetric_matrix(double *Mat, int N);
int processCommandLineArguments(int argc, char **argv, int rank);
int readData(const char *instance);
int readParameters(const char *path, int rank);

/* qap_simuted_annealing.c */
double qap_simulated_annealing(int *H, int k, double *X, int n, int *pent);

#endif /*ADMM_H */
