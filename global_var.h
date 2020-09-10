// NOTE: only included in allocate_free.c !

/********************************************************/
/************ List of all global variables **************/
/********************************************************/
Parameters params;                  
FILE *output;                       // output file
Problem *SP;                        // original problem instance
Problem *PP;                        // subproblem instance
double root_bound;                  // SDP upper bound at root node
double TIME;                        // CPU time
int stopped = 0;                    // true if the algorithm stopped at root node
double diff;			            // difference between basic SDP relaxation and bound with added cutting planes  
/********************************************************/


/********************************************************/
/*************** Specific to node ***********************/
/********************************************************/
/* PRIMAL variables */
double *X;                          // Stores current (psd) X (primal solution). Violated inequalities are computed from X.
double *u;                          // (nonnegative) slack variable for triangle inequalities B(X) + u = e 

/* DUAL variables */
double *y;                          // (free) dual multiplier for diagonal constraints
double *Q;                          // dual (psd) matrix to constraint X >= 0
double *t;                          // dual multiplier for triangle inequality constraints (unconstrained due to slack variable u)
double *s;                          // dual multiplier to constraint u >= 0
                          
double f;                           // objective value of relaxation                      

/* Triangle Inequalities variables */
Triangle_Inequality *Cuts;          // vector (MaxTriIneqAdded) of current triangle inequality constraints
Triangle_Inequality *List;          // vector (params.TriIneq) of new violated triangle inequalities

/* Pentagonal Inequalities variables */
Pentagonal_Inequality *Pent_Cuts;   // vector (MaxPentIneqAdded) of current pentagonal inequality constraints
Pentagonal_Inequality *Pent_List;   // vector (params.PentIneq) of new violated pentagonal inequalities

/* Heptagonal Inequalities variables */
Heptagonal_Inequality *Hepta_Cuts;   // vector (MaxHeptaIneqAdded) of current heptagonal inequality constraints
Heptagonal_Inequality *Hepta_List;   // vector (params.HeptaIneq) of new violated heptagonal inequalities

/* dsyevr function variables (projection onto PSD cone) */
int M = 0;                          // number of eigenvalues
double *W;                          // contains the eigenvalues 
double *Z;                          // contains the eigenvectors
int *ISUPPZ;                        // dim = 2*max(1,M), support of the eigvecs. in Z
double *WORK;                       // dim = LWORK
int *IWORK;                         // dim = LIWORK
int sizeWORK, sizeIWORK;            // sizes of WORK and IWORK
