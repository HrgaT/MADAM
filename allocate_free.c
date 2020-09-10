#include <stdio.h>

#include "admm.h"
#include "global_var.h"

extern BabSolution *BabSol; // global solution of B&B algorithm defined in heap.c

void allocMemory(void) {

    /* 
     * SP, SP->n, SP->L, PP, PP->n and PP->L 
     * are all allocated and defined in readData (process_input.c),
     * before this function is called
     */
    int N = SP->n;

    /* triangle inequalities */
    alloc_vector(Cuts, MaxTriIneqAdded, Triangle_Inequality);
    alloc_vector(List, params.TriIneq, Triangle_Inequality);

    /* pentagonal inequalities */
    alloc_vector(Pent_Cuts, MaxPentIneqAdded, Pentagonal_Inequality);
    alloc_vector(Pent_List, params.PentIneq, Pentagonal_Inequality);

    /* heptagonal inequalities */
    alloc_vector(Hepta_Cuts, MaxHeptaIneqAdded, Heptagonal_Inequality);
    alloc_vector(Hepta_List, params.HeptaIneq, Heptagonal_Inequality);

    /* primal and dual variables */
    alloc_matrix(X, N, double);
    alloc_vector(u, MaxTriIneqAdded + MaxPentIneqAdded + MaxHeptaIneqAdded, double);
    alloc_vector(y, N, double);
    alloc_matrix(Q, N, double);
    alloc_vector(t, MaxTriIneqAdded + MaxPentIneqAdded + MaxHeptaIneqAdded, double);
    alloc_vector(s, MaxTriIneqAdded + MaxPentIneqAdded + MaxHeptaIneqAdded, double);

    /* variables for projection onto PSD cone */
    sizeWORK = 26 * N;
    sizeIWORK = 10 * N;

    alloc_vector(WORK, sizeWORK, double);
    alloc_vector(IWORK, sizeIWORK, int);
    alloc_vector(W, N, double);
    alloc_matrix(Z, N, double);
    alloc_vector(ISUPPZ, 2 * N, int);
}


void freeMemory(void) {

    free(SP->L);
    free(SP);
    free(PP->L);
    free(PP);

    free(Cuts);
    free(List);

    free(Pent_Cuts);
    free(Pent_List);

    free(Hepta_Cuts);
    free(Hepta_List);

    free(X);
    free(u);
    free(y);
    free(Q);
    free(t);
    free(s);
   
    free(W);
    free(Z);
    free(ISUPPZ);
    free(WORK);
    free(IWORK);

    free(BabSol);
}
