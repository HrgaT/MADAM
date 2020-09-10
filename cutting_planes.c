#include <math.h>
#include "admm.h"

extern Parameters params;

extern Triangle_Inequality *Cuts;        
extern Triangle_Inequality *List; 

extern Pentagonal_Inequality *Pent_Cuts;        
extern Pentagonal_Inequality *Pent_List;

extern Heptagonal_Inequality *Hepta_Cuts;        
extern Heptagonal_Inequality *Hepta_List;

extern double *X;

/************************* TRIANGLE INEQUALITIES *************************/

/* evaluate triangle inequality */
double evaluateTriangleInequality(double *XX, int N, int type, int ii, int jj, int kk)  {

    double ineqvalue;

    switch (type) {
        case 1:
            ineqvalue =  -XX[ii+jj*N] - XX[ii+kk*N] - XX[jj+kk*N] - 1.0;
            break;
        case 2:
            ineqvalue =  -XX[ii+jj*N] + XX[ii+kk*N] + XX[jj+kk*N] - 1.0;
            break;
        case 3:
            ineqvalue =   XX[ii+jj*N] - XX[ii+kk*N] + XX[jj+kk*N] - 1.0;
            break;
        default: // case 4 :
            ineqvalue =   XX[ii+jj*N] + XX[ii+kk*N] - XX[jj+kk*N] - 1.0;
    }

    return ineqvalue;
}


/* 
 * Evaluates each triangle inequality (i.e. cut) using the X matrix and
 * fills the Triangle_Inequality array named List, with at most params.TriIneq
 * inequalities that are violated by at least params.violated_TriIneq.
 * It also returns the value of the cut that is violated the most by X.
 */
double getViolated_TriangleInequalities(double *X, int N, Triangle_Inequality *List, int *ListSize) {

    int ListCount;                              // loop index
    int size = 0;                               // number of added cuts
    int LeastViolatedIneq = 0;                  // index of least violated inequality
    double LeastViolatedIneqValue = BIG_NUMBER; // minimum violation
    double maxAllIneq = -BIG_NUMBER;            // maximum violation
    double test_ineqvalue;                      // violation of current cut

    // Loop through all inequalities
    for (int type = 1; type <= 4; ++type) {

        for (int ii = 0; ii < N; ++ii) {
            for (int jj = 0; jj < ii; ++jj) {
                for (int kk = 0; kk < jj; ++kk) {

                    test_ineqvalue = evaluateTriangleInequality(X, N, type, ii, jj, kk);

                    // keep track of the maximum value of test_ineqvalue, i.e. 
                    // current most violated cut value
                    maxAllIneq = (test_ineqvalue > maxAllIneq) ? test_ineqvalue : maxAllIneq;

                    if (test_ineqvalue > params.violated_TriIneq) {

                        // (1) put first params.TriIneq violated inequalities in list and 
                        //     keep track of the least violated inequality
                        if (size < params.TriIneq) {

                            // add ineq to the end of List
                            List[size].type  = type;
                            List[size].i     = ii;
                            List[size].j     = jj;
                            List[size].k     = kk;
                            List[size].value = test_ineqvalue;

                            // update LeastViolatedIneq
                            if (size == 0 || test_ineqvalue < LeastViolatedIneqValue) {
                                LeastViolatedIneqValue = test_ineqvalue;
                                LeastViolatedIneq = size;
                            }

                            ++size;
                        } 
                        else if (test_ineqvalue > LeastViolatedIneqValue) 
                        {
                            // (2) if you find an inequality that is violated more than the 
                            //     least violated inequality, then add it to the list, 
                            //     remove the least violated inequality, and find 
                            //     the new least violated inequality in the list

                            // add ineq to the list, replacing LeastViolatedIneq
                            List[LeastViolatedIneq].type  = type;
                            List[LeastViolatedIneq].i     = ii;
                            List[LeastViolatedIneq].j     = jj;
                            List[LeastViolatedIneq].k     = kk;
                            List[LeastViolatedIneq].value = test_ineqvalue;

                            // update LeastViolatedIneq
                            LeastViolatedIneqValue = BIG_NUMBER;
                            for (ListCount = 0; ListCount < size; ++ListCount) 
                            {
                                // if ineq is less violated
                                if (List[ListCount].value < LeastViolatedIneqValue) 
                                {
                                    LeastViolatedIneqValue = List[ListCount].value;
                                    LeastViolatedIneq = ListCount;
                                }
                            }
                        }
                    }
                } // kk loop
            } // jj loop
        } // ii loop

    } // type loop

    *ListSize = size;

    return maxAllIneq;
}


/* update inequalities: purge old ones and separate new ones */
/* returns the new maximum violation of triangle inequalities */
double updateTriangleInequalities(Problem *PP, double *y, int *NumAdded, int *NumSubtracted) {

    int ineq;                   // index for inequality        
    int yindex;                 // index for dual multipliers          
    int ListCount, ListSize;
    int N = PP->n;
    
    // purge inequalities: decide which cuts that were previously added to remove
    int subtracted = 0;
    yindex = 0;
    int next_ineq = 0;

    for (ineq = 0; ineq < PP->NIneq; ++ineq) {

        // store the dual multiplier
        Cuts[ineq].y = y[yindex];
        ++yindex;

        // remove inequality if dual multiplier is small
        if (Cuts[ineq].y < 1e-5) {
            ++subtracted;
        } 
        else // keep inequality
        {
            Cuts[next_ineq].type  = Cuts[ineq].type;
            Cuts[next_ineq].i     = Cuts[ineq].i;
            Cuts[next_ineq].j     = Cuts[ineq].j;
            Cuts[next_ineq].k     = Cuts[ineq].k;
            Cuts[next_ineq].y     = Cuts[ineq].y;

            ++next_ineq;
        }
    }

    PP->NIneq -= subtracted;


    // separate new triangle inequalities
    double maxAllIneq = getViolated_TriangleInequalities(X, N, List, &ListSize);

    // Add List to Cuts
    int added = 0;
    for (ListCount = 0; ListCount < ListSize; ++ListCount) { 

        // Stop if we have reached the maximum number of cuts we can add
        if (next_ineq == MaxTriIneqAdded)
            break;

        // Check if inequality is already included in Cuts
        int found_ineq = 0;
        for (ineq = 0; ineq < PP->NIneq; ++ineq) {

            if (Cuts[ineq].type == List[ListCount].type &&
                    Cuts[ineq].i == List[ListCount].i &&
                    Cuts[ineq].j == List[ListCount].j &&
                    Cuts[ineq].k == List[ListCount].k) 
            {
                found_ineq = 1;
            }
        }

        // If inequality not already in Cuts, add it to Cuts
        if (!found_ineq) {

            Cuts[next_ineq].type  = List[ListCount].type;
            Cuts[next_ineq].i     = List[ListCount].i;
            Cuts[next_ineq].j     = List[ListCount].j;
            Cuts[next_ineq].k     = List[ListCount].k;
            Cuts[next_ineq].value = List[ListCount].value;
            Cuts[next_ineq].y     = 0.0;     // set dual multipliers of new triangle ineq to 0

            ++next_ineq;
            ++added;
        }
    }
    
    PP->NIneq += added;

    *NumAdded = added;
    *NumSubtracted = subtracted;

    return maxAllIneq;
}


/************************* PENTAGONAL INEQUALITIES *************************/

/* 
 * Separates pentagonal inequalities using the X matrix and simulated annealing heuristic for 
 * QAP and fills the Pentagonal_Inequality array named Pent_List, with at most params.PentIneq
 * inequalities. It also returns the value of the inequlity that is violated the most by X.
 */
double getViolated_PentagonalInequalities(double *X, int N, Pentagonal_Inequality *Pent_List, int *ListSize) {

    int ListCount;                                  // loop index
    int size = 0;                                   // number of added cuts
    int LeastViolatedIneq = 0;                      // index of least violated inequality
    double LeastViolatedIneqValue = -BIG_NUMBER;    // maximum violation
    double minAllIneq = BIG_NUMBER;                 // minimum violation
    double test_ineqvalue;                          // violation of current cut


    // 5 tuple of indeces defining the violated pentagonal inequality
    static int pent[5];

    /* 5x5 matrices that define pentagonal inequalities are stored as rows in H */

    // H1 = ee^T, where e is vector of all ones
    // e[0] = -1; H2 = ee^T
    // e[0] = -1, e[1] = -1; H3 = ee^T
    static int H[3][25] = { {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},                    // H1
                            {1, -1, -1, -1, -1, -1, 1, 1, 1, 1, -1, 1, 1, 1, 1, -1, 1, 1, 1, 1, -1, 1, 1, 1, 1},            // H2
                            {1, 1, -1, -1, -1, 1, 1, -1, -1, -1, -1, -1, 1, 1, 1, -1, -1, 1, 1, 1, -1, -1, 1, 1, 1} };      // H3


    for (int num_trial = 0; num_trial < params.Pent_Trials; ++num_trial) {
        for (int type = 1; type <= 3; ++type) {

            test_ineqvalue = qap_simulated_annealing(&H[type-1][0], 5, X, N, pent);

            // keep track of the minimum value of test_ineqvalue, i.e. 
            // current most violated cut value
            minAllIneq = (test_ineqvalue < minAllIneq) ? test_ineqvalue : minAllIneq;

            if (test_ineqvalue < 1.0) {

                // (1) put first params.PentIneq violated inequalities in list and 
                //     keep track of the least violated inequality
                if (size < params.PentIneq) {

                    // add ineq to the end of Pent_List
                    Pent_List[size].type  = type;
                    Pent_List[size].value = test_ineqvalue;
                    for (int i = 0; i < 5; ++i)
                        Pent_List[size].permutation[i] = pent[i];
                    
                    // update LeastViolatedIneq
                    if (size == 0 || test_ineqvalue > LeastViolatedIneqValue) {
                        LeastViolatedIneqValue = test_ineqvalue;
                        LeastViolatedIneq = size;
                    }

                    ++size;
                } 
                else if (test_ineqvalue < LeastViolatedIneqValue) 
                {
                    // (2) if you find an inequality that is violated more than the 
                    //     least violated inequality, then add it to the list, 
                    //     remove the least violated inequality, and find 
                    //     the new least violated inequality in the list

                    // add ineq to the list, replacing LeastViolatedIneq
                    Pent_List[LeastViolatedIneq].type  = type;
                    Pent_List[LeastViolatedIneq].value = test_ineqvalue;
                    for (int i = 0; i < 5; ++i)
                        Pent_List[LeastViolatedIneq].permutation[i] = pent[i];                    

                    // update LeastViolatedIneq
                    LeastViolatedIneqValue = BIG_NUMBER;
                    for (ListCount = 0; ListCount < size; ++ListCount) 
                    {
                        // if ineq is less violated
                        if (Pent_List[ListCount].value > LeastViolatedIneqValue) 
                        {
                            LeastViolatedIneqValue = Pent_List[ListCount].value;
                            LeastViolatedIneq = ListCount;
                        }
                    }
                }
            }

        }        
    }

    *ListSize = size;

    return minAllIneq;
}



/* update inequalities: purge old ones and separate new ones */
/* returns the new maximum violation of pentagonal inequalities */

double updatePentagonalInequalities(Problem *PP, double *y, int *NumAdded, int *NumSubtracted, int triag) {

    int ineq;                   // index for inequality        
    int yindex;                 // index for dual multipliers          
    int ListCount, ListSize;
    int N = PP->n;



    // purge inequalities: decide which cuts that were previously added to remove
    int subtracted = 0;
    yindex = triag;             // first triag elements of dual vector y belong to triangle inequalities!
    int next_ineq = 0;

    for (ineq = 0; ineq < PP->NPentIneq; ++ineq) {

        // store the dual multiplier
        Pent_Cuts[ineq].y = y[yindex];
        ++yindex;

        // remove inequality if dual multiplier is small
        if (Pent_Cuts[ineq].y < 1e-5) {
            ++subtracted;
        } 
        else // keep inequality
        {
            Pent_Cuts[next_ineq].type  = Pent_Cuts[ineq].type;
            Pent_Cuts[next_ineq].y     = Pent_Cuts[ineq].y;

            for (int i = 0; i < 5; ++i)
                Pent_Cuts[next_ineq].permutation[i] = Pent_Cuts[ineq].permutation[i];

            ++next_ineq;
        }
    }

    PP->NPentIneq -= subtracted;

    // separate new pentagonal inequalities
    double maxAllIneq = getViolated_PentagonalInequalities(X, N, Pent_List, &ListSize);



    // Add List to Cuts
    int added = 0;
    for (ListCount = 0; ListCount < ListSize; ++ListCount) { 

        // Stop if we have reached the maximum number of cuts we can add
        if (next_ineq == MaxPentIneqAdded)
            break;

        // Check if inequality is already included in Pent_Cuts
        int found_ineq = 0;
        for (ineq = 0; ineq < PP->NPentIneq; ++ineq) {

            if (Pent_Cuts[ineq].type == Pent_List[ListCount].type &&
                Pent_Cuts[ineq].permutation[0] == Pent_List[ListCount].permutation[0] &&
                Pent_Cuts[ineq].permutation[1] == Pent_List[ListCount].permutation[1] &&
                Pent_Cuts[ineq].permutation[2] == Pent_List[ListCount].permutation[2] &&
                Pent_Cuts[ineq].permutation[3] == Pent_List[ListCount].permutation[3] &&
                Pent_Cuts[ineq].permutation[4] == Pent_List[ListCount].permutation[4])
            {
                found_ineq = 1;
            }
        }

        // If inequality not already in Pent_Cuts, add it to Pent_Cuts
        if (!found_ineq) {

            Pent_Cuts[next_ineq].type           = Pent_List[ListCount].type;
            Pent_Cuts[next_ineq].permutation[0] = Pent_List[ListCount].permutation[0];
            Pent_Cuts[next_ineq].permutation[1] = Pent_List[ListCount].permutation[1];
            Pent_Cuts[next_ineq].permutation[2] = Pent_List[ListCount].permutation[2];
            Pent_Cuts[next_ineq].permutation[3] = Pent_List[ListCount].permutation[3];
            Pent_Cuts[next_ineq].permutation[4] = Pent_List[ListCount].permutation[4];

            Pent_Cuts[next_ineq].value          = Pent_List[ListCount].value;

            // set dual multipliers of new pentagonal ineq to 0
            Pent_Cuts[next_ineq].y              = 0.0;     

            ++next_ineq;
            ++added;
        }
    }
    
    PP->NPentIneq += added;

    *NumAdded = added;
    *NumSubtracted = subtracted;

    return maxAllIneq;

}


/************************* HEPTAGONAL INEQUALITIES *************************/

/* 
 * Separates heptagonal inequalities using the X matrix and simulated annealing heuristic for 
 * QAP and fills the Heptagonal_Inequality array named Hepta_List, with at most params.HeptIneq
 * inequalities. It also returns the value of the inequlity that is violated the most by X.
 */
double getViolated_HeptagonalInequalities(double *X, int N, Heptagonal_Inequality *Hepta_List, int *ListSize) {

    int ListCount;                                  // loop index
    int size = 0;                                   // number of added cuts
    int LeastViolatedIneq = 0;                      // index of least violated inequality
    double LeastViolatedIneqValue = -BIG_NUMBER;    // maximum violation
    double minAllIneq = BIG_NUMBER;                 // minimum violation
    double test_ineqvalue;                          // violation of current cut


    // 7 tuple of indeces defining the violated heptagonal inequality
    static int hept[7];

    /* 7x7 matrices that define heptagonal inequalities are stored as rows in H */

    // H1 = ee^T, where e is vector of all ones
    // e[0] = -1; H2 = ee^T
    // e[0] = -1, e[1] = -1; H3 = ee^T
    // e[0] = -1, e[1] = -1, e[2] = -1; H4 = ee^T
    static int H[4][49] = { {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},                              // H1
                            {1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1},                  // H2
                            {1, 1, -1, -1, -1, -1, -1, 1, 1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, -1, -1, 1, 1, 1, 1, 1, -1, -1, 1, 1, 1, 1, 1, -1, -1, 1, 1, 1, 1, 1, -1, -1, 1, 1, 1, 1, 1},          // H3 
                            {1, 1, 1, -1, -1, -1, -1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, 1, 1, 1, 1} };    // H4


    for (int num_trial = 0; num_trial < params.Hepta_Trials; ++num_trial) {
        for (int type = 1; type <= 4; ++type) {

            test_ineqvalue = qap_simulated_annealing(&H[type-1][0], 7, X, N, hept);

            // keep track of the minimum value of test_ineqvalue, i.e. 
            // current most violated cut value
            minAllIneq = (test_ineqvalue < minAllIneq) ? test_ineqvalue : minAllIneq;

            if (test_ineqvalue < 1.0) {

                // (1) put first params.PentIneq violated inequalities in list and 
                //     keep track of the least violated inequality
                if (size < params.HeptaIneq) {

                    // add ineq to the end of Pent_List
                    Hepta_List[size].type  = type;
                    Hepta_List[size].value = test_ineqvalue;
                    for (int i = 0; i < 7; ++i)
                        Hepta_List[size].permutation[i] = hept[i];
                    
                    // update LeastViolatedIneq
                    if (size == 0 || test_ineqvalue > LeastViolatedIneqValue) {
                        LeastViolatedIneqValue = test_ineqvalue;
                        LeastViolatedIneq = size;
                    }

                    ++size;
                } 
                else if (test_ineqvalue < LeastViolatedIneqValue) 
                {
                    // (2) if you find an inequality that is violated more than the 
                    //     least violated inequality, then add it to the list, 
                    //     remove the least violated inequality, and find 
                    //     the new least violated inequality in the list

                    // add ineq to the list, replacing LeastViolatedIneq
                    Hepta_List[LeastViolatedIneq].type  = type;
                    Hepta_List[LeastViolatedIneq].value = test_ineqvalue;
                    for (int i = 0; i < 7; ++i)
                        Hepta_List[LeastViolatedIneq].permutation[i] = hept[i];                    

                    // update LeastViolatedIneq
                    LeastViolatedIneqValue = BIG_NUMBER;
                    for (ListCount = 0; ListCount < size; ++ListCount) 
                    {
                        // if ineq is less violated
                        if (Hepta_List[ListCount].value > LeastViolatedIneqValue) 
                        {
                            LeastViolatedIneqValue = Hepta_List[ListCount].value;
                            LeastViolatedIneq = ListCount;
                        }
                    }
                }
            }

        }        
    }

    *ListSize = size;

    return minAllIneq;
}


/* update inequalities: purge old ones and separate new ones */
/* returns the new maximum violation of heptagonal inequalities */

double updateHeptagonalInequalities(Problem *PP, double *y, int *NumAdded, int *NumSubtracted, int hept_index) {

    int ineq;                   // index for inequality        
    int yindex;                 // index for dual multipliers          
    int ListCount, ListSize;
    int N = PP->n;

    // purge inequalities: decide which cuts that were previously added to remove
    int subtracted = 0;
    yindex = hept_index;        // first hept_index elements of dual vector y belong to triangle and pentagonal inequalities!
    int next_ineq = 0;

    for (ineq = 0; ineq < PP->NHeptaIneq; ++ineq) {

        // store the dual multiplier
        Hepta_Cuts[ineq].y = y[yindex];
        ++yindex;

        // remove inequality if dual multiplier is small
        if (Hepta_Cuts[ineq].y < 1e-5) {
            ++subtracted;
        } 
        else // keep inequality
        {
            Hepta_Cuts[next_ineq].type  = Hepta_Cuts[ineq].type;
            Hepta_Cuts[next_ineq].y     = Hepta_Cuts[ineq].y;

            for (int i = 0; i < 7; ++i)
                Hepta_Cuts[next_ineq].permutation[i] = Hepta_Cuts[ineq].permutation[i];

            ++next_ineq;
        }
    }

    PP->NHeptaIneq -= subtracted;

    // separate new heptagonal inequalities
    double maxAllIneq = getViolated_HeptagonalInequalities(X, N, Hepta_List, &ListSize);

    // Add List to Cuts
    int added = 0;
    for (ListCount = 0; ListCount < ListSize; ++ListCount) { 

        // Stop if we have reached the maximum number of cuts we can add
        if (next_ineq == MaxHeptaIneqAdded)
            break;

        // Check if inequality is already included in Hepta_Cuts
        int found_ineq = 0;
        for (ineq = 0; ineq < PP->NHeptaIneq; ++ineq) {

            if (Hepta_Cuts[ineq].type == Hepta_List[ListCount].type &&
                Hepta_Cuts[ineq].permutation[0] == Hepta_List[ListCount].permutation[0] &&
                Hepta_Cuts[ineq].permutation[1] == Hepta_List[ListCount].permutation[1] &&
                Hepta_Cuts[ineq].permutation[2] == Hepta_List[ListCount].permutation[2] &&
                Hepta_Cuts[ineq].permutation[3] == Hepta_List[ListCount].permutation[3] &&
                Hepta_Cuts[ineq].permutation[4] == Hepta_List[ListCount].permutation[4] &&
                Hepta_Cuts[ineq].permutation[5] == Hepta_List[ListCount].permutation[5] &&
                Hepta_Cuts[ineq].permutation[6] == Hepta_List[ListCount].permutation[6])
            {
                found_ineq = 1;
            }
        }

        // If inequality not already in Hepta_Cuts, add it to Hepta_Cuts
        if (!found_ineq) {

            Hepta_Cuts[next_ineq].type           = Hepta_List[ListCount].type;
            Hepta_Cuts[next_ineq].permutation[0] = Hepta_List[ListCount].permutation[0];
            Hepta_Cuts[next_ineq].permutation[1] = Hepta_List[ListCount].permutation[1];
            Hepta_Cuts[next_ineq].permutation[2] = Hepta_List[ListCount].permutation[2];
            Hepta_Cuts[next_ineq].permutation[3] = Hepta_List[ListCount].permutation[3];
            Hepta_Cuts[next_ineq].permutation[4] = Hepta_List[ListCount].permutation[4];
            Hepta_Cuts[next_ineq].permutation[5] = Hepta_List[ListCount].permutation[5];
            Hepta_Cuts[next_ineq].permutation[6] = Hepta_List[ListCount].permutation[6];

            Hepta_Cuts[next_ineq].value          = Hepta_List[ListCount].value;

            // set dual multipliers of new heptagonal ineq to 0
            Hepta_Cuts[next_ineq].y              = 0.0;     

            ++next_ineq;
            ++added;
        }
    }
    
    PP->NHeptaIneq += added;

    *NumAdded = added;
    *NumSubtracted = subtracted;

    return maxAllIneq;

}
