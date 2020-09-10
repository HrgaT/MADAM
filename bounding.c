#include <math.h>

#include "admm.h"

extern Parameters params;
extern FILE *output;
extern int BabPbSize;

extern double TIME;                 

extern Triangle_Inequality *Cuts;               // vector of triangle inequality constraints
extern Pentagonal_Inequality *Pent_Cuts;        // vector of pentagonal inequality constraints
extern Heptagonal_Inequality *Hepta_Cuts;       // vector of heptagonal inequality constraints

extern double *X;                   // current X
extern double *y;                   // dual multiplier for diagonal constraints
extern double *Q;                   // Z matrix in admm method
extern double *s;                   // nonnegative dual multiplier to constraint u >= 0
                                    // used for purging of cuts
extern double f;                    // function value of relaxation

extern double diff;		            // difference between basic SDP relaxation and bound with added cutting planes

/******** main bounding routine calling ADMM method ********/
double SDPbound(BabNode *node, Problem *SP, Problem *PP, int rank) {

    int index;                      // helps to store the fractional solution in the node
    double bound;                   // f + fixedvalue
    double gap;                     // difference between best lower bound and upper bound
    double oldf;                    // stores f from previous iteration 
    int x[BabPbSize];               // vector for heuristic
    double viol3 = 0.0;             // maximum violation of triangle inequalities
    double viol5 = 0.0;             // maximum violation of pentagonal inequalities
    double viol7 = 0.0;             // maximum violation of heptagonal inequalities
    int count = 0;                  // number of iterations (adding and purging of cutting planes)
    int nbit;                       // number of iterations in ADMM method

    int triag;                      // starting index for pentagonal inequalities in vectors t and s 
    int penta;                      // starting index for heptagonal inequalities in vectors t and s

    /* stopping conditions */
    int done = 0;                   
    int giveup = 0;                                   
    int prune = 0;

    // Parameters
    double sigma = params.sigma0;
    double tol = params.tol0;

    // fixed value contributes to the objective value
    double fixedvalue = getFixedValue(node, SP);

    /*** start with no cuts ***/
    // triangle inequalities
    PP->NIneq = 0; 
    int Tri_NumAdded = 0;
    int Tri_NumSubtracted = 0;

    // pentagonal inequalities
    PP->NPentIneq = 0;
    int Pent_NumAdded = 0;
    int Pent_NumSubtracted = 0;

    // heptagonal inequalities
    PP->NHeptaIneq = 0;
    int Hepta_NumAdded = 0;
    int Hepta_NumSubtracted = 0;         

    /* solve basic SDP relaxation with interior-point method */
    ipm_mc_pk(PP->L, PP->n, X, y, Q, &f, 0);

    // store basic SDP bound to compute diff in the root node
    double basic_bound = f + fixedvalue;

    // Store the fractional solution in the node    
    index = 0;
    for (int i = 0; i < BabPbSize; ++i) {
        if (node->xfixed[i]) {
            node->fracsol[i] = (double) node->sol.X[i];
        }
        else {
            // convert x (last column X) from {-1,1} to {0,1}
            node->fracsol[i] = 0.5*(X[(PP->n - 1) + index*PP->n] + 1.0); 
            ++index;
        }
    }

    /* run heuristic */
    for (int i = 0; i < BabPbSize; ++i) {
        if (node->xfixed[i]) {
            x[i] = node->sol.X[i];
        }
        else {
            x[i] = 0;
        }
    }

    runHeuristic(SP, PP, node, x);
    updateSolution(x);

    // upper bound
    bound = f + fixedvalue;

    // check pruning condition
    if ( bound < Bab_LBGet() + 1.0 ) {
        prune = 1;
        goto END;
    }

    // check if cutting planes need to be added     
    if (params.use_diff && (rank != 0) && (bound > Bab_LBGet() + diff + 1.0)) {
        giveup = 1;
        goto END;
    }

    /* separate first triangle inequality */
    viol3 = updateTriangleInequalities(PP, s, &Tri_NumAdded, &Tri_NumSubtracted);

   
    /*** Main loop ***/
    while (!done) {

        // Update iteration counter
        ++count;
        oldf = f;

        // Call ADMM solver
        ADMM_solver(PP, X, y, Q, &sigma, tol, &nbit, 0);

        // increase precision
        tol *= params.scaleTol;

        if (tol < params.minTol)
            tol = params.minTol;

        // upper bound
        bound = f + fixedvalue;

        // prune test
        prune = ( bound < Bab_LBGet() + 1.0 ) ? 1 : 0;
 
        /******** heuristic ********/
        if (!prune) {

            for (int i = 0; i < BabPbSize; ++i) {
                if (node->xfixed[i]) {
                    x[i] = node->sol.X[i];
                }
                else {
                    x[i] = 0;
                }
            }

            runHeuristic(SP, PP, node, x);
            updateSolution(x);

            prune = ( bound < Bab_LBGet() + 1.0 ) ? 1 : 0;
        }
        /***************************/

        // compute gap
        gap = bound - Bab_LBGet();

        /* check if we will not be able to prune the node */
        if (count == params.triag_iter + params.pent_iter + params.hept_iter) {
            if ( (gap - 1.0 > (oldf - f)*(params.max_outer_iter - count)))
                giveup = 1;
        }

        /* check if extra iterations can close the gap */
        if (count == params.max_outer_iter) {
            if ( gap - 1.0 > (oldf - f)*params.extra_iter )
                giveup = 1;
        }
        
        /* max number of iterations reached */
        if (count == params.max_outer_iter + params.extra_iter)
            giveup = 1; 


        /* increase number of pentagonal and heptagonal inequalities (during separation)
         * if the gap is still too big 
         */
        if ((rank != 0) && giveup && !prune) {
            params.Pent_Trials += 60;    // add 3 types * 60 = 180 pentagonal inequalities
            params.Hepta_Trials += 50;  // add 4 types * 50 = 200 heptagonal inequalities
        }

        // purge inactive cutting planes, add new inequalities
        if (!prune && !giveup) {
            
            triag = PP->NIneq;          // save number of triangle and pentagonal inequalities before purging
            penta = PP->NPentIneq;      // --> to know with which index in dual vector gamma, pentagonal
                                        // and heptagonal inequalities start!

            viol3 = updateTriangleInequalities(PP, s, &Tri_NumAdded, &Tri_NumSubtracted);
                      
            /* include pentagonal and heptagonal inequalities */          
            if ( params.include_Pent && (count >= params.triag_iter || viol3 < 0.2) )
                viol5 = updatePentagonalInequalities(PP, s, &Pent_NumAdded, &Pent_NumSubtracted, triag);  

            if ( params.include_Hepta && ( (count >= params.triag_iter + params.pent_iter) || (viol3 < 0.2 && (1 - viol5 < 0.4)) ) )
                viol7 = updateHeptagonalInequalities(PP, s, &Hepta_NumAdded, &Hepta_NumSubtracted, triag + penta);      
        }
        else {               
            Tri_NumAdded = 0;
            Tri_NumSubtracted = 0;
            Pent_NumAdded = 0;
            Pent_NumSubtracted = 0;
            Hepta_NumAdded = 0;
            Hepta_NumSubtracted = 0;
        }

        


        // Test stopping conditions
        done = 
            prune ||                         // can prune the B&B tree 
            giveup ||                        // upper bound to far away from lower bound
            (nbit >= params.ADMM_itermax);   // ADMM reached max iter

        // Store the fractional solution in the node    
        index = 0;
        for (int i = 0; i < BabPbSize; ++i) {
            if (node->xfixed[i]) {
                node->fracsol[i] = (double) node->sol.X[i];
            }
            else {
                // convert x (last column X) from {-1,1} to {0,1}
                node->fracsol[i] = 0.5*(X[(PP->n - 1) + index*PP->n] + 1.0); 
                ++index;
            }
        }

    } // end while loop

    bound = f + fixedvalue;

    // compute difference between basic SDP relaxation and bound with added cutting planes
    if (rank == 0) {
        diff = basic_bound - bound;
    }

    END:
      

    return bound;

}

