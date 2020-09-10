#include "admm.h"

extern Parameters params;
extern FILE *output;
extern int BabPbSize;      

/*
 * Evaluate a specific node.
 * This function computes the upper and lower bounds of a specific node
 * (calls SDP bound function) and returns the upper bound of the node
 */
double Evaluate(BabNode *node, Problem *SP, Problem *PP, int rank) {  

    // create subproblem PP
    createSubproblem(node, SP, PP);

    // compute the SDP upper bound and run heuristic
    double bound = SDPbound(node, SP, PP, rank);    

    return bound;
}



/*
 * Writes subproblem to PP.
 *
 * Computes the subproblem removing the rows and the columns of the 
 * fixed variables upper left corner of SP->L
 * 
 * SP is the original problem
 * PP is the subproblem (some variables are fixed)
 *
 * PP is made from SP
 * Function prepares objective matrix L for model in -1,1 variables:
 * 
 * max x'LX, s.t. x in {-1,1}^(PP->n)
 */
void createSubproblem(BabNode *node, Problem *SP, Problem *PP) {

    // Subproblem size is the number of non-fixed variables in the node
    PP->n = BabPbSize + 1 - countFixedVariables(node);

    /* build objective:
     * Laplacian;
     * z'*L*z = sum_{i != fixed, j != fixed} L_ij*xi*xj (smaller matrix L_bar for subproblem)
              + sum_{i = fixed, j = fixed} L_ij*x_i*xj (getFixedValue)
              + sum_{rows of fixed vertices without fixed entries}  (linear part, that is twice added to diagonal of L_bar)
     */

    /* Laplacian is created by deleting appropriate rows and cols
     * of upper left corner of SP->L 
     */ 
    int index = 0;
    int N = SP->n;
    double row_sum = 0.0;
    
    // rows which are deleted due to fixed variable
    // later add to diagonal
    double fixedRow[PP->n - 1];
    for (int i = 0; i < PP->n - 1; ++i)
        fixedRow[i] = 0.0;

    // counter for fixedRow
    int fixed = 0;      
       
    // last element (lower right corner) is sum
    double sum = 0.0;   


    for (int i = 0; i < BabPbSize; ++i) {
        for (int j = 0; j < BabPbSize; ++j) {
            if (!node->xfixed[i] && !node->xfixed[j]) {     // delete rows and cols of SP->L
                PP->L[index] = SP->L[j + i*N];
                row_sum += PP->L[index];
                ++index;               
            }
            else if ( (node->xfixed[i] && node->sol.X[i] == 1) && !node->xfixed[j]) { // save fixed rows to add to diagonal
                fixedRow[fixed] += SP->L[j + i*N];
                ++fixed;
            }
        }

        if (!node->xfixed[i]) {
            PP->L[index] = row_sum;       // vector part of PP->L (last column)
            ++index;
        }    

        // row scaned, set to 0
        row_sum = 0.0;
        fixed = 0;
    }

    // add last row (copy from last column)
    for (int i = 0; i < PP->n - 1; ++i)
        PP->L[i + (PP->n - 1) * PP->n] = PP->L[PP->n - 1 + i * PP->n];
    
   
    /* LINEAR PART OF PP->L:   add 2x vector fixedRow to diagonal, last col and last row */
    for (int i = 0; i < PP->n - 1; ++i) {
        PP->L[i + i * PP->n] += 2*fixedRow[i];
        PP->L[PP->n - 1 + i * PP->n] += 2*fixedRow[i];
        PP->L[i + (PP->n - 1) * PP->n] += 2*fixedRow[i];

        sum += PP->L[i + (PP->n - 1) * PP->n];
    }

    /* CONSTANT PART OF PP->L:  element (PP->n - 1, PP->n - 1) */
    PP->L[PP->n - 1 + (PP->n - 1) * PP->n] = sum;
    
    /* multiple by 1/4 the whole matrix L */
    double alpha = 0.25;
    int inc = 1;
    int nn = (PP->n) * (PP->n);
    dscal_(&nn,&alpha,PP->L,&inc);
}


/* 
 * Return the fixed value of the node.
 * The fixed value is contribution of the fixed variables to 
 * the objective value.
 */
double getFixedValue(BabNode *node, Problem *SP) {

    int N = SP->n;
    double fixedvalue = 0.0;

    for (int i = 0; i < BabPbSize; ++i) {
        for (int j = 0; j < BabPbSize; ++j) {
            if (node->xfixed[i] && node->xfixed[j]) {
                fixedvalue += SP->L[j + i*N] * node->sol.X[i] * node->sol.X[j];
            }
        }
    }

    return fixedvalue;

}
