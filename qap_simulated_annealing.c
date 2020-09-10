#include <math.h>

#include "admm.h"

/* simulated annealing heuristic to find approximate solution of QAP:
 *      
 *      min <(PXP^T)(1:k,1:k), H> over permutation matrices P
 *
 * H is kxk matrix that determines k-gonal inequalities (pentagonal or heptagonal)
 * ineq stores the most violated k-gonal inequality
 */ 
double qap_simulated_annealing(int *H, int k, double *X, int n, int *ineq) {

    /* copy upper triangle to lower */
    for (int row = 0; row < n; ++row)     
        for (int col = 0; col < row; ++col)    
            X[col + n * row] = X[row + n * col];

    /* parameters */
    int inner_iter = n;
    double reduce_temp = 0.6;
    double increase_iter = 1.1;

    /* initial temperature t = sum(sum(abs(H)) * sum(sum(abs(X))) / (n * (n-1)) */
    double t = k * k;  // added sum(sum(abs(H)) = k^2

    //add sum of diagonal elemets of X, which equals n since diag(X) = e
    double sumX = 0;

    // add upper triangle of X 
    for (int i = 0; i < n; ++i) {
        sumX += X[i + i*n];
        for (int j = i + 1; j < n; ++j) 
            sumX += 2 * fabs(X[j + i*n]);
    }
    
    t *= sumX;
    t /= (n * (n-1));


    // best found cost
    double cost = BIG_NUMBER;

    // permutation of n elements: current and best
    int *perm;
    int *perm_best;
    alloc_vector(perm, n, int);
    alloc_vector(perm_best, n, int);

    // other variables
    double t1;          // current temperature
    int m1;             // current number of inner iterations
    double sol;         // current solution
    double sol_temp;    // temporary solution when applying transposition to permutation
    double delta;       // difference between sol_temp and sol

    double dt1;         // delta/t1;
    double prob;        // probability for acepting worse soution

    int not_done;
    int accept;         // indicator whether to accept new permutation

    // temporary variables
    int random_num, temp, i1, i2;
    double tempValue1, tempValue2;


    /* generate random permutation */
    // a) perm = 0:n-1
    for (int i = 0; i < n; ++i)
        perm[i] = i;

    // b) from end till start do random transpositions of elements
    for (int i = n - 1; i > 0; --i) {
        random_num = (rand() % i);
        // swap entries in perm at position random_num i
        temp = perm[i];
        perm[i] = perm[random_num];
        perm[random_num] = temp;
    }

    /* initialize current values */
    t1 = t;
    m1 = inner_iter;
    

    /* current solution = <H, X> */
    // add diagonal terms (both matrices have ones on diagonal)
    sol = 0;

    // add off-diagonal terms
    for (int i = 0; i < k; ++i) {
        sol += X[perm[i] + perm[i] * n];
        for (int j = i + 1; j < k; ++j) {
            sol += 2 * H[j + i*k] * X[perm[j] + perm[i] * n];
        }
    }




    not_done = 1;            
    /* while things change */
    while (not_done > 0) {

        not_done = 0;

        /* m1 iterations at constant temperature */
        for (int num_it = 0; num_it < m1; ++num_it) {


            /* generate two random numbers:
             * i1 in (0,...,k-1) and i2 (0,...,n-1) 
             * with i1 <= i2
             */
            i1 = rand() % k;
            i2 = rand() % n;
            if (i2 < i1) {
                temp = i1;
                i1 = i2;
                i2 = temp;
            }

            // compute cost after transposing perm(i1) and perm(i2)
            delta = 0.0;
            for (int i = 0; i < k; ++i)
                delta += H[i + k*i1] * (X[perm[i] + n * perm[i2]] - X[perm[i] + n * perm[i1]]);


            if (i2 < k) {
                for (int i = 0; i < k; ++i)
                    delta += H[i + k*i2] * (X[perm[i] + n * perm[i1]] - X[perm[i] + n * perm[i2]]);
            }

            delta *= 2.0;

            tempValue1 = X[perm[i1] + n * perm[i1]] + X[perm[i2] + n * perm[i2]] - 2.0 * X[perm[i2] + n * perm[i1]];

            tempValue2 = 1.0;
            if (i2 < k)
                tempValue2 += -2 * H[i1 + k*i2] + 1.0;
                
            delta += tempValue1 * tempValue2;    

            sol_temp = sol + delta;

            // determine whether the swap is accepted
            if (delta > 0) { // if solution gets worse accept with probability

                dt1 = delta/t1;
                if (dt1 > 5)
                    accept = 0;
                else {
                    prob = exp(-dt1);
                    if ( ((double)rand()/((double)RAND_MAX)) < prob )
                        accept = 1;
                    else 
                        accept = 0;
                }

            }
            else // if solution is better accept
                accept = 1;


            /* do the update if the swap is accepted */
            if (accept) {

                if (fabs(delta) > 0.0001)
                    not_done = 1;

                // make transposition in perm
                temp = perm[i1];
                perm[i1] = perm[i2];
                perm[i2] = temp;

                sol = sol_temp;

                /* if the just found solution is the best one of all, store it */
                if (sol < cost) {
                    cost = sol;
                    for (int i = 0; i < n; ++i)
                        perm_best[i] = perm[i];
                }

            }


        } /* End m1 iterations at constant temperature */

        // reduce temperature
        t1 *= reduce_temp;
        m1 = round(m1 * increase_iter);

    }

    // save final k-gonal inequality
    for (int i = 0; i < k; ++i)
        ineq[i] = perm_best[i];


    free(perm);
    free(perm_best);

    return cost;
}
