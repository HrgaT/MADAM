#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#include "admm.h"
  
  /* defined in heap.c */  
extern Heap *heap;  
extern int BabPbSize;           
extern BabSolution *BabSol;
   
extern Parameters params;
extern Problem *SP;
extern Problem *PP;

extern double root_bound;
extern double TIME;
extern int stopped;

extern int num_workers_used;


/* initialize global lower bound to 0 and global solution vector to zero */
void initializeBabSolution() { 

    BabSolution bs;

    for (int i = 0; i < BabPbSize; ++i) {
        bs.X[i] = 0;
    }

    Bab_LBInit(0, &bs);
}


/************** Initialization: root node and priority queue **************/
int Init_PQ(void) {
    
    int over = 0;
    extern BabNode *BabRoot;

    // Create the root node
    BabRoot = newNode(NULL);

    // increase number of evaluated nodes
    Bab_incEvalNodes();

    // Evaluate root node: compute upper and lower bound 
    root_bound = Evaluate(BabRoot, SP, PP, 0);
    printf("Root node bound: %.2f\n", root_bound);

    // save upper bound
    BabRoot->upper_bound = root_bound;

    /* insert node into the priority queue or prune */
    // NOTE: optimal solution has INTEGER value, i.e. add +1 to lower bound
    if (Bab_LBGet() + 1.0 < BabRoot->upper_bound) {    
        Bab_PQInsert(BabRoot); 
    }
    else {
        // otherwise, intbound <= BabLB, so we can prune
        over = -1;
        free(BabRoot);
    }

    return over;
}


/* Bab function which initializes the problem and allocates the structures */
int Bab_Init(int argc, char **argv, int rank) {

    int read_error = 0;

    // Process the command line arguments
    if ( (read_error = processCommandLineArguments(argc, argv, rank)) )
        return read_error;

    // Seed the random number generator
    srand(2020);

    // Provide B&B with an initial solution
    initializeBabSolution();

    // Allocate the memory
    allocMemory();

    return read_error;
}

/* NOTE: int *sol in functions evaluateSolution and updateSolution have length BabPbSize
 * -> to get objecive multiple with Laplacian that is stored in upper left corner of SP->L
 */
double evaluateSolution(int *sol) {

    double val = 0.0;
    
    for (int i = 0; i < BabPbSize; ++i) {
        for (int j = 0; j < BabPbSize; ++j) {
            val += SP->L[j + i * SP->n] * sol[i] * sol[j];
        }
    }

    return val;
}


/*
 * Only this function can update best solution and value.
 * Returns 1 if success.
 */
int updateSolution(int *x) {
    
    int solutionAdded = 0;
    double sol_value;
    BabSolution solx;

    // Copy x into solx --> because Bab_LBUpd needs BabSolution and not int*
    for (int i = 0; i < BabPbSize; ++i) {
      solx.X[i] = x[i];
    }

    sol_value = evaluateSolution(x); // computes objective value of solx

    /* If new solution is better than the global solution, 
     * then update and print the new solution. */
    if (Bab_LBUpd(sol_value, &solx)) {
        solutionAdded = 1;
    }
    
    return solutionAdded;
}


/* MASTER process main routine */
void master_Bab_Main(Message message, int source, int *busyWorkers, int numbWorkers, int *numbFreeWorkers, MPI_Datatype BabSolutiontype) {

#if 0
    // If the algorithm stops before finding the optimal solution
    if (!stopped && (params.time_limit > 0 && (MPI_Wtime() - TIME) > params.time_limit) ) {
        
        // signal to printFinalOutput that algorihtm stopped early
        stopped = 1;        
    }
#endif

    MPI_Status status;

    switch(message) {

        case IDLE:

            busyWorkers[source] = 0;
            ++(*numbFreeWorkers);
            break;  

        case NEW_VALUE:
        {
        
            // receive best lower bound and corresponding feasible solution
            double g_lowerBound;
            BabSolution solx;
            
            MPI_Recv(&g_lowerBound, 1, MPI_DOUBLE, source, LOWER_BOUND, MPI_COMM_WORLD, &status);
            MPI_Recv(&solx, 1, BabSolutiontype, source, SOLUTION, MPI_COMM_WORLD, &status);  

            if ( Bab_LBUpd(g_lowerBound, &solx) ){
                printf("Feasible solution %.0lf\n", Bab_LBGet());
            }
            
            // send update information back to worker
            g_lowerBound = Bab_LBGet();

            MPI_Send(&g_lowerBound, 1, MPI_DOUBLE, source, LOWER_BOUND, MPI_COMM_WORLD);
            break;       
        }
        case SEND_FREEWORKERS:
        {            
            // get number of requested workers            
            int workers_request;                
            MPI_Recv(&workers_request, 1, MPI_INT, source, FREEWORKER, MPI_COMM_WORLD, &status);
                        
            // compute number of freeworkers
            int num_workers_available = (workers_request < *numbFreeWorkers) ? workers_request : *numbFreeWorkers;
            int available_workers[num_workers_available];

            for(int i = 1, j = 0; (i < numbWorkers) && (j < num_workers_available); ++i)    // master has rank 0 and is not considered
            {
                if(busyWorkers[i] == 0){ // is free
                    available_workers[j] = i;
                    ++j;
                    busyWorkers[i] = 1; // set to busy
                    --(*numbFreeWorkers);
                }
            }
    

            // worker branched subproblem in local queue --> add 2 bab nodes
            Bab_incEvalNodes();
            Bab_incEvalNodes(); 

	        // count current number of busy workers
	        int current_busy = 0;

	        for (int i = 1; i < numbWorkers; ++i) {
		       if (busyWorkers[i] == 1)
		       ++current_busy;
	        } 
	
	        num_workers_used = (current_busy > num_workers_used) ? current_busy : num_workers_used;

	        // send message back
            double g_lowerBound = Bab_LBGet();            
            MPI_Send(&num_workers_available, 1, MPI_INT, source, NUM_FREE_WORKERS, MPI_COMM_WORLD);              
            MPI_Send(available_workers, num_workers_available, MPI_INT, source, FREEWORKER, MPI_COMM_WORLD);
            MPI_Send(&g_lowerBound, 1, MPI_DOUBLE, source, LOWER_BOUND, MPI_COMM_WORLD);     
            break;
        }
    }
    
}

/* WORKER process main routine */
void worker_Bab_Main(MPI_Datatype BabSolutiontype, MPI_Datatype BabNodetype, int rank) {

    Message message;
    MPI_Status status;
    int over = 0;

    // get next subproblem from priority queue
    BabNode *node = Bab_PQPop();

    // save "old" lower bound
    double g_lowerBound = Bab_LBGet();

    /* compute upper bound (SDP bound) and lower bound (via heuristic) for this node */
    node->upper_bound = Evaluate(node, SP, PP, rank);

    // check if better lower bound found --> update info with master
    if (Bab_LBGet() > g_lowerBound){

        message = NEW_VALUE;
        g_lowerBound = Bab_LBGet();

        MPI_Send(&message, 1, MPI_INT, 0, MESSAGE, MPI_COMM_WORLD);
        MPI_Send(&g_lowerBound, 1, MPI_DOUBLE, 0, LOWER_BOUND, MPI_COMM_WORLD);
        MPI_Send(BabSol, 1, BabSolutiontype, 0, SOLUTION, MPI_COMM_WORLD);
        
        MPI_Recv(&g_lowerBound, 1, MPI_DOUBLE, 0, LOWER_BOUND, MPI_COMM_WORLD, &status);

        // update
        BabSolution solx;
        Bab_LBUpd(g_lowerBound, &solx);
    }


    /* if BabLB + 1.0 < child_node->upper_bound, 
     * then we must branch since there could be a better feasible 
     * solution in this subproblem
     */
    if (Bab_LBGet() + 1.0 < node->upper_bound) {

        /***** branch *****/

        // Determine the variable x[ic] to branch on
        int ic = getBranchingVariable(node);

        BabNode *child_node; 
        
        for (int xic = 0; xic <= 1; ++xic) { 

            // Create a new child node from the parent node
            child_node = newNode(node);

            // split on node ic
            child_node->xfixed[ic] = 1;
            child_node->sol.X[ic] = xic;

            /* insert node into the priority queue */
            Bab_PQInsert(child_node);
        }

        // free parent node
        free(node); 

        /************ distribute subproblems ************/

        // leave 1 problem for this worker and the rest is distributed
        int workers_request = heap->used - 1;
        int num_free_workers;
        double g_lowerBound;
        BabSolution solx;

        // check if other subproblems can be send to free workers --> ask master
        message = SEND_FREEWORKERS;
        
        MPI_Send(&message, 1, MPI_INT, 0, MESSAGE, MPI_COMM_WORLD);
        MPI_Send(&workers_request, 1, MPI_INT, 0, FREEWORKER, MPI_COMM_WORLD);
        
        MPI_Recv(&num_free_workers, 1, MPI_INT, 0, NUM_FREE_WORKERS, MPI_COMM_WORLD, &status);
        
        int free_workers[num_free_workers];
        
        MPI_Recv(free_workers, num_free_workers, MPI_INT, 0, FREEWORKER, MPI_COMM_WORLD, &status);
        MPI_Recv(&g_lowerBound, 1, MPI_DOUBLE, 0, LOWER_BOUND, MPI_COMM_WORLD, &status);

        Bab_LBUpd(g_lowerBound, &solx);

        // send subproblems to free workers
        if ( num_free_workers != 0 ) {// free workers found
      
            for (int i = 0; i < num_free_workers; ++i){

                // get next subproblem from queue and send it
                node = Bab_PQPop();

                // send subproblem to free worker
                MPI_Send(&over, 1, MPI_INT, free_workers[i], OVER, MPI_COMM_WORLD);
                MPI_Send(&g_lowerBound, 1, MPI_DOUBLE, free_workers[i], LOWER_BOUND, MPI_COMM_WORLD);
                MPI_Send(node, 1, BabNodetype, free_workers[i], PROBLEM, MPI_COMM_WORLD);

                free(node);
            }    
        }                            

    }
    else {
        // otherwise, intbound <= BabLB, so we can prune
        free(node);
    }

}  




/* print solution 0-1 vector */
void printSolution(FILE *file) {

    fprintf(file, "Solution = ( ");
    for (int i = 0; i < BabPbSize; ++i) {
        if (BabSol->X[i] == 1) {
            fprintf(file, "%d ", i + 1);
        }
    }
    fprintf(file, ")\n");
}


/* print final output */
void printFinalOutput(FILE *file, int num_nodes) {

    // Best solution found
    double best_sol = Bab_LBGet();

    fprintf(file, "\nNodes = %d\n", num_nodes);
    
    // normal termination
    if (!stopped) {
        fprintf(file, "Root node bound = %.2lf\n", root_bound);
        fprintf(file, "Maximum value = %.0lf\n", best_sol);
        
    } else { // B&B stopped early
        fprintf(file, "TIME LIMIT REACHED.\n");
        fprintf(file, "Root node bound = %.2lf\n", root_bound); 
        fprintf(file, "Best value = %.0lf\n", best_sol);
    }

    printSolution(file);
    fprintf(file, "Time = %.2f s\n\n", MPI_Wtime() - TIME);
}


/* Bab function called at the end of the execution.
 * This function frees the memory allocated by the program. */
void Bab_End(void) {
    freeMemory();   
}


/*
 * getBranchingVariable function used in the Bab_GenChild routine to determine
 * which variable x[ic] to branch on.
 *
 * node: the current node of the branch-and-bound search tree
 */
int getBranchingVariable(BabNode *node) {

    int ic = -1;  // x[ic] is the variable to branch on
    double maxValue, minValue;

    /* 
     * Choose the branching variable x[ic] based on params.branchingStrategy
     */
    if (params.branchingStrategy == LEAST_FRACTIONAL) {
        // Branch on the variable x[ic] that has the least fractional value
        maxValue = -BIG_NUMBER;
        for (int i = 0; i < BabPbSize; ++i) {
            if (!(node->xfixed[i]) && fabs(0.5 - node->fracsol[i]) > maxValue) {
                ic = i;
                maxValue = fabs(0.5 - node->fracsol[ic]);
            }
        }
    }
    else if (params.branchingStrategy == MOST_FRACTIONAL) {
        // Branch on the variable x[ic] that has the most fractional value
        minValue = BIG_NUMBER;
        for (int i = 0; i < BabPbSize; ++i) {
            if (!(node->xfixed[i]) && fabs(0.5 - node->fracsol[i]) < minValue) {
                ic = i;
                minValue = fabs(0.5 - node->fracsol[ic]);
            }
        }
    }
    else {
        fprintf(stderr, "Error: Wrong value for params.branchingStrategy\n");
        MPI_Abort(MPI_COMM_WORLD,10);
    }

    return ic;
}


/* Count the number of fixed variables */
int countFixedVariables(BabNode *node) {
    
    int numFixedVariables = 0;

    for (int i = 0; i < BabPbSize; ++i) {
        if (node->xfixed[i]) {
            ++numFixedVariables;
        }
    }

    return numFixedVariables;
}
