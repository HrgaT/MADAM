#include <stdio.h>
#include <mpi.h>

#include "admm.h"  

#define HEAP_SIZE 1000000
extern Heap *heap;
extern double diff;

extern Parameters params;
extern double TIME;
extern FILE *output;

int num_workers_used = 0;

int main(int argc, char **argv) {

    /*******************************************************
    *********** BRANCH & BOUND: PARALLEL ALGORITHM ********
    ******************************************************/

    // number of processes = master + workers
    int numbWorkers;

    // rank of each process: from 0 to numWorkers-1
    int rank;

    // MPI Start: start parallel environment
    MPI_Init(&argc, &argv);

    // get number of proccesses and corresponding ranks
    MPI_Comm_size(MPI_COMM_WORLD, &numbWorkers);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;

    if (rank == 0)
	   printf("Number of cores: %d\n", numbWorkers);
	
    /***** user defined MPI struct: for sending and receiving *****/
    // (1) for BabSolution
    MPI_Datatype BabSolutiontype;
    MPI_Datatype type1[1] = { MPI_INT };
    int blocklen1[1] = { NMAX };
    MPI_Aint disp1[1];
    disp1[0] = offsetof(BabSolution, X);
    MPI_Type_create_struct(1, blocklen1, disp1, type1, &BabSolutiontype);
    MPI_Type_commit(&BabSolutiontype);

    // (2) for BabNode
    MPI_Datatype BabNodetype;
    MPI_Datatype type2[5] = { MPI_INT, BabSolutiontype, MPI_DOUBLE, MPI_INT, MPI_INT };
    int blocklen2[5] = { NMAX, 1, NMAX, 1, 1 };
    MPI_Aint disp2[5];
    disp2[0] = offsetof(BabNode, xfixed);
    disp2[1] = offsetof(BabNode, sol);
    disp2[2] = offsetof(BabNode, fracsol);
    disp2[3] = offsetof(BabNode, level);
    disp2[4] = offsetof(BabNode, upper_bound);
    MPI_Type_create_struct(5, blocklen2, disp2, type2, &BabNodetype);
    MPI_Type_commit(&BabNodetype);
    /***********************************/

    // Start the timer
    TIME = MPI_Wtime();

    // type of message
    Message message;

    // tag to FINISH set to false
    int over = 0;

    // helper variables
    BabNode *node;
    double g_lowerBound;

    /* each process allocates its local priority queue */
    heap = Init_Heap(HEAP_SIZE);

    /* every process reads params and initializes B&B solution,
     * only master process creates output file, reads input graph
     * and broadcast it */
    int read_error = Bab_Init(argc, argv, rank);
	
    if (read_error)
        goto FINISH;
	
	

    /******************** MASTER PROCESS ********************/
    if (rank == 0)
    {

        // only master evaluates the root node
        // and places it in priority queue if not able to prune
        over = Init_PQ();

	printf("Initial lower bound: %.0lf\n", Bab_LBGet());    
	
		// broadcast diff
	if (params.use_diff)
	    MPI_Bcast(&diff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);	

        // broadcast lower bound to others or -1 to exit
	MPI_Bcast(&over, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if ( (over == -1) || params.root) {          
            goto FINISH;
        }
        else {
            g_lowerBound = Bab_LBGet();
            MPI_Bcast(&g_lowerBound, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }


        // array of busy workers: 0 = free, 1 = busy
        // only master is busy
        int busyWorkers[numbWorkers];
        busyWorkers[0] = 1;
        for (int i = 1; i < numbWorkers; ++i)
            busyWorkers[i] = 0;

        int numbFreeWorkers = numbWorkers - 1;
        int source;


        /***** branch root node and send to workers *****/
        node = Bab_PQPop();

        // Determine the variable x[ic] to branch on
        int ic = getBranchingVariable(node);

        // send two nodes to workers 1 and 2
        BabNode *child_node;
        int worker;
        
        for (int xic = 0; xic <= 1; ++xic) { 

            // Create a new child node from the parent node
            child_node = newNode(node);

            // split on node ic
            child_node->xfixed[ic] = 1;
            child_node->sol.X[ic] = xic;

            // increment the number of explored nodes
            Bab_incEvalNodes();

            worker = xic + 1;
            busyWorkers[worker] = 1;
            --numbFreeWorkers;

            MPI_Send(&over, 1, MPI_INT, worker, OVER, MPI_COMM_WORLD);
	    MPI_Send(&g_lowerBound, 1, MPI_DOUBLE, worker, LOWER_BOUND, MPI_COMM_WORLD);
            MPI_Send(child_node, 1, BabNodetype, worker, PROBLEM, MPI_COMM_WORLD);

            free(child_node);
        }

        // free parent nodes
        free(node);    

	    num_workers_used = 2;

	
        /************* MAIN LOOP for master **************/
        do {

            /*** wait for messages: extract source from status ***/
            MPI_Recv(&message, 1, MPI_INT, MPI_ANY_SOURCE, MESSAGE, MPI_COMM_WORLD, &status);
            source = status.MPI_SOURCE;

            master_Bab_Main(message, source, busyWorkers, numbWorkers, &numbFreeWorkers, BabSolutiontype);

        } while ( numbFreeWorkers != numbWorkers - 1 );
        /*************************************************/

        // send over messages to the workers
        over = 1;
        for(int i = 1; i < numbWorkers; ++i) {
            MPI_Send(&over, 1, MPI_INT, i, OVER, MPI_COMM_WORLD);
        }

    }
     /******************** WORKER PROCESS ********************/
    else
    {
		// receive diff
	if (params.use_diff)
	    MPI_Bcast(&diff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// receive over (stop or continue)
	MPI_Bcast(&over, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // receive lower bound
	if (over == -1 || params.root )   // root node is pruned
	    goto FINISH;
	else
            MPI_Bcast(&g_lowerBound, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// update lower bound
	BabSolution solx;
	Bab_LBUpd(g_lowerBound, &solx);
        
        /************* MAIN LOOP for worker **************/
        do {

            // wait for info: stop (from master) or receive new subproblem from other worker
            MPI_Recv(&over, 1, MPI_INT, MPI_ANY_SOURCE, OVER, MPI_COMM_WORLD, &status);

            if (!over) {

		alloc(node, BabNode);

                // receive subproblem from master or other worker
		        MPI_Recv(&g_lowerBound, 1, MPI_DOUBLE, MPI_ANY_SOURCE, LOWER_BOUND, MPI_COMM_WORLD, &status);
                MPI_Recv(node, 1, BabNodetype, MPI_ANY_SOURCE, PROBLEM, MPI_COMM_WORLD, &status);
		
		        // update
		        Bab_LBUpd(g_lowerBound, &solx);

                // start local queue
                Bab_PQInsert(node);

                while(!isPQEmpty()){

                    worker_Bab_Main(BabSolutiontype, BabNodetype, rank);
                }    

                message = IDLE;
                MPI_Send(&message, 1, MPI_INT, 0, MESSAGE, MPI_COMM_WORLD);
            }
            
        } while (over != 1);

        //free(node);
    }

    FINISH:

    /* Print results to the standard output and to the output file */
    if (rank == 0) {
        printFinalOutput(stdout,Bab_numEvalNodes());
        printFinalOutput(output,Bab_numEvalNodes());
	fprintf(output, "Number of cores: %d\n", numbWorkers);
	fprintf(output, "Maximum number of workers used: %d\n", num_workers_used);
	printf("Maximum number of workers used: %d\n", num_workers_used);
        fclose(output);
    }

    /* free memory */
    Bab_End();

    free(heap->data);
    free(heap);

    // MPI finish
    MPI_Finalize();

    return 0;
}
