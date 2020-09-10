#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>   // for numbering of output files
#include <mpi.h>

#include "admm.h"

extern FILE *output;
extern Parameters params;
extern Problem *SP;             
extern Problem *PP;            
extern int BabPbSize;

// macro to handle the errors in the input reading
#define READING_ERROR(file,cond,message)\
        if ((cond)) {\
            fprintf(stderr, "\nError: "#message"\n");\
            fclose(file);\
            return 1;\
        }


void print_symmetric_matrix(double *Mat, int N) {

    double val;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            val = (i >= j) ? Mat[i + j*N] : Mat[j + i*N];
            printf("%24.16e", val);
        }
        printf("\n");
    }
}        

int processCommandLineArguments(int argc, char **argv, int rank) {

    int read_error = 0;

    if (argc != 3) {
        if (rank == 0)
            fprintf(stderr, "Usage: ./biqbin file.rudy file.params\n");
        read_error = 1;
        return read_error;
    }

    /***** only master process creates output file and reads the whole graph *****/

    // Control the command line arguments
    if (rank == 0) {

        // Create the output file
        char output_path[200];
        sprintf(output_path, "%s.output", argv[1]);

        // Check if the file already exists, if so aappend _<NUMBER> to the end of the output file name
        struct stat buffer;
        int counter = 1;
        
        while (stat(output_path, &buffer) == 0)
            sprintf(output_path, "%s.output_%d", argv[1], counter++);

        output = fopen(output_path, "w");
        if (!output) {
            fprintf(stderr, "Error: Cannot create output file.\n");
            read_error = 1;
            MPI_Bcast(&read_error, 1, MPI_INT, 0, MPI_COMM_WORLD);
            return read_error;
        }

        // Read the input file instance
        read_error = readData(argv[1]);

        // bcast first read_error then whole graph
        MPI_Bcast(&read_error, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (read_error)
            return read_error;
        else {
            MPI_Bcast(&(SP->n), 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(SP->L, SP->n * SP->n, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
        }
            
    }
    else {

        MPI_Bcast(&read_error, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (read_error) 
            return read_error;    

        // allocate memory for original problem SP and subproblem PP
        alloc(SP, Problem);
        alloc(PP, Problem);

        MPI_Bcast(&(SP->n), 1, MPI_INT, 0, MPI_COMM_WORLD);

        // allocate memory for objective matrices for SP and PP
        alloc_matrix(SP->L, SP->n, double);
        alloc_matrix(PP->L, SP->n, double);

        MPI_Bcast(SP->L, SP->n * SP->n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // IMPORTANT: last node is fixed to 0
        // --> BabPbSize is one less than the size of problem SP
        BabPbSize = SP->n - 1; // num_vertices - 1;
        PP->n = SP->n;

        int N2 = SP->n * SP->n;
        int incx = 1;
        int incy = 1;
        dcopy_(&N2, SP->L, &incx, PP->L, &incy);
    }    
    
    // Read the parameters from a user file
    read_error = readParameters(argv[2], rank);
    if (read_error)
        return read_error;

    /* adjust parameters */
    // change number of added cutting planes per iteration to n*10
    if (params.adjust_TriIneq)
        params.TriIneq = SP->n * 10;


    if (rank == 0) {
        // print parameters to output file
        fprintf(output, "Parameters:\n");
#define P(type, name, format, def_value)\
            fprintf(output, "%20s = "format"\n", #name, params.name);
        PARAM_FIELDS
#undef P
    }

    return read_error;
}



/* Read parameters contained in the file given by the argument */
int readParameters(const char *path, int rank) {

    FILE* paramfile;
    char s[128];            // read line
    char param_name[50];

    // Initialize every parameter with its default value
#define P(type, name, format, def_value)\
    params.name = def_value;
    PARAM_FIELDS
#undef P

    // open parameter file
    if ( (paramfile = fopen(path, "r")) == NULL ) {
        if (rank == 0)
            fprintf(stderr, "Error: parameter file %s not found.\n", path);
        return 1;
    }
    
    while (!feof(paramfile)) {
        if ( fgets(s, 120, paramfile) != NULL ) {
        
            // read parameter name
            sscanf(s, "%[^=^ ]", param_name);

            // read parameter value
#define P(type, name, format, def_value)\
            if(strcmp(#name, param_name) == 0)\
            sscanf(s, "%*[^=]="format"\n", &(params.name));
            PARAM_FIELDS
#undef P

        }
    }
    fclose(paramfile);   
   
    return 0;
}


/*** read graph file ***/
int readData(const char *instance) {

    // open input file
    FILE *f = fopen(instance, "r");
    if (f == NULL) {
        fflush(stdout);
        fprintf(stderr, "Error: problem opening input file %s\n", instance);
        return 1;
    }
    printf("Input file: %s\n", instance);	
    fprintf(output,"Input file: %s\n", instance);

    int num_vertices;
    int num_edges;

    READING_ERROR(f, fscanf(f, "%d %d \n", &num_vertices, &num_edges) != 2,
                  "Problem reading number of vertices and edges");
    READING_ERROR(f, num_vertices <= 0, "Number of vertices has to be positive");

    // OUTPUT information on instance
    fprintf(stdout, "\nGraph has %d vertices and %d edges.\n", num_vertices, num_edges);
    fprintf(output, "\nGraph has %d vertices and %d edges.\n", num_vertices, num_edges);

    // read edges and store them in matrix Adj
    // NOTE: last node is fixed to 0
    int i, j;
    double weight;


    // Adjacency matrix Adj: allocate and set to 0 
    double *Adj;
    alloc_matrix(Adj, num_vertices, double);

    for (int edge = 0; edge < num_edges; ++edge) {
        
        READING_ERROR(f, fscanf(f, "%d %d %lf \n", &i, &j, &weight) != 3,
                      "Problem reading edges of the graph"); 

        READING_ERROR(f, ((i < 1 || i > num_vertices) || (j < 1 || j > num_vertices)),
                      "Problem with edge. Vertex not in range");  
        
        Adj[ num_vertices * (j - 1) + (i - 1) ] = weight;
        Adj[ num_vertices * (i - 1) + (j - 1) ] = weight;      
    }   

    fclose(f);

    // allocate memory for original problem SP and subproblem PP
    alloc(SP, Problem);
    alloc(PP, Problem);

    // size of matrix L
    SP->n = num_vertices;                 

    // allocate memory for objective matrices for SP and PP
    alloc_matrix(SP->L, SP->n, double);
    alloc_matrix(PP->L, SP->n, double);


    // IMPORTANT: last node is fixed to 0
    // --> BabPbSize is one less than the size of problem SP
    BabPbSize = SP->n - 1; // num_vertices - 1;
    PP->n = SP->n;
    

    /********** construct SP->L from Adj **********/
    /*
     * SP->L = [ Laplacian,  Laplacian*e; (Laplacian*e)',  e'*Laplacian*e]
     */
    // NOTE: we multiply with 1/4 afterwards when subproblems PP are created!
    //       (in function createSubproblem)
    // NOTE: Laplacian is stored in upper left corner of L

    // (1) construct vec Adje = Adj*e 
    double *Adje;
    alloc_vector(Adje, num_vertices, double);

    for (int ii = 0; ii < num_vertices; ++ii) {
        for (int jj = 0; jj < num_vertices; ++jj) {
            Adje[ii] += Adj[jj + ii * num_vertices];
        }
    }

    // (2) construct Diag(Adje)
    double *tmp;
    alloc_matrix(tmp, num_vertices, double);
    Diag(tmp, Adje, num_vertices);

    // (3) fill upper left corner of L with Laplacian = tmp - Adj,
    //     vector parts and constant part      
    double sum_row = 0.0;
    double sum = 0.0;

    // NOTE: skip last vertex!!
    for (int ii = 0; ii < num_vertices; ++ii) {            
        for (int jj = 0; jj < num_vertices; ++jj) {

            // matrix part of L
            if ( (ii < num_vertices - 1) && (jj < num_vertices - 1) ) {
                SP->L[jj + ii * num_vertices] = tmp[jj + ii * num_vertices] - Adj[jj + ii * num_vertices]; 
                sum_row += SP->L[jj + ii * num_vertices];       
            }
            // vector part of L
            else if ( (jj == num_vertices - 1) && (ii != num_vertices - 1)  ) {
                SP->L[jj + ii * num_vertices] = sum_row;
                sum += sum_row;
            }
            // vector part of L
            else if ( (ii == num_vertices - 1) && (jj != num_vertices - 1)  ) {
                SP->L[jj + ii * num_vertices] = SP->L[ii + jj * num_vertices];
            }
            // constant term in L
            else { 
                SP->L[jj + ii * num_vertices] = sum;
            }
        }
        sum_row = 0.0;
    } 

    int N2 = SP->n * SP->n;
    int incx = 1;
    int incy = 1;
    dcopy_(&N2, SP->L, &incx, PP->L, &incy);


    free(Adj);
    free(Adje);
    free(tmp);  

    return 0;
}
