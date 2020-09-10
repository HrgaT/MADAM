/* Max-Heap data structure implementation in C */
/* Used for priority queue for B&B algorithm */

#include "admm.h"

/* definitions of global variables for priority queue */
int BabPbSize;                      // number of nodes in original graph - 1 (last variable fixed to 0)
BabSolution *BabSol;                // global solution of B&B algorithm
static double BabLB;                // global lower bound (use double since int may overflow!)
static int Bab_numNodes = 0;        // number of B&B nodes
BabNode *BabRoot;                   // root node
Heap *heap = NULL;                  // heap is allocated as array of BabNode*


double Bab_LBGet(void) { return BabLB; }
int Bab_numEvalNodes(void) { return Bab_numNodes; }
void Bab_incEvalNodes(void) { ++Bab_numNodes; }


/* Function that determines priority of the BaBNode-s.
 * Priority is based on upper bound: 
 * takes node with higher upper bound (worst bound) first
 *
 * Returns 1 if node1 has bigger priority than node2 and
 *        -1 if other way around
 */
inline int compare_Nodes(const BabNode* node1, const BabNode* node2) {

    return ( (node1->upper_bound > node2->upper_bound) ? 1 : -1 );           
}


inline void swap_entries(int i, int j) {

    BabNode** data = heap->data;
    BabNode* t;

    t       = data[i];
    data[i] = data[j];
    data[j] = t;
}

/* heapify down from root */
static void heapify_down(int current) {

    BabNode** data = heap->data;  
    int child = 2 * current + 1;  // left child

    /* 
     * place element in root in correct position to maintain heap
     */
    if (child + 1 < heap->used) {// right child check
      if (compare_Nodes(data[child + 1], data[child]) > 0)
         child++;
    }
     
    while(child < heap->used && compare_Nodes(data[current], data[child]) < 0)
    {
        swap_entries(current, child);

        current = child;
        child   = 2 * current + 1;
      
        if (child + 1 < heap->used)
            if (compare_Nodes(data[child + 1], data[child]) > 0)
                child++;
    }
}

/* heapify up from last node*/
static void heapify_up(int current) {

    BabNode** data = heap->data;
    int parent = (current-1) / 2;

    while(current > 0 && compare_Nodes(data[parent], data[current]) < 0)
    {
        swap_entries(current, parent);
        current = parent;
        parent  = (current-1) / 2;
    }
}


// initializes heap for storing B&B subproblems via BabNode*
Heap* Init_Heap(int size) {

    Heap *heap;
    alloc(heap, Heap);

    heap->size = size;
    heap->used = 0;
    alloc_vector(heap->data, heap->size, BabNode*);

    return heap;
}


int isPQEmpty(void) {
   return heap->used == 0;
}


BabNode* Bab_PQPop(void) {

   /* safe root, swap it with last node and heapify */     
   BabNode *node = heap->data[0];

   heap->data[0] = NULL;
   heap->used--;

   swap_entries(0, heap->used);
   
   heapify_down(0);
      
   return node;
}


void Bab_PQInsert(BabNode *node) {
   
    // check heap size
    if (heap->size == heap->used) {
        puts("\nERROR: Maximum size of heap reached.\n");
        exit(1);
    }
   
   /* place new node at the end of heap and heapify */
   heap->data[heap->used] = node;
   heap->used++;

   heapify_up(heap->used - 1);
}


/*
 * Create a new B&B node.
 *
 * If parentNode == NULL, it will create the root node.
 * Otherwise, the new node will be a child of parentNode.
 */
// NOTE: Bab_GenChild will place created child node in priority queue
BabNode* newNode(BabNode *parentNode) {

    // allocate memory for the new child node
    BabNode *node = (BabNode *) malloc(sizeof(BabNode));
    if (node == NULL) {
        fprintf(stderr, "Error: Not enough memory for creating new node.\n");
        exit(1);
    }

    // copy the solution information from the parent node
    for (int i = 0; i < BabPbSize; ++i) {
        if (parentNode == NULL) {
            node->xfixed[i] = 0;
            node->sol.X[i] = 0;
        }
        else {
            node->xfixed[i] = parentNode->xfixed[i];
            node->sol.X[i] = (node->xfixed[i]) ? parentNode->sol.X[i] : 0;
        }
    }

    // child is one level deeper than parent
    node->level = (parentNode == NULL) ? 0 : parentNode->level + 1;

    return node;
}


/* Allocate and initialize global lower bound and solution vector */
void Bab_LBInit(double lowerBound, BabSolution *bs) {

    BabSol = (BabSolution *) malloc(sizeof(BabSolution));
    if (BabSol == NULL) {
        fprintf(stderr, "Not enough memory for BabSol.\n");
        exit(1);
    }
    *BabSol = *bs;
    BabLB = lowerBound;
}


/* If new solution is better than the global solution, update the solution */
int Bab_LBUpd(double new_LB, BabSolution *bs) {

    if (new_LB > BabLB) {
        BabLB = new_LB;
        *BabSol = *bs;
        return 1;
    }
    return 0;
}
