#### MADAM makefile ####

# Directories
OBJ = Obj
SUITESPARSE = SuiteSparse

# Compilers
CC = mpicc

INCLUDES = -I./$(SUITESPARSE)/include
LIB 	 = -L./$(SUITESPARSE)/lib -lcholmod -lamd -lcolamd -lsuitesparseconfig -lccolamd -lcamd $(SUITESPARSE)/lib/libmetis.so
LINALG 	 = -lopenblas -lm
OPTI     = -O3 -ffast-math -fexceptions -fPIC -fno-common

# binary
BINS =  admm

# Objects
BBOBJS = $(OBJ)/admm.o $(OBJ)/allocate_free.o $(OBJ)/bab_functions.o \
	 $(OBJ)/bounding.o $(OBJ)/cutting_planes.o \
         $(OBJ)/evaluate.o $(OBJ)/heap.o $(OBJ)/ipm_mc_pk.o \
         $(OBJ)/heuristic.o $(OBJ)/main.o $(OBJ)/operators.o \
         $(OBJ)/process_input.o $(OBJ)/qap_simulated_annealing.o

# All objects
OBJS = $(BBOBJS)

CFLAGS = $(OPTI) -Wall -W -pedantic 


#### Rules ####

.PHONY : all clean run

# Default rule is to create all binaries #
all: $(BINS)

# Rules for binaries #
$(BINS) : $(OBJS)
	$(CC) -o $@ $^ $(INCLUDES) $(LIB) $(OPTI) $(LINALG)  


# code rules 
$(OBJ)/%.o : %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

run : $(BINS)
	./admm rudy/g05_60.0 params
		

# Clean rule #
clean :
	rm $(BINS) $(OBJS)