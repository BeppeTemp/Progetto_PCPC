#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mpi.h"

#define MASTER 0

void main() {
    int world_size, rank;
    double start, end;

    // Inizializzazione MPI
    MPI_Status status;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Esecuzione
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    int *test = malloc(sizeof(int));

    while (1){
        test =  realloc(test, sizeof(int) * 20);
    }

    free(test);
    
    MPI_Barrier(MPI_COMM_WORLD);

    end = MPI_Wtime();
    MPI_Finalize();

    if (rank == MASTER) printf("\n\n🕒 Time in ms = %f\n", end - start);
}