#define LAB4_EXTEND

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "Lab4_IO.h"
#include "timer.h"

#define EPSILON 0.00001
#define DAMPING_FACTOR 0.85

int main(int argc, char* argv[]) {
    int rank, size;
    struct node *nodehead;
    int nodecount;
    double *r_global, *r_local, *r_prev_global;
    double damping_const;
    int i, j;
    int iterationcount;
    int local_start, local_end, local_count;
    double norm, local_norm;
    double start, end;
    FILE *fp;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        fp = fopen("data_input_meta", "r");
        if (!fp) {
            printf("Error opening the data_input_meta file.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        fscanf(fp, "%d", &nodecount);
        fclose(fp);
    }

    // Broadcast nodecount to all processes
    MPI_Bcast(&nodecount, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (node_init(&nodehead, 0, nodecount)) {
        MPI_Abort(MPI_COMM_WORLD, 2);
    }

    damping_const = (1.0 - DAMPING_FACTOR) / nodecount;

    // Allocate memory
    r_global = malloc(nodecount * sizeof(double));
    r_prev_global = malloc(nodecount * sizeof(double));

    // Each process works on [local_start, local_end)
    local_count = nodecount / size + (rank < nodecount % size ? 1 : 0);
    local_start = (nodecount / size) * rank + (rank < nodecount % size ? rank : nodecount % size);
    local_end = local_start + local_count;
    r_local = malloc(local_count * sizeof(double));

    for (i = 0; i < nodecount; i++)
        r_global[i] = 1.0 / nodecount;

    iterationcount = 0;
    MPI_Barrier(MPI_COMM_WORLD);
    GET_TIME(start);

    do {
        iterationcount++;
        for (i = 0; i < nodecount; i++)
            r_prev_global[i] = r_global[i];

        for (i = local_start; i < local_end; i++) {
            double sum = 0.0;
            for (j = 0; j < nodehead[i].num_in_links; j++) {
                int inbound = nodehead[i].inlinks[j];
                sum += r_prev_global[inbound] / nodehead[inbound].num_out_links;
            }
            r_local[i - local_start] = damping_const + DAMPING_FACTOR * sum;
        }

        // Gather all updated ranks into r_global
        // How many to gather, wher to gather 
        int *recvcounts = malloc(size * sizeof(int));
        int *location = malloc(size * sizeof(int));
        for (i = 0; i < size; i++) {
            recvcounts[i] = nodecount / size + (i < nodecount % size ? 1 : 0);
            location[i] = (nodecount / size) * i + (i < nodecount % size ? i : nodecount % size);
        }

        MPI_Allgatherv(r_local, local_count, MPI_DOUBLE,
                       r_global, recvcounts, location, MPI_DOUBLE,
                       MPI_COMM_WORLD);

        //local norm is same as norm, just easy to synchronize this way
        local_norm = rel_error(r_global, r_prev_global, nodecount);
        MPI_Allreduce(&local_norm, &norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        free(recvcounts);
        free(location);

    } while (norm >= EPSILON);

    MPI_Barrier(MPI_COMM_WORLD);
    GET_TIME(end);

    if (rank == 0) {
        Lab4_saveoutput(r_global, nodecount, end - start);
    }

    node_destroy(nodehead, nodecount);
    free(r_global); free(r_prev_global); free(r_local);
    MPI_Finalize();
    return 0;
}
