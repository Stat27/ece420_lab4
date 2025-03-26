//
// Created by kid1412 on 3/23/25.
//
/*
Serial Implementation of Lab 4
*/

#define LAB4_EXTEND

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Lab4_IO.h"
#include "timer.h"

#define EPSILON 0.00001
#define DAMPING_FACTOR 0.85

int main (int argc, char* argv[]){
    // instantiate variables
    struct node *nodehead;
    int nodecount;
    double *r, *r_pre;
    int i, j;
    int iterationcount;
    double start, end;
    FILE *ip;
    /* INSTANTIATE MORE VARIABLES IF NECESSARY */
    double damping_const;
    double norm = 0.0;
    GET_TIME(start);

    // load data
    if ((ip = fopen("data_input_meta","r")) == NULL) {
        printf("Error opening the data_input_meta file.\n");
        return 253;
    }
    fscanf(ip, "%d\n", &nodecount);
    fclose(ip);
   GET_TIME(start);

    if (node_init(&nodehead, 0, nodecount)) return 254;

    // initialize variables
    r = malloc(nodecount * sizeof(double));
    r_pre = malloc(nodecount * sizeof(double));

    // (1 - d) * 1/N
    damping_const = (1.0 - DAMPING_FACTOR) / nodecount;

    // initial pagerank
    iterationcount = 0;
    for ( i = 0; i < nodecount; ++i)
        r[i] = 1.0 / nodecount;
    /* INITIALIZE MORE VARIABLES IF NECESSARY */

    // core calculation
    do{
        ++iterationcount;

        for (i = 0; i < nodecount; ++i) {
            r_pre[i] = r[i];
        }

        for (i = 0; i < nodecount; i++) {
            double sum = 0.0;
            for (j = 0; j < nodehead[i].num_in_links; j++) {
                int inbound = nodehead[i].inlinks[j];
                sum += r_pre[inbound] / nodehead[inbound].num_out_links;
            }
            r[i] = damping_const + DAMPING_FACTOR * sum;
        }

        norm = rel_error(r, r_pre, nodecount);
        /* IMPLEMENT ITERATIVE UPDATE */

    }while(norm >= EPSILON);

    GET_TIME(end);

    Lab4_saveoutput(r, nodecount, end - start);

    // post processing
    node_destroy(nodehead, nodecount);
    free(r); free(r_pre);
    return 0;
}