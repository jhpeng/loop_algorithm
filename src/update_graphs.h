#ifndef update_graphs_h
#define update_graphs_h

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>
#define NDEBUG
#include <assert.h>
#include <math.h>

#include <gsl/gsl_rng.h>

#include "configuration.h"
#include "graph.h"

typedef struct insertion_plan{
    int size;
    int max_nspin;
    int ntau;
    double* taus;
    int* accept;
    int* sigma;
} insertion_plan;

void remove_all_graphs_with_no_kink(
                kinks** ks,                     //bounch of kinks
                bond** bd,                      //bounch of bond
                int nsite,                      //total number of site
                int nbond);                     //total number of bond

void update_graph_user_friendly(
                kinks** ks, 
                bond** bd, 
                int nsite, 
                int nbond, 
                int max_nspin, 
                double beta, 
                gsl_rng* rng);
#endif
