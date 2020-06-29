#ifndef lattice_h
#define lattice_h

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <gsl/gsl_rng.h>

#include "configuration.h"
#include "graph.h"
#include "update_graphs.h"
#include "loop.h"


typedef struct lattice{
    int dim;
    int* shape;
    int max_nspin;
    int nsite;
    int nbond;
    kinks** ks;
    bond** bd;
    double beta;
    gsl_rng* rng;
} lattice;

lattice* lattice_afm_heisenberg_2d_uniform(
            int lx, 
            int ly, 
            int seed);

void lattice_set_beta(
            lattice* lat, 
            double beta);

void lattice_monte_carlo_sweep(
            lattice* lat, 
            int nsweep);

void lattice_free(lattice* lat);

#endif
