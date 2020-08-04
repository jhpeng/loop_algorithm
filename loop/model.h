#ifndef model_h
#define model_h

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "chain.h"
#include "table.h"
#include "insert.h"
#include "loops.h"
#include "cluster.h"


typedef struct model{
    int nsite;
    int nbond;
    double beta;
    int* type;
    int* bond2site;     // site_id = bond2site[bond_id*NSPIN_MAX+spin_id]
    double* weight;
} model;

#endif 
