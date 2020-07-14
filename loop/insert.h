#ifndef insert_h
#define insert_h

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "chain.h"
#include "table.h"

void insert_horizontal_graph(
                chain* c1, 
                chain* c2, 
                table* t, 
                double w, 
                double beta, 
                gsl_rng* rng);

#endif
