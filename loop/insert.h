#ifndef insert_h
#define insert_h

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#define NDEBUG
#include <assert.h>

#include "chain.h"
#include "table.h"

#define INSERT_MAX 100000

void insert_horizontal_graph(
                chain* c1, 
                chain* c2, 
                table* t, 
                double w, 
                double beta, 
                gsl_rng* rng);

void insert_triangular_cut_graph(
                chain** c, 
                table* t, 
                double w, 
                double beta, 
                gsl_rng* rng);

void insert_triangular_graph(
                chain** c, 
                table* t, 
                double w, 
                double beta, 
                gsl_rng* rng);

#endif
