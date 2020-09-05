#ifndef loops_h
#define loops_h

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#define NDEBUG
#include <assert.h>

#include "chain.h"
#include "table.h"

void loops_update_table(table* t);

void loops_update_chain(
                chain* c, 
                table* t, 
                gsl_rng* rng);

void loops_link_vertex(
                chain* c, 
                table* t);

void loops_traverse(
                table* t, 
                gsl_rng* rng);

#endif
