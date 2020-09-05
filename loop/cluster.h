#ifndef cluster_h
#define cluster_h

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#define NDEBUG
#include <assert.h>

#include "chain.h"
#include "table.h"

void cluster_update_table(table* t);

void cluster_update_chain(chain* c, table* t);

void cluster_link_vertex(chain* c, table* t);

void cluster_traverse(table* t, gsl_rng* rng);

#endif
