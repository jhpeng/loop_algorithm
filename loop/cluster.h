#ifndef cluster_h
#define cluster_h

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "chain.h"
#include "table.h"

void cluster_update_table(table* t);

void claster_update_chain(chain* c, table* t);

void cluster_link_vertex(chain* c, table* t);

#endif
