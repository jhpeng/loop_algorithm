#ifndef update_graphs_h
#define update_graphs_h

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <assert.h>
#include <math.h>

#include <gsl/gsl_rng.h>

#include "configuration.h"
#include "graph.h"

void generate_graphs_with_uniform_dist(
                kinks** ks,                     //bounch of kinks 
                bond** bd,                      //bounch of bond
                int bond_id,                    //choose a bond to generate graphs
                double beta,                    //inverse temperature of the system
                double** taus,                  //working space of uniform dist. 
                size_t* size,                   //size of working space
                gsl_rng* rng,                   //gsl random number generater
                int (*rule_nspin)(int*,int*));  //the rule for nspin graph

void remove_all_graphs_with_no_kink(
                kinks** ks,                     //bounch of kinks
                bond** bd,                      //bounch of bond
                int nbond);                     //total number of bond
#endif
