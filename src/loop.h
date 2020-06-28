#ifndef loop_h
#define loop_h

#include <stdio.h>
#include <stdlib.h>
#define NDEBUG
#include <assert.h>

#include <gsl/gsl_rng.h>

#include "configuration.h"
#include "graph.h"


typedef struct link{
    size_t size;
    size_t nbond;
    size_t max_type_id;
    int max_nspin;
    int* index;             //  max_nspin*(type_id*nbond+bond_id)+spin_id
} link;

link* link_alloc(
                size_t nbond, 
                size_t max_type_id, 
                int max_nspin);

void link_realloc(
                link* lk, 
                size_t max_type_id);

void link_free(link* lk);

int link_check_init(link* lk);

void link_measure_size(
                bond** bd, 
                size_t nbond, 
                size_t* max_type_id, 
                int* mas_nspin);

void loop_construct_outer_link(
                link* lk, 
                kinks** ks, 
                bond** bd, 
                size_t nsite, 
                size_t nbond);

void loop_construct_inner_link(
                link* lk, 
                bond** bd, 
                size_t nbond);

void loop_cluster_identify(
                link* outer_lk, 
                link* inner_lk, 
                gsl_rng* rng);

void loop_cluster_flip(
                kinks** ks, 
                bond** bd, 
                link* outer_lk, 
                int nbond);

void loop_cluster_sigma_i(
                kinks** ks, 
                int nsite, 
                gsl_rng* rng);

int loop_cluster_check_periodic(
                kinks** ks, 
                int nsite);

void loop_cluster_update_user_friendly(
                kinks** ks, 
                bond** bd, 
                int nsite, 
                int nbond, 
                gsl_rng* rng);

#endif
