#ifndef graph_h
#define graph_h

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <assert.h>

typedef struct graph{
    int nspin;
    int* link;
    int* rule;
} graph;

/*============= define 4-spin graph here ==============*/
//
//             ^
//             |       2|        |3
//             |        |********|
//             |        |--------|
//             |        |--------|
//    temporal |        |********|
//             |       0|        |1
//
//                    -------------->
//                     spatial
/*=====================================================*/

int rule_4spin(int* rule, int* s);

extern graph GRAPH_DIAG;
extern graph GRAPH_HORI;


/*=====================================================*/

typedef struct bond{
    int nspin;
    int* site_id;
    int ngraph;
    graph** graphs;
    double* weight;
    size_t size;
    size_t ntype;
    int* types;
    double* tau;
    int* kink_id;
} bond;

bond* bond_alloc(
                size_t size, 
                int nspin, 
                int ngraph);

void bond_free(
                bond* bd);

void bond_memcpy(
                bond* dest, 
                const bond* src);

void bond_set_site_id(
                bond* bd, 
                int nspin, 
                const int* site_id);

void bond_set_graph(
                bond* bd, 
                int ngraph, 
                graph** graphs, 
                const double* weight);

int bond_check_available_id(
                bond* bd);

int bond_get_nspin(
                bond* bd);

int bond_get_site_id(
                bond* bd, 
                int site);

int bond_get_ngraph(
                bond* bd);

graph* bond_get_graph(
                bond* bd, 
                int type);

double bond_get_weight(
                bond* bd, 
                int type);

size_t bond_get_size(
                bond* bd);

size_t bond_get_ntype(
                bond* bd);

int bond_get_type(
                bond* bd, 
                int type_id);

double bond_get_tau(
                bond* bd, 
                int type_id);

int bond_get_kink_id(
                bond* bd, 
                int type_id, 
                int site);

int bond_insert_graph(
                bond* bd, 
                int type, 
                double tau, 
                int* kink_id);

void bond_remove_graph(
                bond* bd, 
                int type_id);

void bond_print_state(bond* bd);

#endif
