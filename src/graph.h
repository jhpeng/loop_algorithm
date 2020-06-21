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
    int ntype;
    graph** type;
    double* weight;
    size_t size;
    size_t ngraph;
    int* graphs;
    double* tau;
    int* kink_id;
} bond;

bond* bond_alloc(
                size_t size, 
                int nspin, 
                int ntype);

void bond_free(
                bond* bd);

void bond_memcpy(
                bond* dest, 
                const bond* src);

void bond_set_site_id(
                bond* bd, 
                int nspin, 
                const int* site_id);

void bond_set_graph_type(
                bond* bd, 
                int ntype, 
                graph** type, 
                const double* weight);

int bond_check_available_id(
                bond* bd);

int bond_insert_graph(
                bond* bd, 
                int type, 
                double tau, 
                int* kink_id);

void bond_remove_graph(
                bond* bd, 
                int graph_id);
#endif
