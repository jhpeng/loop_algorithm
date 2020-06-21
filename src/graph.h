#ifndef graph_h
#define graph_h

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>

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

int rule_4spin(int* rule, int s0, int s1, int s2, int s3);

//diagonal graph  
int GRAPH_LINK_DIAG[4]  = {3,2,1,0};
int GRAPH_RULE_DIAG[16] = {1,0,0,0,
                           0,0,1,0,
                           0,1,0,0,
                           0,0,0,1};

graph GRAPH_DIAG = {4,GRAPH_LINK_DIAG,GRAPH_RULE_DIAG};

//horizontal graph  
int GRAPH_LINK_HORI[4]  = {1,0,3,2};
int GRAPH_RULE_HORI[16] = {0,0,0,0,
                           0,1,1,0,
                           0,1,1,0,
                           0,0,0,0};

graph GRAPH_HORI = {4,GRAPH_LINK_HORI,GRAPH_RULE_HORI};


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

int bond_insert_graph(
                bond* bd, 
                int type, 
                double tau, 
                int* kink_id);

void bond_remove_graph(
                bond* bd, 
                int graph_id);
#endif
