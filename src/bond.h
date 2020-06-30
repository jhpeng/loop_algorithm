#ifndef bond_h
#define bond_h

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "linked_list.h"
#include "graph.h"

typedef struct bond{
    int nspin;
    int* site_id;
    int ngraph;
    double* weight;
    graph** graphs;
    int Nnode;
    bond_node* first;
} bond;

#endif
