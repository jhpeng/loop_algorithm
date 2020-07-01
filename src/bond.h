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

bond* bond_alloc(
                int ngraph, 
                int nspin);

void bond_free(
                bond* b);

void bond_set_site_id(
                bond* b, 
                int* site_id, 
                int nspin);

int bond_get_site_id(
                bond* b, 
                int spin_id);

void bond_set_graphs(
                bond* b, 
                graph** graphs, 
                int ngraph);

graph* bond_get_graphs(
                bond* b, 
                int graph_id);

void bond_set_weight(
                bond* b, 
                double* weight, 
                int ngraph);

double bond_get_weight(
                bond* b, 
                int graph_id);

int bond_get_nspin(bond* b);

int bond_get_ngraph(bond* b);

bond_node* bond_create_new_node(
                bond* b, 
                int graph_id);

void bond_remove_node(
                bond* b, 
                bond_node* bnode);

#endif
