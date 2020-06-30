#include "bond.h"

bond* bond_alloc(int ngraph, int nspin){
    bond* b = (bond*)malloc(sizeof(bond));
    b->nspin   = nspin;
    b->site_id = (int*)malloc(sizeof(int)*nspin);
    b->ngraph  = ngraph;
    b->weight  = (double*)malloc(sizeof(double)*ngraph);
    b->graphs  = (graph**)malloc(sizeof(graph*)*ngraph);
    b->Nnode   = 0;
    b->first   = NULL;

    return b;
}

void bond_free(bond* b){
    free(b->site_id);
    free(b->weight);
    free(b->graphs);

    b->first = bond_node_remove_all(b->first);
    free(b);
}

void bond_set_site_id(bond* b, int* site_id, int nspin){
    assert(b!=NULL);

    for(int i=0;i<nspin;++i)
        b->site_id[i] = site_id[i];
}

int bond_get_site_id(bond* b, int spin_id){
    assert(b!=NULL);
    assert(spin_id<(b->nspin));

    return b->site_id[spin_id];
}

void bond_set_graphs(bond* b, graph** graphs, int ngraph){
    assert(b!=NULL);
    
    for(int i=0;i<ngraph;++i)
        b->graphs[i] = graphs[i];
}

graph* bond_get_graphs(bond* b, int graph_id){
    assert(b!=NULL);
    assert(graph_id<(b->ngraph));

    return b->graphs[graph_id];
}

void bond_set_weight(bond* b, double* weight, int ngraph){
    assert(b!=NULL);

    for(int i=0;i<ngraph;++i)
        b->weight[i] = weight[i];
}

double bond_get_weight(bond* b, int graph_id){
    assert(b!=NULL);
    assert(graph_id<(b->ngraph));

    return b->weight[graph_id];
}

int bond_get_nspin(bond* b){
    assert(b!=NULL);

    return b->nspin;
}

int bond_get_ngraph(bond* b){
    assert(b!=NULL);

    return b->ngraph;
}

bond_node* bond_create_new_node(bond* b, int graph_id){
    assert(b!=NULL);
    assert(graph_id<bond_get_ngraph(b));

    b->first = bond_node_insertFirst(b->first,graph_id,bond_get_nspin(b));
    b->Nnode++;

    return b->first;
}

void bond_remove_node(bond* b, bond_node* bnode){
    assert(b!=NULL);
    
    b->first = bond_node_remove(b->first,bnode);
    b->Nnode--;
}
