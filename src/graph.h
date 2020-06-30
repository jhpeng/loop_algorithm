#ifndef graph_h
#define graph_h

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

typedef struct graph{
    int nspin;
    int* link;
    int* rule;
} graph;

int rule_4spin(int* rule, int* s);

extern graph GRAPH_DIAG;
extern graph GRAPH_HORI;

#endif
