#ifndef site_h
#define site_h

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "linked_list.h"

typedef struct site{
    int sigma;
    int Nnode;
    site_node* first;
} site;

site* site_alloc();

void site_free(site* s);

void site_set_sigma(
                site* s, 
                int sigma);

int site_get_sigma(site* s);

site_node* site_get_site_node_before_tau(
                site* s, 
                double tau);

#endif
