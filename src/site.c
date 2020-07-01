#include "site.h"

site* site_alloc(){
    site* s = (site*)malloc(sizeof(site));

    if(s==NULL){
        printf("site_alloc : memory allocate fail!\n");
        exit(1);
    }

    s->Nnode = 0;
    s->first = NULL;

    return s;
}

void site_free(site* s){
    s->first = site_node_remove_all(s->first);
    free(s);
}

void site_set_sigma(site* s, int sigma){
    s->sigma = sigma;
}

int site_get_sigma(site* s){
    assert(s!=NULL);

    return s->sigma;
}

site_node* site_get_site_node_before_tau(site* s, double tau){
    if(s->first==NULL) return NULL;

    site_node* node;
    node = s->first;

    while((node->tau)<tau){
        if((node->next)==NULL) return node;
        node = node->next;
    }
    
    return node->prev;
}
