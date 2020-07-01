#include "linked_list.h"

bond_node* bond_node_alloc(int nspin){
    bond_node* bn = (bond_node*)malloc(sizeof(bond_node));

    if(bn==NULL){
        printf("bond_node_alloc : alocate memory fail!\n");
        exit(1);
    }

    bn->nspin = nspin;
    bn->leg   = (int*)malloc(sizeof(int)*2*nspin);
    bn->sigma = (int*)malloc(sizeof(int)*2*nspin);
    bn->termination = 0;
    bn->prev = NULL;
    bn->next = NULL;

    if(bn->leg==NULL || bn->sigma==NULL){
        printf("bond_node_alloc : alocate memory fail!\n");
        exit(1);
    }

    return bn;
}

void bond_node_free(bond_node* bn){
    free(bn->leg);
    free(bn->sigma);
    free(bn);
}

site_node* site_node_alloc(){
    site_node* sn = (site_node*)malloc(sizeof(site_node));
    sn->prev = NULL;
    sn->next = NULL;
    sn->cont = NULL;

    if(sn==NULL){
        printf("site_node_alloc : alocate memory fail!\n");
        exit(1);
    }

    return sn;
}

void site_node_free(site_node* sn){
    free(sn);
}

bond_node*  bond_node_insertFirst(bond_node* first, int graph_id, int nspin){
    bond_node* bn = bond_node_alloc(nspin);
    bn->graph_id = graph_id;
    if(first!=NULL){
        assert(first->prev==NULL);
        first->prev = bn;
        bn->next = first;
    }

    return bn;
}

bond_node* bond_node_remove(bond_node* first, bond_node* node){
    assert(node!=NULL);
    if(node->next==NULL && node->prev==NULL){
        bond_node_free(node);
        return NULL;
    }
    else if(node->prev==NULL){
        bond_node* next = node->next;
        next->prev = NULL;
        bond_node_free(node);
        return next;
    }
    else if(node->next==NULL){
        bond_node* prev = node->prev;
        prev->next = NULL;
        bond_node_free(node);
        return first;
    }
    else{
        bond_node* next = node->next;
        bond_node* prev = node->prev;
        next->prev = prev;
        prev->next = next;
        bond_node_free(node);
        return first;
    }
}

bond_node* bond_node_remove_all(bond_node* first){
    while(first!=NULL){
        first = bond_node_remove(first,first);
    }

    return NULL;
}

void bond_node_print_all(bond_node* first){
    bond_node* now;
    now = first;
    while(now!=NULL){
        printf("%d\n",now->graph_id);
        now = now->next;
    }
}

site_node* site_node_insertAfter(site_node* prev, double tau, int sigma, int spin_id, bond_node* bnode){
    site_node* snode = site_node_alloc();
    snode->tau = tau;
    snode->sigma[0] = sigma;
    snode->sigma[1] = sigma;
    snode->spin_id  = spin_id;
    snode->cont = bnode;

    if(prev==NULL){
        return snode;
    }
    else if(prev->next==NULL){
        prev->next  = snode;
        snode->prev = prev;
    }
    else{
        site_node* next = prev->next;
        prev->next = snode;
        next->prev = snode;
        snode->next = next;
        snode->prev = prev;
    }

    return snode;
}

site_node* site_node_insertBefore(site_node* next, double tau, int sigma, int spin_id, bond_node* bnode){
    site_node* snode = site_node_alloc();
    snode->tau = tau;
    snode->sigma[0] = sigma;
    snode->sigma[1] = sigma;
    snode->spin_id  = spin_id;
    snode->cont = bnode;

    if(next==NULL){
    }
    else if(next->prev==NULL){
        snode->next = next;
        next->prev = snode;
    }
    else{
        site_node* prev = next->prev;
        prev->next = snode;
        next->prev = snode;
        snode->next = next;
        snode->prev = prev;
    }
    return snode;
}

site_node* site_node_remove(site_node* first, site_node* node){
    assert(node!=NULL);
    if(node->prev==NULL && node->next==NULL){
        site_node_free(node);
        return NULL;
    }
    else if(node->prev==NULL){
        site_node* next = node->next;
        next->prev = NULL;
        site_node_free(node);
        return next;
    }
    else if(node->next==NULL){
        site_node* prev = node->prev;
        prev->next = NULL;
        site_node_free(node);
        return first;
    }
    else{
        site_node* prev = node->prev;
        site_node* next = node->next;
        prev->next = next;
        next->prev = prev;
        site_node_free(node);
        return first;
    }
}

site_node* site_node_remove_all(site_node* first){
    while(first!=NULL){
        first = site_node_remove(first,first);
    }

    return NULL;
}

void linked_list_test(){
    bond_node* first = NULL;
    first = bond_node_insertFirst(first,0,2);
    first = bond_node_insertFirst(first,0,2);
    first = bond_node_insertFirst(first,1,2);
    first = bond_node_insertFirst(first,1,2);
    first = bond_node_insertFirst(first,1,2);

    first = bond_node_remove_all(first);

    if(first==NULL) printf("NULL!\n");

    bond_node_print_all(first);

    site_node* sn = NULL;

    for(int i=0;i<100000;++i){
        sn = site_node_insertBefore(sn,0.9,1,0,NULL);
        sn = site_node_insertBefore(sn,0.8,1,0,NULL);
        sn = site_node_insertBefore(sn,0.7,1,0,NULL);
        sn = site_node_insertBefore(sn,0.6,1,0,NULL);
        sn = site_node_insertBefore(sn,0.5,1,0,NULL);
    }
    
    sn = site_node_remove_all(sn);
}
