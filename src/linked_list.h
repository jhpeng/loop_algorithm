#ifndef linked_list_h
#define linked_list_h

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

typedef struct bond_nodeT{
    int graph_id;               //graph_id for specific bond
    int nspin;                  //n-spin interaction in quantum system
    int* leg;                   //index of linked vertex  (size=2*nspin)
    int* sigma;                 //sigma for each spin     (size=2*nspin)
    int termination;            //temination signal for site_node (0->keep) (1->teminate)
    struct bond_nodeT* next;    //pointer to next bond_node
    struct bond_nodeT* prev;    //pointer to previous bond_node
} bond_node;

typedef struct site_nodeT{
    double tau;                 //position on the imaginary time
    int sigma[2];               //sigma before and after kink
    int spin_id;                //the spin_id for bond_node
    struct bond_nodeT* cont;    //pointer to the bond_node
    struct site_nodeT* next;    //pointer to next site_node
    struct site_nodeT* prev;    //pointer to previous site_node
} site_node;


//operator for bond_node
bond_node* bond_node_alloc(
                int nspin);

void bond_node_free(
                bond_node* bn);

bond_node*  bond_node_insertFirst(
                bond_node* first, 
                int graph_id, 
                int nspin);

bond_node* bond_node_remove(
                bond_node* first, 
                bond_node* node);

bond_node* bond_node_remove_all(
                bond_node* first);

void bond_node_print_all(
                bond_node* first);


//operator for site_node
site_node* site_node_alloc();

void site_node_free(
                site_node* sn);

site_node* site_node_insertAfter(
                site_node* first, 
                site_node* prev, 
                double tau, 
                int sigma, 
                int spin_id, 
                bond_node* bnode);

site_node* site_node_insertBefore(
                site_node* first, 
                site_node* next, 
                double tau, 
                int sigma, 
                int spin_id, 
                bond_node* bnode);

site_node* site_node_remove(
                site_node* first, 
                site_node* node);

site_node* site_node_remove_all(
                site_node* first);
#endif
