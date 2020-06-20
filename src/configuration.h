#ifndef configuration_h 
#define configuration_h 

#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <gsl/gsl_sort.h>

typedef struct kinks{
    size_t size;
    size_t nkink;
    int sigma_i;
    int* active;
    int* sigma_b;
    int* sigma_a;
    int* bond_id;
    int* graph_id;
    double* tau;
    size_t* sort;
} kinks;

kinks* kinks_alloc(
                int size);

void kinks_free(
                kinks* ks);

void kinks_memcpy(
                kinks* dest, 
                const kinks* src);

void kinks_set_sigma_i(
                kinks* ks, 
                int sigma_i);

int kinks_sigma_from_tau(
                const kinks* ks, 
                double tau);

int kinks_insert(
                kinks* ks, 
                int bond_id, 
                int graph_id, 
                int sigma, 
                double tau);

void kinks_remove(
                kinks* ks, 
                int kink_id);

void kinks_sort_index_with_tau(
                kinks* ks);

#endif
