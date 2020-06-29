#ifndef configuration_h 
#define configuration_h 

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>
//#define NDEBUG
#include <assert.h>

#include <gsl/gsl_sort.h>

typedef struct kinks{
    size_t size;
    size_t nkink;
    int sigma_i;
    int* active;
    int* sigma_b;
    int* sigma_a;
    int* bond_id;
    int* type_id;
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

void kinks_flip_sigma_b(
                kinks* ks, 
                int kink_id);

void kinks_flip_sigma_a(
                kinks* ks, 
                int kink_id);

int kinks_sigma_from_tau(
                const kinks* ks, 
                double tau);

int kinks_check_available_id(
                kinks* ks);

int kinks_check_no_kink(
                kinks* ks, 
                int kink_id);

size_t kinks_get_size(kinks* ks);

size_t kinks_get_nkink(kinks* ks);

int kinks_get_sigma_i(kinks* ks);

int kinks_get_active(kinks* ks, size_t kink_id);

int kinks_get_sigma_b(kinks* ks, size_t kink_id);

int kinks_get_sigma_a(kinks* ks, size_t kink_id);

int kinks_get_bond_id(kinks* ks, size_t kink_id);

int kinks_get_type_id(kinks* ks, size_t kink_id);

double kinks_get_tau(kinks* ks, size_t kink_id);

size_t kinks_get_sort(kinks* ks, size_t rank);

int kinks_insert(
                kinks* ks, 
                int bond_id, 
                int type_id, 
                int sigma, 
                double tau);

void kinks_remove(
                kinks* ks, 
                int kink_id);

void kinks_sort_index_with_tau(
                kinks* ks);

void kinks_print_state(kinks* ks);
#endif
