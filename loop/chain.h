#ifndef chain_h
#define chain_h

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

typedef struct kink{
    double tau;
    int state[2];
    int spin_id;
    uint64_t key;
} kink;

typedef struct chain{
    int flag;
    int size;
    int n;
    kink* node;
    int state;
} chain;

chain* chain_alloc(int size);

void chain_free(chain* c);

void chain_realloc(chain* c, int size);

void chain_insert(
                chain* c, 
                double* tau, 
                uint64_t* key, 
                int n, 
                int spin_id);
#endif
