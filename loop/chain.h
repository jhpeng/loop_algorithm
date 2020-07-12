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

#endif
