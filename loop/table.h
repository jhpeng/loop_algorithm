#ifndef table_h
#define table_h

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <assert.h>

#define NSPIN_MAX 4

extern uint64_t table_statistic_index;

typedef struct item{
    uint64_t key;
    int type;
    int nspin;
    int state[2*NSPIN_MAX];
    uint64_t link_key[2*NSPIN_MAX];
    int link_spin[2*NSPIN_MAX];
} item;

typedef struct table{
    int size;
    int n;
    uint64_t key;
    item* list;
} table;


table* table_alloc(int scale);

void table_free(table* t);

void table_realloc(table* t);

int table_hash(table* t, uint64_t key);

uint64_t table_generate_key(table* t);

item* table_search_from_key(table* t, uint64_t key);

void table_print_state(table* t);

#endif
