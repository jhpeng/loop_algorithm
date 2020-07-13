#ifndef table_h
#define table_h

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#define NSPIN_MAX 8

typedef struct item{
    uint64_t key;
    int flag;
    int type;
    int nspin;
    int state[NSPIN_MAX];
    uint64_t link_key[NSPIN_MAX];
    int link_spin[NSPIN_MAX];
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

#endif
