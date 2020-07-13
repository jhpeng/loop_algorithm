#include "table.h"

table* table_alloc(int scale){
    int size = 1<<scale;
    
    table* t = (table*)malloc(sizeof(table));
    if(t==NULL){
        printf("table_alloc : memory allocating fail!\n");
        exit(1);
    }
    t->list = (item*)malloc(sizeof(item)*size);
    if(t->list==NULL){
        printf("table_alloc : memory allocating fail!\n");
        exit(1);
    }
    
    t->size = size;
    t->n = 0;
    t->key = 0;
    for(int i=0;i<size;++i)
        t->list[i].key = UINT64_MAX;

    return t;
}

void table_free(table* t){
    free(t->list);
    free(t);
}

void table_realloc(table* t){
    int size = 2*(t->size);

    item* it = (item*)malloc(sizeof(item)*size);
    if(it==NULL){
        printf("table_realloc : memory allocating fail!\n");
        exit(1);
    }

    for(int i=0;i<size;++i)
        it[i].key = UINT64_MAX;

    for(int i=0;i<(t->size);++i)
        it[i] = t->list[i];

    free(t->list);
    t->list = it;
    t->size = size;
}

int table_hash(table* t, uint64_t key){
    return (int)(key&(t->size-1));
}

uint64_t table_generate_key(table* t){
    int i,check=1;
    while(check){
        ++(t->key);
        i = table_hash(t,t->key);
        if(t->list[i].key==UINT64_MAX) check=0;
    }

    return t->key;
}

item* table_search_from_key(table* t, uint64_t key){
    int i = table_hash(t,key);

    return &(t->list[i]);
}

int table_test(){
    int scale = 13;

    table* t = table_alloc(scale);
    table_realloc(t);

    for(int i=0;i<10;++i){
        uint64_t u64 = table_generate_key(t);

        printf("%llu\n",u64);
    }   

    table_free(t);

    return 0;
}
