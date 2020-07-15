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
    for(int i=0;i<size;++i){
        t->list[i].key = UINT64_MAX;
        t->list[i].type  = 0;
        t->list[i].nspin = 0;
        for(int j=0;j<(2*NSPIN_MAX);++j){
            t->list[i].state[j] = 0;
            t->list[i].link_key[j] = UINT64_MAX;
            t->list[i].link_spin[j] = -1;
        }
    }

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

    for(int i=0;i<size;++i){
        it[i].key = UINT64_MAX;
        it[i].type  = 0;
        it[i].nspin = 0;
        for(int j=0;j<(2*NSPIN_MAX);++j){
            it[i].state[j] = 0;
            it[i].link_key[j] = UINT64_MAX;
            it[i].link_spin[j] = -1;
        }
    }

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

void table_print_state(table* t){
    int size = t->size;
    int n = t->n;
    uint64_t key = t->key;

    printf("---------------table--------------\n");
    printf("size = %d\n",size);
    printf("n    = %d\n",n);
    printf("new key = %" PRIu64 "\n",key);

    int type,nspin;
    printf("(key,flag,type,nspin)\n");
    for(int i=0;i<size;++i){
        if(t->list[i].key!=UINT64_MAX){
            key   = t->list[i].key;
            type  = t->list[i].type;
            nspin = t->list[i].nspin;
            printf("(%" PRIu64 ", %d, %d) ",key,type,nspin);

            for(int j=0;j<2*nspin;++j) 
                printf(" %" PRIu64 " ",t->list[i].link_key[j]);

            printf("\n");
        }
    }
    printf("----------------------------------\n");
}

int table_test(){
    int scale = 13;

    table* t = table_alloc(scale);
    table_realloc(t);

    for(int i=0;i<10;++i){
        uint64_t u64 = table_generate_key(t);

        printf("%" PRIu64 "\n",u64);
    }   

    table_free(t);

    return 0;
}
