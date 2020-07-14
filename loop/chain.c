#include "chain.h"

chain* chain_alloc(int size){
    chain* c = (chain*)malloc(sizeof(chain));
    if(c==NULL){
        printf("chain_alloc : memory allocating fail!\n");
        exit(1);
    }
    c->flag = 0;
    c->size = size;
    c->n = 0;
    c->node = (kink*)malloc(2*size*sizeof(kink));
    if(c->node==NULL){
        printf("chain_alloc : memory allocating fail!\n");
        exit(1);
    }

    for(int i=0;i<2*size;++i)
        (c->node)[i].key = UINT64_MAX;

    return c;
}

void chain_free(chain* c){
    free(c->node);
    free(c);
}

void chain_realloc(chain* c, int size){
    if(size>(c->size)){
        kink* ks = (kink*)malloc(2*size*sizeof(kink));
        if(ks==NULL){
            printf("chain_realloc : memory allocating fail!\n");
            exit(1);
        }

        for(int i=0;i<2*size;++i)
            ks[i].key = UINT64_MAX;

        for(int i=0;i<(c->n);++i)
            ks[(c->flag)*size+i] = (c->node)[(c->flag)*(c->size)+i];

        free(c->node);
        c->size = size;
        c->node = ks;
    }
}

void chain_insert(chain* c, double* tau, uint64_t* key, int n, int spin_id){
    assert((n+c->n)<=(c->size));

    int flag = c->flag;
    int size = c->size;

    int i=0;
    int j=0;
    int k=0;
    int m;

    int state = c->state;
    while((i<n) && (j<(c->n))){
        if((c->node[flag*size+j]).key==UINT64_MAX){
            ++j;
        }
        else if(tau[i] < (c->node[flag*size+j]).tau){
            m = (flag^1)*size+k;
            c->node[m].tau = tau[i];
            c->node[m].state[0] = state;
            c->node[m].state[1] = state;
            c->node[m].key = key[i];
            c->node[m].spin_id = spin_id;
            ++i;
            ++k;
        }
        else if(tau[i] > (c->node[flag*size+j]).tau){
            c->node[(flag^1)*size+k] = c->node[flag*size+j];
            state = c->node[flag*size+j].state[1];
            ++j;
            ++k;
        }
    }

    if(i==n){
        for(;j<(c->n);++j){
            c->node[(flag^1)*size+k] = c->node[flag*size+j];
            ++k;
        }
    }
    else if(j==(c->n)){
        for(;i<n;++i){
            m = (flag^1)*size+k;
            c->node[m].tau = tau[i];
            c->node[m].state[0] = state;
            c->node[m].state[1] = state;
            c->node[m].key = key[i];
            c->node[m].spin_id = spin_id;
            ++k;
        }
    }
    
    c->flag = flag^1;
    c->n = k;
}

void chain_print_state(chain* c){
    int flag = c->flag;
    int size = c->size;
    int n = c->n;
    int state = c->state;

    printf("---------------------------------\n");
    printf("flag = %d\n",flag);
    printf("size = %d\n",size);
    printf("n    = %d\n",n);
    printf("initial state = %d\n",state);
    printf("(tau,s0,s1,spin_id,key)\n");
    
    double tau;
    int s0,s1,spin_id;
    uint64_t key;
    for(int i=0;i<n;++i){
        tau = c->node[flag*size+i].tau;
        s0 = c->node[flag*size+i].state[0];
        s1 = c->node[flag*size+i].state[1];
        spin_id = c->node[flag*size+i].spin_id;
        key = c->node[flag*size+i].key;
        printf("%d (%.4f, %d, %d, %d, %lu)\n",i,tau,s0,s1,spin_id,key);
    }
    printf("---------------------------------\n");
}

int chain_test(){
    int size = 1000;
    chain* c = chain_alloc(size);
    c->state = 1;

    double tau1[5] = {1,3,5,7,9};
    uint64_t key1[5] = {1,2,3,4,5};

    double tau2[3] = {2,4,6};
    uint64_t key2[3] = {6,7,8};

    chain_insert(c,tau1,key1,5,0);
    //for(int i=0;i<c->n;++i)
    //    printf("%.4f\n",c->node[(c->flag)*(c->size)+i].tau);

    chain_insert(c,tau2,key2,3,1);
    for(int i=0;i<c->n;++i){
        double tau  = c->node[(c->flag)*(c->size)+i].tau;
        uint64_t key     = c->node[(c->flag)*(c->size)+i].key;
        int spin_id = c->node[(c->flag)*(c->size)+i].spin_id;
        printf("%.4f %lu %d\n",tau,key,spin_id);
    }
    printf("\n%d\n",c->n);

    chain_free(c);

    printf("%d\n",1<<10);
    
    return 0;
}
