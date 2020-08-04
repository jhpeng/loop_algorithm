#include "insert.h"

static double* insert_tau=NULL;
static int insert_size=0;

static void generate_uniform_dist(double* tau, int *n, int size, double w, double beta, gsl_rng* rng){
    double dis = gsl_rng_uniform_pos(rng);
    double temp=-log(dis)/w;
    int m=0;

    while(temp<beta && m<size){
        tau[m] = temp;
        dis = gsl_rng_uniform_pos(rng);
        temp = temp - log(dis)/w;

        ++m;
    }

    *n = m;
}

void insert_horizontal_graph(chain* c1, chain* c2, table* t, double w, double beta, gsl_rng* rng){
    if(insert_tau==NULL){
        int size = 10000;
        insert_tau = (double*)malloc(sizeof(double)*size);
        insert_size = size;

        if(insert_tau==NULL){
            printf("insert_horizontal_graph : allocating memory fail!\n");
            exit(1);
        }
    }

    // generate uniform distribution to insert_tau
    int n=0;
    generate_uniform_dist(insert_tau,&n,insert_size,w,beta,rng);
    if(n==insert_size){
        insert_tau = (double*)realloc(insert_tau,sizeof(double)*2*insert_size);
        insert_size *= 2;

        if(insert_tau==NULL){
            printf("insert_horizontal_graph : allocating memory fail!\n");
            exit(1);
        }

        generate_uniform_dist(insert_tau,&n,insert_size,w,beta,rng);
    }

    // create the workong space
    double* stau = (double*)malloc(sizeof(double)*n);
    int* ss1 = (int*)malloc(sizeof(int)*n);
    int* ss2 = (int*)malloc(sizeof(int)*n);
    int n1 = c1->n+1;
    int n2 = c2->n+1;
    int* s1 = (int*)malloc(sizeof(int)*n1);
    int* s2 = (int*)malloc(sizeof(int)*n2);
    double* tau1 = (double*)malloc(sizeof(double)*n1);
    double* tau2 = (double*)malloc(sizeof(double)*n2);

    for(int ii=0;ii<(n1-1);++ii){
        int flag = c1->flag;
        int size = c1->size;
        s1[ii] = (c1->node[flag*size+ii]).state[0];
        tau1[ii] = (c1->node[flag*size+ii]).tau;
    }
    for(int ii=0;ii<(n2-1);++ii){
        int flag = c2->flag;
        int size = c2->size;
        s2[ii] = (c2->node[flag*size+ii]).state[0];
        tau2[ii] = (c2->node[flag*size+ii]).tau;
    }
    s1[n1-1] = c1->state;
    s2[n2-1] = c2->state;
    tau1[n1-1] = beta;
    tau2[n2-1] = beta;

    // check the condition for inserting the graph
    double tau;
    int i=0;
    int j=0;
    int m=0;
    int k;
    for(k=0;k<n;++k){
        tau = insert_tau[k];
        while(tau1[i]<tau) ++i;
        while(tau2[j]<tau) ++j;

        if((s1[i]!=s2[j] && tau1[i]!=tau) && tau2[j]!=tau){
            stau[m] = tau;
            ss1[m] = s1[i];
            ss2[m] = s2[j];
            ++m;
        }
    }

    // check memory space and buffer
    while(2*(m+t->n)>(t->size)) table_realloc(t);
    if((m+c1->n)>(c1->size)) chain_realloc(c1,(m+c1->n));
    if((m+c2->n)>(c2->size)) chain_realloc(c2,(m+c2->n));

    uint64_t* key = (uint64_t*)malloc(sizeof(uint64_t)*m);
    for(k=0;k<m;++k){
        key[k] = table_generate_key(t);
        i = table_hash(t,key[k]);

        (t->list[i]).key = key[k];
        (t->list[i]).type  = 1;
        (t->list[i]).nspin = 2;

        (t->list[i]).state[0] = ss1[k];
        (t->list[i]).state[1] = ss2[k];
        (t->list[i]).state[2] = ss1[k];
        (t->list[i]).state[3] = ss2[k];
    }

    chain_insert(c1,stau,key,m,0);
    chain_insert(c2,stau,key,m,1);
    t->n += m;

    // free the working space
    free(stau);
    free(ss1);
    free(ss2);
    free(s1);
    free(s2);
    free(tau1);
    free(tau2);
    free(key);
}

int test_insert(){
    int nc = 2048;
    int scale = 16;
    double w = 1.0;
    double beta = 10;
    int seed = 21203;

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    chain* c1 = chain_alloc(nc);
    chain* c2 = chain_alloc(nc);
    table* t = table_alloc(scale);

    c1->state = 1;
    c2->state = -1;


    insert_horizontal_graph(c1,c2,t,w,beta,rng);
    insert_horizontal_graph(c1,c2,t,w,beta,rng);
    insert_horizontal_graph(c1,c2,t,w,beta,rng);
    insert_horizontal_graph(c1,c2,t,w,beta,rng);

    chain_print_state(c1);
    chain_print_state(c2);
    table_print_state(t);
    

    gsl_rng_free(rng);
    chain_free(c1);
    chain_free(c2);
    table_free(t);

    free(insert_tau);

    return 0;
}
