#include "cluster.h"

void cluster_update_table(table* t){
    int i,j,nspin,check;
    uint64_t key;

    int size = t->size;
    int n = t->n;

    int nspin_max2 = 2*NSPIN_MAX;

    for(i=0;i<size;++i){
        key = t->list[i].key;
        if(key!=UINT64_MAX){
            nspin = t->list[i].nspin;
            check = 1;
            for(j=0;j<nspin;++j){
                if((t->list[i].link_spin[j])>=2*nspin_max2)
                    (t->list[i].state[j]) *= -1;

                if((t->list[i].link_spin[j+nspin])>=2*nspin_max2)
                    (t->list[i].state[j+nspin]) *= -1;
                    
                if((t->list[i].state[j])!=(t->list[i].state[j+nspin]))
                    check=0;

                t->list[i].link_spin[j] = -1;
                t->list[i].link_spin[j+nspin] = -1;
            }

            if(check){
                t->list[i].key = UINT64_MAX;
                n--;
            }
        }
    }

    t->n = n;
}

void claster_update_chain(chain* c, table* t){
    int i,j,spin_id,nspin;
    uint64_t key;

    int size = c->size;
    int n    = c->n;
    int flag = c->flag;
    item* it;

    for(i=0;i<n;++i){
        j = flag*size+i;
        key = c->node[j].key;
        spin_id = c->node[j].spin_id;
        it = table_search_from_key(t,key);

        nspin = it->nspin;
        c->node[j].state[0] = it->state[spin_id];
        c->node[j].state[1] = it->state[spin_id+nspin];

        if((it->key)==UINT64_MAX){
            c->node[j].key = UINT64_MAX;
            assert(c->node[j].state[0]==c->node[j].state[1]);
        }
    }

    if(n>0) c->state = c->node[flag*size].state[0];
}

void cluster_link_vertex(chain* c, table* t){
    int size = c->size;
    int n    = c->n;
    int flag = c->flag;

    int id[2];
    int nspin[2];
    int spin_id[2];
    uint64_t key[2];
    item* it[2];
    if(n>0){
        for(int i=0;i<n;++i){
            id[0] = flag*size+i;
            id[1] = flag*size+(i+1)%n;
            key[0] = c->node[id[0]].key;
            key[1] = c->node[id[1]].key;
            spin_id[0] = c->node[id[0]].spin_id;
            spin_id[1] = c->node[id[1]].spin_id;

            it[0] = table_search_from_key(t,key[0]);
            it[1] = table_search_from_key(t,key[1]);

            nspin[0] = it[0]->nspin;
            nspin[1] = it[1]->nspin;

            it[0]->link_key[spin_id[0]+nspin[0]] = key[1];
            it[1]->link_key[spin_id[1]] = key[0];
            it[0]->link_spin[spin_id[0]+nspin[0]] = spin_id[1];
            it[1]->link_spin[spin_id[1]] = spin_id[0]+nspin[0];
        }
    }
}

//type : 4
//triangular cut graph
static int cluster_twolink[16] = {0,1,2,3,
                                  0,1,2,3,
                                  0,1,2,3,
                                  0,1,2,3};

//type : 5
//triangular cut graph
static int cluster_tricut[64] = {0,1,2,4,5,6,-1,-1,
                                 0,1,2,4,5,6,-1,-1,
                                 0,1,2,4,5,6,-1,-1,
                                 3,-1,-1,-1,-1,-1,-1,-1,
                                 0,1,2,4,5,6,-1,-1,
                                 0,1,2,4,5,6,-1,-1,
                                 0,1,2,4,5,6,-1,-1,
                                 7,-1,-1,-1,-1,-1,-1,-1};

//type : 6
//triangular graph
static int cluster_triang[36] = {0,1,2,3,4,5,
                                 0,1,2,3,4,5,
                                 0,1,2,3,4,5,
                                 0,1,2,3,4,5,
                                 0,1,2,3,4,5,
                                 0,1,2,3,4,5};

static uint64_t* cluster_key;
static int* cluster_spin;
static int cluster_size=0;

static void cluster_clustering(table* t, int item_id, int spin_id, gsl_rng* rng){
    assert(item_id<t->size);
    assert(t->list[item_id].key!=UINT64_MAX);

    if(cluster_size==0){
        cluster_size = t->size*NSPIN_MAX*2;
        cluster_key  = (uint64_t*)malloc(sizeof(uint64_t)*cluster_size);
        cluster_spin = (int*)malloc(sizeof(int)*cluster_size);
    }
    else if(cluster_size<t->size*NSPIN_MAX*2){
        cluster_size = t->size*NSPIN_MAX*2;
        cluster_key  = (uint64_t*)realloc(cluster_key,sizeof(uint64_t)*cluster_size);
        cluster_spin = (int*)realloc(cluster_spin,sizeof(int)*cluster_size);
    }

    int nspin_max2 = 2*NSPIN_MAX;

    if(t->list[item_id].link_spin[spin_id]<nspin_max2){
        uint64_t key_now = t->list[item_id].key;
        int spin_now = spin_id;
        int spin_next;
        int type,nspin,i;
        int* cluster;

        int flag=1;
        if(gsl_rng_uniform_pos(rng)<0.5) flag=2;

        int n=0;
        int m=0;
        while(m<=n){
            type  = t->list[item_id].type;
            nspin = t->list[item_id].nspin;


            if(type==4) cluster = cluster_twolink;
            else if(type==5) cluster = cluster_tricut;
            else if(type==6) cluster = cluster_triang;
            else{
                printf("cluster_clustering : no this type %d !\n",type);
                exit(1);
            }

            for(i=0;i<2*nspin;++i){
                spin_id = cluster[2*nspin*spin_now+i];
                if(spin_id!=-1){
                    cluster_key[n]  = t->list[item_id].key;
                    cluster_spin[n] = spin_id;
                    ++n;

                    t->list[item_id].link_spin[spin_id] += nspin_max2*flag;
                }
            }

            int check=1;
            while(m<n){
                item_id = table_hash(t,cluster_key[m]);
                key_now = t->list[item_id].link_key[cluster_spin[m]];
                spin_now = t->list[item_id].link_spin[cluster_spin[m]];
                spin_now -= nspin_max2*flag;
                ++m;


                item_id = table_hash(t,key_now);
                spin_next = t->list[item_id].link_spin[spin_now];
                if(spin_next<nspin_max2){
                    check=0;
                    break;
                }
            }

            if(check)
                break;
        }
    }
}

void cluster_traverse(table* t, gsl_rng* rng){
    int item_id,spin_id,nspin;
    int size = t->size;

    for(item_id=0;item_id<size;++item_id){
        if(t->list[item_id].key!=UINT64_MAX){
            nspin = t->list[item_id].nspin;
            for(spin_id=0;spin_id<2*nspin;++spin_id){
                if((t->list[item_id].link_spin[spin_id]<2*NSPIN_MAX) 
                    && (t->list[item_id].link_spin[spin_id]>=0))
                    cluster_clustering(t,item_id,spin_id,rng);
            }
        }
    }
}
