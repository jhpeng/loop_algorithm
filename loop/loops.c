#include "loops.h"

void loops_update_table(table* t){
    int i,j,nspin,check;
    uint64_t key;

    int size = t->size;
    int n = t->n;

    for(i=0;i<size;++i){
        key = t->list[i].key;
        if(key!=UINT64_MAX){
            nspin = t->list[i].nspin;
            check = 1;
            for(j=0;j<nspin;++j){
                if((t->list[i].link_spin[j])==-2)
                    (t->list[i].state[j]) *= -1;

                if((t->list[i].link_spin[j+nspin])==-2)
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

void loops_update_chain(chain* c, table* t, gsl_rng* rng){
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
        }
    }

    if(n>0) c->state = c->node[0].state[0];
    else{
        if(gsl_rng_uniform_pos(rng)<0.5) c->state *=-1;
    }
}

void loops_link_vertex(chain* c, table* t){
    int size = c->size;
    int n    = c->n;
    int flag = c->flag;

    int id[2];
    int nspin;
    int spin_id[2];
    uint64_t key[2];
    item* it[2];
    if(n>0){
        for(int i=0;i<(n-1);++i){
            id[0] = flag*size+i;
            id[1] = flag*size+(i+1)%n;
            key[0] = c->node[id[0]].key;
            key[1] = c->node[id[1]].key;
            spin_id[0] = c->node[id[0]].spin_id;
            spin_id[1] = c->node[id[1]].spin_id;

            it[0] = table_search_from_key(t,key[0]);
            it[1] = table_search_from_key(t,key[1]);

            nspin = it[0]->nspin;

            it[0]->link_key[spin_id[0]+nspin] = key[1];
            it[1]->link_key[spin_id[1]] = key[0];
            it[0]->link_spin[spin_id[0]+nspin] = spin_id[1];
            it[1]->link_spin[spin_id[1]] = spin_id[0]+nspin;
        }
    }
}

// type : 1
// horizontal graph
static int rule_hori[4] = {1,0,3,2};

static void loops_looping(table* t, int item_id, int spin_id, gsl_rng* rng){
    assert(item_id<t->size);
    assert(t->list[item_id].key!=UINT64_MAX);
    if(t->list[item_id].link_spin[spin_id]>=0){
        uint64_t key_now,key_next;
        int spin_now,spin_next,type;

        int flag=-1;
        if(gsl_rng_uniform_pos(rng)<0.5) flag=-2;

        key_now = t->list[item_id].key;
        spin_now = spin_id;

        while(t->list[item_id].link_spin[spin_now]>=0){
            type = t->list[item_id].type;
            t->list[item_id].link_spin[spin_now] = flag;
            if(type==1){
                key_next  = t->list[item_id].link_key[rule_hori[spin_now]];
                spin_next = t->list[item_id].link_spin[rule_hori[spin_now]];
                t->list[item_id].link_spin[rule_hori[spin_now]] = flag;
            }
            //printf("%d %lu %lu\n",item_id,key_now,key_next);

            key_now  = key_next;
            spin_now = spin_next;
            item_id = table_hash(t,key_now);
        }
    }
}

void loops_traverse(table* t, gsl_rng* rng){
    int item_id,spin_id,nspin;
    int size = t->size;
    for(item_id=0;item_id<size;++item_id){
        if(t->list[item_id].key!=UINT64_MAX){
            nspin = t->list[item_id].nspin;
            for(spin_id=0;spin_id<2*nspin;spin_id+=2){
                loops_looping(t,item_id,spin_id,rng);
            }
        }
    }
}

#include "insert.h"
int main(){
    int nc = 2048;
    int scale = 16;
    double w = 1.0;
    double beta = 10;
    int seed = 21203;
    int nsweep = 100;
    int nx=8;

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    chain* c[nx];
    for(int i=0;i<nx;++i){
        c[i] = chain_alloc(nc);
        c[i]->state = 2*(i%2)-1;
    }

    table* t = table_alloc(scale);

    for(int i=0;i<nsweep;++i){
        for(int j=0;j<nx;++j)
            insert_horizontal_graph(c[j],c[(j+1)%nx],t,w,beta,rng);
        
        for(int j=0;j<nx;++j)
            loops_link_vertex(c[j],t);

        loops_traverse(t,rng);

        loops_update_table(t);

        for(int j=0;j<nx;++j)
            loops_update_chain(c[j],t,rng);

        for(int j=0;j<nx;++j)
            printf("%d ",c[j]->state);

        printf("\n");
            //chain_print_state(c[j]);
        //table_print_state(t);
    }
    

    gsl_rng_free(rng);
    for(int i=0;i<nx;++i)
        chain_free(c[i]);
    table_free(t);

    return 0;
}
