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

void loops_update_chain(chain* c, table* t){
    int i=0;
    int j,spin_id,nspin;
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
}

#include "insert.h"
int main(){
    int nc = 2048;
    int scale = 16;
    double w = 1.0;
    double beta = 10;
    int seed = 21203;
    int nsweep = 1000;

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    chain* c1 = chain_alloc(nc);
    chain* c2 = chain_alloc(nc);
    table* t = table_alloc(scale);

    c1->state = 1;
    c2->state = -1;


    for(int i=0;i<nsweep;++i){
        insert_horizontal_graph(c1,c2,t,w,beta,rng);
        
        loops_update_table(t);
        loops_update_chain(c1,t);
        loops_update_chain(c2,t);

        //chain_print_state(c1);
        //chain_print_state(c2);
        //table_print_state(t);
    }
    

    gsl_rng_free(rng);
    chain_free(c1);
    chain_free(c2);
    table_free(t);

    return 0;
}
