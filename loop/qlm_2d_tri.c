#include "qlm_2d_tri.h"

static void gauss_law_A(chain** c, table* t, model* m, int x, int y){
    int ai,aj,ak,bi,bj,bk;
    int block_id;
    double beta = m->beta;

    if(m->nsite!=2*x*y){
        printf("gauss_law_A : nsite should be 2*x*y!\n");
        exit(1);
    }

    int ix,iy;
    int sb[3];
    uint64_t key;
    item* it;
    for(block_id=0;block_id<x*y;++block_id){
        ix = block_id%x;
        iy = block_id/x;

        ai = iy*x+ix;
        aj = ((iy+1)%y)*x+ix;
        ak = ((iy+1)%y)*x+(ix+1)%x;

        bi = iy*x+ix               + x*y;
        bj = iy*x+(ix+1)%x         + x*y;
        bk = ((iy+1)%y)*x+(ix+1)%x + x*y;

        sb[0] = c[bi]->state;
        sb[1] = c[bj]->state;
        sb[2] = c[bk]->state;

        if(!((sb[0]==sb[1]) && (sb[1]==sb[2]))){
            if(sb[0]==sb[1]){
                key = table_generate_key(t);
                chain_insert(c[aj],&beta,&key,1,0);
                chain_insert(c[ak],&beta,&key,1,1);

                it = table_search_from_key(t,key);
                it->key   = key;
                it->type  = 4;
                it->nspin = 2;
                it->state[0] = c[aj]->state;
                it->state[1] = c[aj]->state;
                it->state[2] = c[ak]->state;
                it->state[3] = c[ak]->state;

                for(int j=0;j<4;++j){
                    it->link_key[j]  = UINT64_MAX;
                    it->link_spin[j] = -1;
                }

                ++t->n;
            }
            else if(sb[1]==sb[2]){
                key = table_generate_key(t);
                chain_insert(c[ai],&beta,&key,1,0);
                chain_insert(c[aj],&beta,&key,1,1);

                it = table_search_from_key(t,key);
                it->key   = key;
                it->type  = 4;
                it->nspin = 2;
                it->state[0] = c[ai]->state;
                it->state[1] = c[ai]->state;
                it->state[2] = c[aj]->state;
                it->state[3] = c[aj]->state;

                for(int j=0;j<4;++j){
                    it->link_key[j]  = UINT64_MAX;
                    it->link_spin[j] = -1;
                }

                ++t->n;
            }
            else if(sb[2]==sb[0]){
                key = table_generate_key(t);
                chain_insert(c[ai],&beta,&key,1,0);
                chain_insert(c[ak],&beta,&key,1,1);

                it = table_search_from_key(t,key);
                it->key   = key;
                it->type  = 4;
                it->nspin = 2;
                it->state[0] = c[ai]->state;
                it->state[1] = c[ai]->state;
                it->state[2] = c[ak]->state;
                it->state[3] = c[ak]->state;

                for(int j=0;j<4;++j){
                    it->link_key[j]  = UINT64_MAX;
                    it->link_spin[j] = -1;
                }

                ++t->n;
            }
        }
    }
}

static void gauss_law_B(chain** c, table* t, model* m, int x, int y){
    int ai,aj,ak,bi,bj,bk;
    int block_id;
    double beta = m->beta;

    if(m->nsite!=2*x*y){
        printf("gauss_law_B : nsite should be 2*x*y!\n");
        exit(1);
    }

    int ix,iy;
    int sa[3];
    uint64_t key;
    item* it;
    for(block_id=0;block_id<x*y;++block_id){
        ix = block_id%x;
        iy = block_id/x;

        ai = iy*x+ix;
        aj = ((iy+1)%y)*x+ix;
        ak = ((iy+1)%y)*x+(ix+1)%x;

        bi = iy*x+ix               + x*y;
        bj = iy*x+(ix+1)%x         + x*y;
        bk = ((iy+1)%y)*x+(ix+1)%x + x*y;

        sa[0] = c[ai]->state;
        sa[1] = c[aj]->state;
        sa[2] = c[ak]->state;

        if(!((sa[0]==sa[1]) && (sa[1]==sa[2]))){
            if(sa[0]==sa[1]){
                key = table_generate_key(t);
                chain_insert(c[bj],&beta,&key,1,0);
                chain_insert(c[bk],&beta,&key,1,1);

                it = table_search_from_key(t,key);
                it->key   = key;
                it->type  = 4;
                it->nspin = 2;
                it->state[0] = c[bj]->state;
                it->state[1] = c[bj]->state;
                it->state[2] = c[bk]->state;
                it->state[3] = c[bk]->state;

                for(int j=0;j<4;++j){
                    it->link_key[j]  = UINT64_MAX;
                    it->link_spin[j] = -1;
                }

                ++t->n;
            }
            else if(sa[1]==sa[2]){
                key = table_generate_key(t);
                chain_insert(c[bi],&beta,&key,1,0);
                chain_insert(c[bj],&beta,&key,1,1);

                it = table_search_from_key(t,key);
                it->key   = key;
                it->type  = 4;
                it->nspin = 2;
                it->state[0] = c[bi]->state;
                it->state[1] = c[bi]->state;
                it->state[2] = c[bj]->state;
                it->state[3] = c[bj]->state;

                for(int j=0;j<4;++j){
                    it->link_key[j]  = UINT64_MAX;
                    it->link_spin[j] = -1;
                }

                ++t->n;
            }
            else if(sa[2]==sa[0]){
                key = table_generate_key(t);
                chain_insert(c[bi],&beta,&key,1,0);
                chain_insert(c[bk],&beta,&key,1,1);

                it = table_search_from_key(t,key);
                it->key   = key;
                it->type  = 4;
                it->nspin = 2;
                it->state[0] = c[bi]->state;
                it->state[1] = c[bi]->state;
                it->state[2] = c[bk]->state;
                it->state[3] = c[bk]->state;

                for(int j=0;j<4;++j){
                    it->link_key[j]  = UINT64_MAX;
                    it->link_spin[j] = -1;
                }

                ++t->n;
            }
        }
    }
}

static void update_A(chain** c, table* t, model* m, int x, int y, gsl_rng* rng){
    int i,j,k,l,bond_id;
    int nsite = m->nsite;
    double beta = m->beta;
    double weight;

    int nsq = x*y;

    chain* c_temp[4];
    for(bond_id=nsq;bond_id<2*nsq;++bond_id){
        //int type = m->type[bond_id];
        i = m->bond2site[bond_id*NSPIN_MAX+0];
        j = m->bond2site[bond_id*NSPIN_MAX+1];
        k = m->bond2site[bond_id*NSPIN_MAX+2];
        l = m->bond2site[bond_id*NSPIN_MAX+3];
        weight = m->weight[bond_id];


        c_temp[0] = c[i];
        c_temp[1] = c[j];
        c_temp[2] = c[k];
        c_temp[3] = c[l];

        insert_triangular_cut_graph(c_temp,t,weight,beta,rng);
    }
    for(bond_id=2*nsq;bond_id<3*nsq;++bond_id){
        //int type = m->type[bond_id];
        i = m->bond2site[bond_id*NSPIN_MAX+0];
        j = m->bond2site[bond_id*NSPIN_MAX+1];
        k = m->bond2site[bond_id*NSPIN_MAX+2];
        weight = m->weight[bond_id];


        c_temp[0] = c[i];
        c_temp[1] = c[j];
        c_temp[2] = c[k];

        insert_triangular_graph(c_temp,t,weight,beta,rng);
    }

    gauss_law_A(c,t,m,x,y);

    for(i=0;i<nsq;++i)
        cluster_link_vertex(c[i],t);

    cluster_traverse(t,rng);

    cluster_update_table(t);

    for(i=0;i<nsite;++i)
        claster_update_chain(c[i],t);
}

static void update_B(chain** c, table* t, model* m, int x, int y, gsl_rng* rng){
    int i,j,k,l,bond_id;
    int nsite = m->nsite;
    double beta = m->beta;
    double weight;

    int nsq = nsite/2;

    chain* c_temp[4];
    for(bond_id=0;bond_id<nsq;++bond_id){
        //int type = m->type[bond_id];
        i = m->bond2site[bond_id*NSPIN_MAX+0];
        j = m->bond2site[bond_id*NSPIN_MAX+1];
        k = m->bond2site[bond_id*NSPIN_MAX+2];
        l = m->bond2site[bond_id*NSPIN_MAX+3];
        weight = m->weight[bond_id];


        c_temp[0] = c[i];
        c_temp[1] = c[j];
        c_temp[2] = c[k];
        c_temp[3] = c[l];

        insert_triangular_cut_graph(c_temp,t,weight,beta,rng);
    }
    for(bond_id=3*nsq;bond_id<4*nsq;++bond_id){
        //int type = m->type[bond_id];
        i = m->bond2site[bond_id*NSPIN_MAX+0];
        j = m->bond2site[bond_id*NSPIN_MAX+1];
        k = m->bond2site[bond_id*NSPIN_MAX+2];
        weight = m->weight[bond_id];


        c_temp[0] = c[i];
        c_temp[1] = c[j];
        c_temp[2] = c[k];

        insert_triangular_graph(c_temp,t,weight,beta,rng);
    }

    gauss_law_B(c,t,m,x,y);

    for(i=nsq;i<2*nsq;++i)
        cluster_link_vertex(c[i],t);

    cluster_traverse(t,rng);

    cluster_update_table(t);

    for(i=0;i<nsite;++i)
        claster_update_chain(c[i],t);
}

model* generate_QLM_2d_triangular(int x, int y, double beta, double lambda){
    int nsite = 2*x*y;
    int nbond = 4*x*y;

    double w1 = 1.0;
    double w2 = lambda;

    int* type = (int*)malloc(sizeof(int)*nbond);
    int* bond2site = (int*)malloc(sizeof(int)*nbond*NSPIN_MAX);
    double* weight = (double*)malloc(sizeof(double)*nbond);

    int i,j,k,ix,iy,bond_id;
    for(bond_id=0;bond_id<x*y;++bond_id){
        type[bond_id] = 5;
        weight[bond_id] = w1;

        ix = bond_id%x;
        iy = bond_id/x;

        i = iy*x+ix;
        j = iy*x+(ix+1)%x;
        k = ((iy+1)%y)*x+(ix+1)%x;

        bond2site[bond_id*NSPIN_MAX+0] = i;
        bond2site[bond_id*NSPIN_MAX+1] = j;
        bond2site[bond_id*NSPIN_MAX+2] = k;
        bond2site[bond_id*NSPIN_MAX+3] = j+x*y;
    }
    for(bond_id=0;bond_id<x*y;++bond_id){
        type[bond_id+x*y] = 5;
        weight[bond_id+x*y] = w1;

        ix = bond_id%x;
        iy = bond_id/x;

        i = iy*x+ix;
        j = ((iy+1)%y)*x+ix;
        k = ((iy+1)%y)*x+(ix+1)%x;

        bond2site[(bond_id+x*y)*NSPIN_MAX+0] = i+x*y;
        bond2site[(bond_id+x*y)*NSPIN_MAX+1] = j+x*y;
        bond2site[(bond_id+x*y)*NSPIN_MAX+2] = k+x*y;
        bond2site[(bond_id+x*y)*NSPIN_MAX+3] = j;
    }
    for(bond_id=0;bond_id<x*y;++bond_id){
        type[bond_id+2*x*y] = 6;
        weight[bond_id+2*x*y] = w2;

        ix = bond_id%x;
        iy = bond_id/x;

        i = iy*x+ix;
        j = iy*x+(ix+1)%x;
        k = ((iy+1)%y)*x+(ix+1)%x;

        bond2site[(bond_id+2*x*y)*NSPIN_MAX+0] = i;
        bond2site[(bond_id+2*x*y)*NSPIN_MAX+1] = j;
        bond2site[(bond_id+2*x*y)*NSPIN_MAX+2] = k;
    }
    for(bond_id=0;bond_id<x*y;++bond_id){
        type[bond_id+3*x*y] = 6;
        weight[bond_id+3*x*y] = w2;

        ix = bond_id%x;
        iy = bond_id/x;

        i = iy*x+ix;
        j = ((iy+1)%y)*x+ix;
        k = ((iy+1)%y)*x+(ix+1)%x;

        bond2site[(bond_id+3*x*y)*NSPIN_MAX+0] = i+x*y;
        bond2site[(bond_id+3*x*y)*NSPIN_MAX+1] = j+x*y;
        bond2site[(bond_id+3*x*y)*NSPIN_MAX+2] = k+x*y;
    }


    model* m = (model*)malloc(sizeof(model));

    m->nsite = nsite;
    m->nbond = nbond;
    m->beta = beta;
    m->type = type;
    m->bond2site = bond2site;
    m->weight = weight;

    return m;
}

#include <string.h>

int main(int argc, char** argv){
    int x = atoi(argv[1]);
    int y = atoi(argv[2]);
    double lambda = atof(argv[3]);
    double beta = atof(argv[4]);
    int ntherm = atoi(argv[5]);
    int nsweep = atoi(argv[6]);
    int seed = atoi(argv[7]);

    model* m = generate_QLM_2d_triangular(x,y,beta,lambda);

    if(0){
    int site[4];
    for(int i=0;i<m->nbond/2;++i){
        site[0] = m->bond2site[i*NSPIN_MAX+0];
        site[1] = m->bond2site[i*NSPIN_MAX+1];
        site[2] = m->bond2site[i*NSPIN_MAX+2];
        site[3] = m->bond2site[i*NSPIN_MAX+3];

        printf("bond_id=%d \t",i);
        for(int j=0;j<4;++j){
            if(!(site[j]<x*y)){
                site[j] -= x*y;
                printf("B %d  ",site[j]);
            }
            else{
                printf("A %d  ",site[j]);
            }
        }
        int type = m->type[i];
        double weight = m->weight[i];
        printf("%d %.3f",type,weight);
        printf("\n");
    }
    for(int i=m->nbond/2;i<m->nbond;++i){
        site[0] = m->bond2site[i*NSPIN_MAX+0];
        site[1] = m->bond2site[i*NSPIN_MAX+1];
        site[2] = m->bond2site[i*NSPIN_MAX+2];

        printf("bond_id=%d \t",i);
        for(int j=0;j<3;++j){
            if(!(site[j]<x*y)){
                site[j] -= x*y;
                printf("B %d  ",site[j]);
            }
            else{
                printf("A %d  ",site[j]);
            }
        }
        int type = m->type[i];
        double weight = m->weight[i];
        printf("%d %.3f",type,weight);
        printf("\n");
    }
    }

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    chain* c[m->nsite];
    for(int i=0;i<m->nsite;++i){
       c[i] = chain_alloc(2000);
       c[i]->state = 1;
    }

    table* t = table_alloc(20);

    for(int i=0;i<ntherm;++i){
        update_A(c,t,m,x,y,rng);
        update_B(c,t,m,x,y,rng);

        for(int j=0;j<m->nsite;++j){
            int n = c[j]->n;
            for(int k=0;n>k;++k){
                int size = c[j]->size;
                int flag = c[j]->flag;
                int s0 = c[j]->node[flag*size+(k+1)%n].state[0];
                int s1 = c[j]->node[flag*size+k].state[1];

                if(s0!=s1){
                    printf("Vailoate the periodic boundary condition!\n");
                    exit(1);
                }
            }

            //chain_print_state(c[j]);
        }
    }

    for(int i=0;i<nsweep;++i){
        update_A(c,t,m,x,y,rng);
        update_B(c,t,m,x,y,rng);

        int Ma=0;
        int Mb=0;
        for(int iy=0;iy<y;++iy){
            for(int ix=0;ix<x;++ix){
                Ma += c[iy*x+ix]->state;
            }
            for(int ix=0;ix<x;++ix){
                Mb += c[iy*x+ix+x*y]->state;
            }
        }

        printf("%d %d\n",Ma,Mb);
    }

    return 0;
}
