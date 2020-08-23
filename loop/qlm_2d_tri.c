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

#if 0
    int s0,s1,nspin,spin_id,size,flag;
    uint64_t key;
    item* it;
    chain* ct;
    for(i=0;i<nsq;++i){
        ct = c[i+nsq];
        for(j=0;j<ct->n;++j){
            size = ct->size;
            flag = ct->flag;
            key     = ct->node[size*flag+j].key;
            spin_id = ct->node[size*flag+j].spin_id;

            it = table_search_from_key(t,key);
            nspin = it->nspin;

            s0 = it->link_spin[spin_id      ];
            s1 = it->link_spin[spin_id+nspin];

            if(s0!=-1 || s1!=-1){
                printf("update_A : someting wrong %d %d\n",s0,s1);
                exit(1);
            }
        }
    }
#endif

    cluster_traverse(t,rng);

    cluster_update_table(t);

    //for(i=0;i<nsq;++i)
    for(i=0;i<2*nsq;++i)
        claster_update_chain(c[i],t);
}

static void update_B(chain** c, table* t, model* m, int x, int y, gsl_rng* rng){
    int i,j,k,l,bond_id;
    double beta = m->beta;
    double weight;

    int nsq = x*y;

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

#if 0
    int s0,s1,nspin,spin_id,size,flag;
    uint64_t key;
    item* it;
    chain* ct;
    for(i=0;i<nsq;++i){
        ct = c[i];
        for(j=0;j<ct->n;++j){
            size = ct->size;
            flag = ct->flag;
            key     = ct->node[size*flag+j].key;
            spin_id = ct->node[size*flag+j].spin_id;

            it = table_search_from_key(t,key);
            nspin = it->nspin;

            s0 = it->link_spin[spin_id      ];
            s1 = it->link_spin[spin_id+nspin];

            if(s0!=-1 || s1!=-1){
                printf("update_B : someting wrong %d %d\n",s0,s1);
                exit(1);
            }
        }
    }
#endif

    cluster_traverse(t,rng);

    cluster_update_table(t);

    //for(i=nsq;i<2*nsq;++i)
    for(i=0;i<2*nsq;++i)
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

#include <gsl/gsl_sort.h>
#include <string.h>
double* mtau;
size_t* msite;
size_t* msort;
int msize=0;
void qlm_measurement(chain** c, table* t, model* m, int x, int y, double lambda, int seed){
    int xy = m->nsite/2;
    int i,j,n,size,flag,s0,s1;

    int* sigma = (int*)malloc(sizeof(int)*m->nsite);
    int Ma=0;
    int Mb=0;
    size = 0;
    for(i=0;i<m->nsite;++i){
        size += c[i]->n;
        sigma[i] = c[i]->state;
        if(i<xy) Ma+=c[i]->state;
        else Mb+=c[i]->state;
    }

    if(msize==0){
        msize = size;
        mtau  = (double*)malloc(sizeof(double)*msize);
        msite = (size_t*)malloc(sizeof(size_t)*msize);
        msort = (size_t*)malloc(sizeof(size_t)*msize);
    }
    else if(size>msize){
        msize = size;
        mtau  = (double*)realloc(mtau ,sizeof(double)*msize);
        msite = (size_t*)realloc(msite,sizeof(size_t)*msize);
        msort = (size_t*)realloc(msort,sizeof(size_t)*msize);
    }

    int t5a = 0;
    int t5b = 0;
    int t6a = 0;
    int t6b = 0;
    uint64_t key;
    int spin_id,type;
    item* it;
    for(i=0;i<m->nsite;++i){
        size = c[i]->size;
        flag = c[i]->flag;
        for(j=0;j<c[i]->n;++j){
            spin_id = c[i]->node[size*flag+j].spin_id;
            key = c[i]->node[size*flag+j].key;
            if(spin_id==0 && key!=UINT64_MAX){
                it = table_search_from_key(t,key);
                type = it->type;

                if(type==6 && i<xy) ++t6a;
                if(type==6 && i>=xy) ++t6b;
            }

            if(spin_id==3 && key!=UINT64_MAX){
                if(i<xy) ++t5a;
                if(i>=xy) ++t5b;
            }
        }
    }

    n=0;
    for(i=0;i<m->nsite;++i){
        for(j=0;c[i]->n>j;++j){
            size = c[i]->size;
            flag = c[i]->flag;

            s0 = c[i]->node[flag*size+j].state[0];
            s1 = c[i]->node[flag*size+j].state[1];
            if(s0!=s1){
                mtau[n]  = c[i]->node[flag*size+j].tau;
                msite[n] = (size_t)i;
                ++n;
            }
        }
    }

    gsl_sort_index(msort,mtau,1,n);

    double Ma2=0;
    double Mb2=0;
    for(j=0;j<n;++j){
        i = msort[j];
        if(msite[i]<xy)
            Ma += -sigma[msite[i]]*2;
        else
            Mb += -sigma[msite[i]]*2;

        sigma[msite[i]] *= -1;

        Ma2 += Ma*Ma;
        Mb2 += Mb*Mb;
    }

    Ma2 = Ma2/n;
    Mb2 = Mb2/n;

    char fname[128];
    sprintf(fname,"data/qlm_x_%d_y_%d_beta_%.1f_lambda_%.2f_seed_%d_.txt",x,y,m->beta,lambda,seed);
    FILE* myfile = fopen(fname,"a");
    fprintf(myfile,"%d %d %.10e %.10e %d\n",Ma,Mb,Ma2,Mb2,n);
    fclose(myfile);

    free(sigma);
}


int main(int argc, char** argv){
    int x;
    int y;
    double lambda;
    double beta;
    int ntherm;
    int nsweep;
    int seed;

    if(argc<7){
        x = 8;
        y = 8;
        lambda = 0.5;
        beta = 20;
        ntherm = 0;
        nsweep = 100000;
        seed = 1;
    }
    else{
        x = atoi(argv[1]);
        y = atoi(argv[2]);
        lambda = atof(argv[3]);
        beta = atof(argv[4]);
        ntherm = atoi(argv[5]);
        nsweep = atoi(argv[6]);
        seed = atoi(argv[7]);
    }

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
    double dis = gsl_rng_uniform_pos(rng);
    for(int i=0;i<x*y;++i){
       c[i] = chain_alloc(2000);
       if(dis<0.5) c[i]->state = 1;
       else c[i]->state = -1;
    }
    dis = gsl_rng_uniform_pos(rng);
    for(int i=x*y;i<2*x*y;++i){
       c[i] = chain_alloc(2000);
       if(dis<0.5) c[i]->state = 1;
       else c[i]->state = -1;
    }


    table* t = table_alloc(10);

    for(int i=0;i<ntherm;++i){
        update_A(c,t,m,x,y,rng);
        update_B(c,t,m,x,y,rng);

        int ng = 0;
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

                s0 = c[j]->node[flag*size+k].state[0];
                s1 = c[j]->node[flag*size+k].state[1];
                if(s0!=s1){
                    ++ng;
                }
            }

            //chain_print_state(c[j]);
        }
        //table_print_state(t);

        //printf("average cut in temporal direction (per site) : %.3f\n",((double)ng)/m->nsite);
    }

    for(int i=0;i<nsweep;++i){
        update_A(c,t,m,x,y,rng);
        update_B(c,t,m,x,y,rng);
        //update_A(c,t,m,x,y,rng);
        //update_B(c,t,m,x,y,rng);

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

        //printf("%d %d\n",Ma,Mb);

        qlm_measurement(c,t,m,x,y,lambda,seed);
    }

    return 0;
}

