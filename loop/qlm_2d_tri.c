#include "qlm_2d_tri.h"

static void update_A(chain** c, table* t, model* m, gsl_rng* rng){
    int i,j,k,l,bond_id;
    int nsite = m->nsite;
    double beta = m->beta;
    double weight;

    int nsq = nsite/2;

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

    for(i=0;i<nsq;++i)
        cluster_link_vertex(c[i],t);

    cluster_traverse(t,rng);

    cluster_update_table(t);

    for(i=0;i<nsite;++i)
        claster_update_chain(c[i],t);
}

static void update_B(chain** c, table* t, model* m, gsl_rng* rng){
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
       c[i] = chain_alloc(100);
       c[i]->state = 1;
    }

    table* t = table_alloc(20);

    for(int i=0;i<ntherm;++i){
        update_A(c,t,m,rng);
        update_B(c,t,m,rng);

        for(int j=0;j<m->nsite;++j){
            int n = c[j]->n;
            if(n>0){
                int size = c[j]->size;
                int flag = c[j]->flag;
                int s0 = c[j]->node[flag*size].state[0];
                int s1 = c[j]->node[flag*size+n-1].state[1];

                if(s0!=s1){
                    printf("Vailoate the periodic boundary condition!\n");
                    exit(1);
                }
            }
        }
    }

    for(int i=0;i<nsweep;++i){
        update_A(c,t,m,rng);
        update_B(c,t,m,rng);

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

