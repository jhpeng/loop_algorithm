#include "model.h"

void model_afm_heisenderg_update(chain** c,table* t, model* m, gsl_rng* rng){
    int i,j,bond_id,site_id,type;
    int nsite = m->nsite;
    int nbond = m->nbond;
    double beta = m->beta;
    double weight;

    for(bond_id=0;bond_id<nbond;++bond_id){
        type = m->type[bond_id];
        if(type==1){
            i = m->bond2site[bond_id*NSPIN_MAX+0];
            j = m->bond2site[bond_id*NSPIN_MAX+1];
            weight = m->weight[bond_id];

            insert_horizontal_graph(c[i],c[j],t,weight,beta,rng);
        }
    }

    for(site_id=0;site_id<nsite;++site_id)
        loops_link_vertex(c[site_id],t);

    int nspin;
    for(i=0;i<t->size;++i){
        if(t->list[i].key!=UINT64_MAX){
            nspin = t->list[i].nspin;
            for(j=0;j<2*nspin;++j)
                assert(t->list[i].link_spin[j]>=0);
        }
    }

    loops_traverse(t,rng);

    for(i=0;i<t->size;++i){
        if(t->list[i].key!=UINT64_MAX){
            nspin = t->list[i].nspin;
            for(j=0;j<2*nspin;++j)
                assert(t->list[i].link_spin[j]<0);
        }
    }

    loops_update_table(t);

    for(site_id=0;site_id<nsite;++site_id)
        loops_update_chain(c[site_id],t,rng);

    //for(int i=0;i<m->nsite;++i)
    //    chain_print_state(c[i]);

    //table_print_state(t);
}

model* model_generate_afm_heisenberg_isotropy(int x, int y, double beta){
    int nsite = x*y;
    int nbond = 2*nsite;

    int* type = (int*)malloc(sizeof(int)*nbond);
    int* bond2site = (int*)malloc(sizeof(int)*nbond*NSPIN_MAX);
    double* weight = (double*)malloc(sizeof(double)*nbond);

    int i,j,ix,iy;
    for(int bond_id=0;bond_id<nbond;++bond_id){
        type[bond_id] = 1;
        weight[bond_id] = 0.5;

        if(bond_id<nsite){
            ix = bond_id%x;
            iy = bond_id/x;

            i = iy*x+ix;
            j = iy*x+(ix+1)%x;

            bond2site[bond_id*NSPIN_MAX+0] = i;
            bond2site[bond_id*NSPIN_MAX+1] = j;
        }
        else{
            ix = (bond_id%nsite)%x;
            iy = (bond_id%nsite)/x;

            i = iy*x+ix;
            j = ((iy+1)%y)*x+ix;

            bond2site[bond_id*NSPIN_MAX+0] = i;
            bond2site[bond_id*NSPIN_MAX+1] = j;
        }
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

void model_QLM_triangular_2d_update(chain** c, table* t, model* m, gsl_rng* rng){
    int i,j,k,l,bond_id,type;
    int nsite = m->nsite;
    int nbond = m->nbond;
    double beta = m->beta;
    double weight;

    chain* c_temp[4];
    for(bond_id=0;bond_id<nbond;++bond_id){
        type = m->type[bond_id];
        if(type==5){
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
        else if(type==6){
            i = m->bond2site[bond_id*NSPIN_MAX+0];
            j = m->bond2site[bond_id*NSPIN_MAX+1];
            k = m->bond2site[bond_id*NSPIN_MAX+2];
            weight = m->weight[bond_id];

            c_temp[0] = c[i];
            c_temp[1] = c[j];
            c_temp[2] = c[k];

            insert_triangular_graph(c_temp,t,weight,beta,rng);
        }
    }

    for(i=0;i<nsite;++i)
        cluster_link_vertex(c[i],t);

    cluster_traverse(t,rng);

    cluster_update_table(t);

    for(i=0;i<nsite;++i)
        claster_update_chain(c[i],t);

}

model* model_generate_QLM_triangular_2d(int x, int y, double beta){
    int nsite = 2*x*y;
    int nbond = 4*x*y;

    double W = 1.0;

    int* type = (int*)malloc(sizeof(int)*nbond);
    int* bond2site = (int*)malloc(sizeof(int)*nbond*NSPIN_MAX);
    double* weight = (double*)malloc(sizeof(double)*nbond);

    int i,j,k,ix,iy,bond_id;
    for(bond_id=0;bond_id<x*y;++bond_id){
        type[bond_id] = 5;
        weight[bond_id] = W;

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
        weight[bond_id+x*y] = W;

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
        weight[bond_id+2*x*y] = W;

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
        weight[bond_id+3*x*y] = W;

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

int main(){
    int x = 16;
    int y = 16;
    double beta = 40.0;
    int seed = 2913;

    model* m = model_generate_QLM_triangular_2d(x,y,beta);

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

    for(int i=0;i<2000;++i){
        model_QLM_triangular_2d_update(c,t,m,rng);

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

    for(int i=0;i<20000;++i){
        model_QLM_triangular_2d_update(c,t,m,rng);

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

int model_afm_heisenberg_test(int argc, char** argv){
    int x=atoi(argv[1]);
    int y=atoi(argv[2]);
    double beta = atof(argv[3]);
    int nthermal = 20000;
    int nblock = 1000;
    int nsweep = 1000;
    int nc = (int)beta;
    int scale = 16;
    int seed = 984293;

    model* m  = model_generate_afm_heisenberg_isotropy(x,y,beta);

    chain* c[m->nsite];
    for(int i=0;i<m->nsite;++i){
        c[i] = chain_alloc(nc);
        c[i]->state = 2*(i%2)-1;
    }

    table* t = table_alloc(scale);

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    double* sfactor = (double*)malloc(sizeof(double)*(m->nsite));
    for(int i=0;i<x;++i){
        for(int j=0;j<y;++j){
            sfactor[j*x+i] = ((i+j)%2)*2-1;
        }
    }


    for(int i=0;i<nthermal;++i){
        model_afm_heisenderg_update(c,t,m,rng);
    }

    double mz2,mz1,mz;
    double ms2,ms1,ms;
    double ng;
    for(int i_block=0;i_block<nblock;++i_block){
        mz2 = 0;
        mz1 = 0;
        ms2 = 0;
        ms1 = 0;
        ng = 0;
        for(int i_sweep=0;i_sweep<nsweep;++i_sweep){
            model_afm_heisenderg_update(c,t,m,rng);
            mz = 0;
            ms = 0;
            for(int i=0;i<m->nsite;++i){
                mz += c[i]->state;
                ms += (c[i]->state)*sfactor[i];
            }

            mz1 += mz;
            mz2 += mz*mz;
            ms1 += fabs(ms);
            ms2 += ms*ms;

            int nng=0;
            for(int i=0;i<m->nsite;++i)
                nng += (double)c[i]->n;

            ng += (double)nng;
        }

        mz2 = mz2/nsweep*0.25;
        mz1 = mz1/nsweep*0.5;
        ms2 = ms2/nsweep*0.25/(m->nsite)/(m->nsite);
        ms1 = ms1/nsweep*0.5/(m->nsite);
        ng = ng/nsweep*0.5;

        printf("%.8f %.8f %.8f %.8f %.3f\n",mz2*beta/(m->nsite),mz2,ms1,ms2,ng);
    }

    return 0;
}
