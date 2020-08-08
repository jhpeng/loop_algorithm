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
