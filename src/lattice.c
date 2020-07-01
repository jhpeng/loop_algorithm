#include "lattice.h"

lattice* lattice_afm_heisenberg_2d_uniform(int lx, int ly, int seed){
    if((lx%2)!=0 || (ly%2)!=0){
        printf("lattice_afm_heisenberg_2d_uniform : lx and ly should be even number!\n");
    }

    int dim = 2;
    int nspin = 4;
    int nsite = lx*ly;
    int nbond = 2*lx*ly;
    int init_size = 1000;

    lattice* lat = (lattice*)malloc(sizeof(lattice));
    lat->shape = (int*)malloc(sizeof(int)*dim);
    lat->ks = (kinks**)malloc(sizeof(kinks*)*nsite);
    lat->bd = (bond** )malloc(sizeof(bond* )*nbond);

    lat->rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(lat->rng,seed);

    lat->dim  = dim;
    lat->max_nspin = nspin;
    lat->shape[0] = lx;
    lat->shape[1] = ly;
    lat->nsite = nsite;
    lat->nbond = nbond;

    int x,y,i,j;
    for(x=0;x<lx;++x){
        for(y=0;y<ly;++y){
            i = y*lx+x;
            lat->ks[i] = kinks_alloc(init_size);
            kinks_set_sigma_i(lat->ks[i],((x+y)%2)*2-1);
        }
    }

    int site_id[2];
    graph* g[1] = {&GRAPH_HORI};
    double w[1] = {1.0};
    for(x=0;x<lx;++x){
        for(y=0;y<ly;++y){
            i = y*lx+x;

            lat->bd[i] = bond_alloc(init_size,nspin,1);
            j = y*lx+(x+1)%lx;
            site_id[0] = i;
            site_id[1] = j;
            bond_set_site_id(lat->bd[i],4,site_id);
            bond_set_graph(lat->bd[i],1,g,w);


            lat->bd[i+nsite] = bond_alloc(init_size,nspin,1);
            j = ((y+1)%ly)*lx+x;
            site_id[0] = i;
            site_id[1] = j;
            bond_set_site_id(lat->bd[i+nsite],4,site_id);
            bond_set_graph(lat->bd[i+nsite],1,g,w);

        }
    }

    return lat;
}

void lattice_set_beta(lattice* lat, double beta){
    lat->beta = beta;
}

void lattice_monte_carlo_sweep(lattice* lat, int nsweep){
    int i;
    int nsite = lat->nsite;
    int nbond = lat->nbond;
    int nspin = lat->max_nspin;
    double beta = lat->beta;

    for(i=0;i<nsweep;++i){
        update_graph_user_friendly(lat->ks,lat->bd,nsite,nbond,nspin,beta,lat->rng);
        loop_cluster_update_user_friendly(lat->ks,lat->bd,nsite,nbond,lat->rng);
    }
}

void lattice_free(lattice* lat){
    free(lat->shape);
    for(int i=0;i<(lat->nsite);++i){
        kinks_free(lat->ks[i]);
    }

    for(int i=0;i<(lat->nbond);++i){
        bond_free(lat->bd[i]);
    }

    gsl_rng_free(lat->rng);
    free(lat->ks);
    free(lat->bd);
    free(lat);
}
