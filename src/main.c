#include "configuration.h"
#include "graph.h"
#include "update_graphs.h"
#include "loop.h"


int main(int argc, char** argv){
    int i;
    size_t x=64;
    size_t init_size = 100;
    double beta = 10;
    int nsweep=1000;
    int seed=3979274;

    kinks* ks[x];
    bond*  bd[x];

//  setting lattice
    int nsite = x;
    int nbond = x;
    for(i=0;i<x;++i){
        ks[i] = kinks_alloc(init_size);
        kinks_set_sigma_i(ks[i],1);
    
        int site_id[2] = {i,(i+1)%x};
        graph* g[2] = {&GRAPH_HORI,&GRAPH_DIAG};
        double w[2] = {1.0,1.0};
        bd[i] = bond_alloc(init_size,4,2);
        bond_set_site_id(bd[i],4,site_id);
        bond_set_graph(bd[i],2,g,w);
    }

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);
    for(int is=0;is<nsweep;++is){
        update_graph_user_friendly(ks,bd,nsite,nbond,4,beta,rng);

        //printf("\n\nsweep = %d\n",is);
        //for(i=0;i<x;++i){
        //    kinks_print_state(ks[i]);
        //}
        //for(i=0;i<x;++i){
        //    bond_print_state(bd[i]);
        //}

        loop_cluster_update_user_friendly(ks,bd,nsite,nbond,rng);
        
        //printf("\n\n");
        //for(i=0;i<x;++i){
        //    kinks_print_state(ks[i]);
        //}

    }


//  free memory
    for(i=0;i<x;++i){
        kinks_free(ks[i]);
        bond_free(bd[i]);
    }
    
    return 0;
}
