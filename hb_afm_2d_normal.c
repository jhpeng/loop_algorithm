#include "configuration.h"
#include "graph.h"
#include "update_graphs.h"
#include "loop.h"
#include "lattice.h"


double mean(double* x, int nx){
    double sum=0;
    for(int i=0;i<nx;++i)
        sum += x[i];

    return sum/nx;
}

double calc_mz_square(lattice* lat){
    int nsite = lat->nsite;
    double mz = 0;
    for(int i=0;i<nsite;++i)
        mz += kinks_get_sigma_i(lat->ks[i]);

    return mz*mz*0.25;
}

int main(int argc, char** argv){
    int lx = 8;
    int ly = 8;
    double beta = 30.0;
    int seed = 983894;
    int nthermal = 2000;
    int block_size = 2000;
    int nblock = 20;

    // setting observable
    double* mz2 = (double*)malloc(sizeof(double)*block_size);
    
    lattice* lat = lattice_afm_heisenberg_2d_uniform(lx,ly,seed);
    lattice_set_beta(lat,beta);

    lattice_monte_carlo_sweep(lat,nthermal);
    for(int i_block=0;i_block<nblock;++i_block){
        for(int i=0;i<block_size;++i){
            lattice_monte_carlo_sweep(lat,1);

            mz2[i]  = calc_mz_square(lat);
        }

        double usus_m = mean(mz2,block_size)*(lat->beta)/(lat->nsite);
        double mz2d_m = mean(mz2,block_size)/(lat->nsite)/(lat->nsite);

        printf("%lf %lf\n",mz2d_m,usus_m);
    }

    lattice_free(lat);
    free(mz2);
}
