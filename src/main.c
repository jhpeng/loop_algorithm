#include "configuration.h"
#include "graph.h"
#include "update_graphs.h"
#include "loop.h"
#include "lattice.h"


int main(int argc, char** argv){
    int lx = 16;
    int ly = 16;
    double beta = 20.0;
    int seed = 983894;
    int nsweep = 2000;
    
    lattice* lat = lattice_afm_heisenberg_2d_uniform(lx,ly,seed);
    lattice_set_beta(lat,beta);

    lattice_monte_carlo_sweep(lat,nsweep);
    lattice_free(lat);
}