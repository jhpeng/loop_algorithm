#include "model.h"

void model_afm_heisenderg_evolution(chain** c,table* t, model* m, gsl_rng* rng){
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

    loops_traverse(t,rng);
    loops_update_table(t);

    for(site_id=0;site_id<nsite;++site_id)
        loops_update_chain(c[site_id],t,rng);
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
        weight[bond_id] = 1.0;

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
