#include "qlm_2d_tri.h"

double* mtau;
size_t* msite;
size_t* msort;
int msize=0;

#if 1
#define gauss_law
#endif

/**
 * @brief Calculate the charge for a given highvar configuration.
 * 
 * This function calculates the charge for a given highvar configuration h by comparing
 * the values of h at different indices and updating the charge accordingly.
 *
 * @param h Pointer to an integer array representing the highvar configuration.
 * @return Calculated charge for the given highvar configuration.
 */
static int highvar2charge(int* h){
    int charge=0;

    if(h[0]==h[1]) charge+=1;
    else charge-=1;
    if(h[2]==h[3]) charge+=1;
    else charge-=1;
    if(h[4]==h[5]) charge+=1;
    else charge-=1;
    if(h[1]==h[2]) charge-=1;
    else charge+=1;
    if(h[3]==h[4]) charge-=1;
    else charge+=1;
    if(h[5]==h[0]) charge-=1;
    else charge+=1;

    return charge;
}

/**
 * @brief Calculate the X flux for a given highvar configuration.
 *
 * This function calculates the X flux for a given highvar configuration h by comparing
 * the values of h at different indices and updating the xflux accordingly.
 *
 * @param h Pointer to an integer array representing the highvar configuration.
 * @return Calculated X flux for the given highvar configuration.
 */
static double highvar2Xflux(int* h){
    double xflux=0.0;
    double factor=0.70710678118;

    if(h[0]==h[1]) xflux-=factor;
    else xflux+=factor;
    if(h[1]==h[2]) xflux+=1.0;
    else xflux-=1.0;
    if(h[3]==h[4]) xflux-=factor;
    else xflux+=factor;
    if(h[4]==h[5]) xflux+=1.0;
    else xflux-=1.0;

    return xflux;
}

/**
 * @brief Calculate the Y flux for a given highvar configuration.
 *
 * This function calculates the Y flux for a given highvar configuration h by comparing
 * the values of h at different indices and updating the yflux accordingly.
 *
 * @param h Pointer to an integer array representing the highvar configuration.
 * @return Calculated Y flux for the given highvar configuration.
 */
static double highvar2Yflux(int* h){
    double yflux=0.0;
    double factor=0.70710678118;

    if(h[0]==h[1]) yflux-=factor;
    else yflux+=factor;
    if(h[2]==h[3]) yflux+=1.0;
    else yflux-=1.0;
    if(h[3]==h[4]) yflux-=factor;
    else yflux+=factor;
    if(h[5]==h[0]) yflux+=1.0;
    else yflux-=1.0;

    return yflux;
}

/**
 * @brief Count the charge in a given state and print the result.
 * 
 * This function counts the charge in a given state for each block by calculating
 * the highvar configuration and calling highvar2charge function. The result is printed.
 *
 * @param state Pointer to an integer array representing the state.
 * @param x Width of the lattice.
 * @param y Height of the lattice.
 * @param distance Distance between charges.
 */
void counting_charge(int* state, int x, int y, int distance){
    int ai,aj,ak,bi,bj,bk;
    int block_id;
    int ix,iy,i;

    int h[6];
    int charge;


    printf("---------------- counting charge ------------------\n");
    for(block_id=0;block_id<x*y;++block_id){
        ix = block_id%x;
        iy = block_id/x;

        ai = iy*x+ix;
        aj = ((iy+1)%y)*x+ix;
        ak = ((iy+1)%y)*x+(ix+1)%x;

        bi = iy*x+ix               + x*y;
        bj = iy*x+(ix+1)%x         + x*y;
        bk = ((iy+1)%y)*x+(ix+1)%x + x*y;

        h[0] = state[ai];
        h[1] = state[bi];
        h[2] = state[aj];
        h[3] = state[bk];
        h[4] = state[ak];
        h[5] = state[bj];

        charge = highvar2charge(h);

        //if(block_id%x==x-1) printf("%d\n",charge);
        //else if(charge>=0) printf("%d  ",charge);
        //else printf("%d ",charge);


        for(i=0;i<6;++i) h[i] = (h[i]+1)/2;
        if(charge!=0)
            printf("block_id=%d state=(%d %d %d %d %d %d) charge=%d\n",block_id,h[0],h[1],h[2],h[3],h[4],h[5],charge);
    }
}

/**
 * @brief Print the charge distribution for a given state.
 *
 * This function prints the charge distribution for a given state by calculating
 * the highvar configuration and calling highvar2charge, highvar2Xflux, and highvar2Yflux functions.
 *
 * @param state Pointer to an integer array representing the state.
 * @param x Width of the lattice.
 * @param y Height of the lattice.
 */
void print_charge(int* state, int x, int y){
    int ai,aj,ak,bi,bj,bk;
    int block_id;
    int ix,iy;

    int h[6];
    int charge;
    double xflux,yflux;


    printf("---------------- printing charge ------------------\n");
    for(block_id=0;block_id<x*y;++block_id){
        ix = block_id%x;
        iy = block_id/x;

        ai = iy*x+ix;
        aj = ((iy+1)%y)*x+ix;
        ak = ((iy+1)%y)*x+(ix+1)%x;

        bi = iy*x+ix               + x*y;
        bj = iy*x+(ix+1)%x         + x*y;
        bk = ((iy+1)%y)*x+(ix+1)%x + x*y;

        h[0] = state[ai];
        h[1] = state[bi];
        h[2] = state[aj];
        h[3] = state[bk];
        h[4] = state[ak];
        h[5] = state[bj];

        charge = highvar2charge(h);
        xflux  = highvar2Xflux(h);
        yflux  = highvar2Yflux(h);

        if(charge>0) printf("+ ");
        else if(charge<0) printf("- ");
        else if(0){
            if(yflux==0){
                if(xflux>0) printf("> ");
                else if(xflux<0) printf("< ");
                else printf("  ");
            }
            else if(xflux==0) {
                if(yflux>0) printf("v ");
                else if(yflux<0)printf("^ ");
                else printf("  ");
            }
            else if(xflux*yflux>0) printf("%c ",(char)32);
            else if(xflux*yflux<0) printf("%c ",(char)32);
        }
        else printf("0 ");

        if(block_id%x==x-1) printf("\n");
    }
}

/**
 * @brief Initialize a state with charge 1.
 *
 * This function initializes a state with charge 1 by setting the state array values based on
 * the given lattice dimensions and distance between charges.
 *
 * @param state Pointer to an integer array representing the state.
 * @param x Width of the lattice.
 * @param y Height of the lattice.
 * @param distance Distance between charges.
 */
void initial_state_with_charge1(int* state, int x, int y, int distance){
    int ak,bj;
    int ix,iy;

    if(distance > x/2){
        printf("initial_state_with_charge : the distance can not larger than x/2!\n");
        exit(1);
    }

    for(iy=0;iy<y;iy++){
        for(ix=0;ix<x;ix++){
            state[iy*x+ix]     = 1;
            state[iy*x+ix+x*y] = 1;
        }
    }

    iy = y/2-1;
    for(ix=(x-distance)/2;ix<(x+distance)/2;ix++){
        ak = ((iy+1)%y)*x+(ix+1)%x;

        bj = iy*x+(ix+1)%x + x*y;

        state[ak] = -1;
        state[bj] = -1;
    }

    counting_charge(state,x,y,distance);
}

/**
 * @brief Initialize a state with charge 2.
 *
 * This function initializes a state with charge 2 by setting the state array values based on
 * the given lattice dimensions and distance between charges.
 *
 * @param state Pointer to an integer array representing the state.
 * @param x Width of the lattice.
 * @param y Height of the lattice.
 * @param distance Distance between charges.
 */
void initial_state_with_charge2(int* state, int x, int y, int distance){
    int ak,bk;
    int ix,iy;

    if(distance > x/2){
        printf("initial_state_with_charge : the distance can not larger than x/2!\n");
        exit(1);
    }

    for(iy=0;iy<y;iy++){
        for(ix=0;ix<x;ix++){
            state[iy*x+ix]     = 1;
            state[iy*x+ix+x*y] = -1;
        }
    }

    iy = y/2-1;
    for(ix=(x-distance)/2;ix<(x+distance)/2+1;ix++){
        ak = ((iy+1)%y)*x+(ix+1)%x;

        bk = ((iy+1)%y)*x+(ix+1)%x + x*y;

        state[ak] *= -1;
        state[bk] *= -1;
    }

    state[ak] *= -1;

    counting_charge(state,x,y,distance);
}

/**
 * @brief Apply Gauss's law on sub-lattice A.
 *
 * This function applys the Gauss's law on sub-lattice A for a given model and lattice
 * configuration. It modifies the chain and table structures by inserting new items
 * and updating their properties.
 *
 * @param c  A 2D array of pointers to chain structures representing the lattice configuration.
 * @param t  A pointer to a table structure storing lattice item information.
 * @param m  A pointer to a model structure containing the lattice model parameters.
 * @param x  The number of columns in the lattice.
 * @param y  The number of rows in the lattice.
 */
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
    int h[6];
    int charge,type;
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

        h[0] = c[ai]->state;
        h[1] = c[bi]->state;
        h[2] = c[aj]->state;
        h[3] = c[bk]->state;
        h[4] = c[ak]->state;
        h[5] = c[bj]->state;

        charge=highvar2charge(h);
        type=4;
        if(charge!=0) type=3;

        if(!((sb[0]==sb[1]) && (sb[1]==sb[2]))){
            if(sb[0]==sb[1]){
                key = table_generate_key(t);
                chain_insert(c[aj],&beta,&key,1,0);
                chain_insert(c[ak],&beta,&key,1,1);

                it = table_search_from_key(t,key);
                it->key   = key;
                it->type  = type;
                it->nspin = 2;
                it->state[0] = c[aj]->state;
                it->state[2] = c[aj]->state;
                it->state[1] = c[ak]->state;
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
                it->type  = type;
                it->nspin = 2;
                it->state[0] = c[ai]->state;
                it->state[2] = c[ai]->state;
                it->state[1] = c[aj]->state;
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
                it->type  = type;
                it->nspin = 2;
                it->state[0] = c[ai]->state;
                it->state[2] = c[ai]->state;
                it->state[1] = c[ak]->state;
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

/**
 * @brief Apply Gauss's law on sub-lattice B.
 *
 * This function applys the Gauss's law on sub-lattice B for a given model and lattice
 * configuration. It modifies the chain and table structures by inserting new items
 * and updating their properties.
 *
 * @param c  A 2D array of pointers to chain structures representing the lattice configuration.
 * @param t  A pointer to a table structure storing lattice item information.
 * @param m  A pointer to a model structure containing the lattice model parameters.
 * @param x  The number of columns in the lattice.
 * @param y  The number of rows in the lattice.
 */
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
    int h[6];
    int charge,type;
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

        h[0] = c[ai]->state;
        h[1] = c[bi]->state;
        h[2] = c[aj]->state;
        h[3] = c[bk]->state;
        h[4] = c[ak]->state;
        h[5] = c[bj]->state;

        charge=highvar2charge(h);
        type=4;
        if(charge!=0) type=3;

        if(!((sa[0]==sa[1]) && (sa[1]==sa[2]))){
            if(sa[0]==sa[1]){
                key = table_generate_key(t);
                chain_insert(c[bj],&beta,&key,1,0);
                chain_insert(c[bk],&beta,&key,1,1);

                it = table_search_from_key(t,key);
                it->key   = key;
                it->type  = type;
                it->nspin = 2;
                it->state[0] = c[bj]->state;
                it->state[2] = c[bj]->state;
                it->state[1] = c[bk]->state;
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
                it->type  = type;
                it->nspin = 2;
                it->state[0] = c[bi]->state;
                it->state[2] = c[bi]->state;
                it->state[1] = c[bj]->state;
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
                it->type  = type;
                it->nspin = 2;
                it->state[0] = c[bi]->state;
                it->state[2] = c[bi]->state;
                it->state[1] = c[bk]->state;
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

/**
 * @brief Updates the sub-lattice A of the quantum loop model (QLM).
 *
 * The function updates the sub-lattice A of the QLM on a 2D triangular lattice based on
 * given model parameters and lattice configuration. It updates chain structures,
 * inserts triangular cut graphs, and performs cluster updates.
 *
 * @param c A 2D array of pointers to chain structures representing the lattice configuration.
 * @param t A pointer to a table structure storing lattice item information.
 * @param m A pointer to a model structure containing the lattice model parameters.
 * @param x The number of columns in the lattice.
 * @param y The number of rows in the lattice.
 * @param rng A pointer to a GSL random number generator.
 */
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

#ifdef gauss_law
    gauss_law_A(c,t,m,x,y);
#endif

    for(i=0;i<nsq;++i)
        cluster_link_vertex(c[i],t);

    cluster_traverse(t,rng);

    cluster_update_table(t);

    for(i=0;i<2*nsq;++i)
        cluster_update_chain(c[i],t);

    for(i=0;i<nsq;++i){
        if(c[i]->n==0){
            if(gsl_rng_uniform_pos(rng)<0.5)
                c[i]->state *= -1;
        }
    }
}

/**
 * @brief Updates the sub-lattice B of the quantum loop model (QLM).
 *
 * The function updates the sub-lattice B of the QLM on a 2D triangular lattice based on
 * given model parameters and lattice configuration. It updates chain structures,
 * inserts triangular cut graphs, and performs cluster updates.
 *
 * @param c A 2D array of pointers to chain structures representing the lattice configuration.
 * @param t A pointer to a table structure storing lattice item information.
 * @param m A pointer to a model structure containing the lattice model parameters.
 * @param x The number of columns in the lattice.
 * @param y The number of rows in the lattice.
 * @param rng A pointer to a GSL random number generator.
 */
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

#ifdef gauss_law
    gauss_law_B(c,t,m,x,y);
#endif

    for(i=nsq;i<2*nsq;++i)
        cluster_link_vertex(c[i],t);

    cluster_traverse(t,rng);

    cluster_update_table(t);

    //for(i=nsq;i<2*nsq;++i)
    for(i=0;i<2*nsq;++i)
        cluster_update_chain(c[i],t);

    for(i=nsq;i<2*nsq;++i){
        if(c[i]->n==0){
            if(gsl_rng_uniform_pos(rng)<0.5)
                c[i]->state *= -1;
        }
    }
}

/**
 * @brief Generates a 2D triangular quantum loop model (QLM) with given parameters.
 *
 * This function creates a 2D triangular QLM with the specified dimensions, beta, and lambda.
 * It initializes the model structure, including bond types, bond-to-site mapping, and weights.
 *
 * @param x The number of columns in the lattice.
 * @param y The number of rows in the lattice.
 * @param beta The inverse temperature parameter for the model.
 * @param lambda The loop weight parameter for the model.
 * @return A pointer to the generated model structure.
 */
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

/**
 * @brief Compute the local energy density for a given site.
 *
 * This function calculates the local energy density for a given site
 * based on the input data provided. It is used as a helper function
 * in other calculations.
 *
 * @param c Pointer to an array of pointers to chains, size should be 4.
 * @param lambda The coupling constant used in the calculations.
 * @param beta The inverse temperature used in the calculations.
 * @return The computed local energy density for the given site.
 */
double local_energy_density(chain** c, double lambda, double beta){
    int state[3];
    int i,j,n,site,ref;
    int flag,size,s0,s1;

    double tau,tau1,stau;

    n=0;
    for(i=0;i<3;++i){
        size = c[i]->size;
        flag = c[i]->flag;

        state[i] = c[i]->state;

        for(j=0;j<c[i]->n;++j){
            s0 = c[i]->node[flag*size+j].state[0];
            s1 = c[i]->node[flag*size+j].state[1];
            if(s0!=s1){
                mtau[n]  = c[i]->node[flag*size+j].tau;
                msite[n] = i;
                ++n;
            }
        }
    }

    gsl_sort_index(msort,mtau,1,n);

    ref=0;
    if((state[0]==state[1]) && (state[1]==state[2])){
        tau1=0;
        ref=1;
    }

    stau=0;
    for(i=0;i<n;++i){
        j = msort[i];

        tau = mtau[j];
        site = msite[j];

        state[site] *= -1;

        if(ref){
            stau += (tau-tau1);
            ref=0;
        }
        else{
            if((state[0]==state[1]) && (state[1]==state[2])){
                tau1=tau;
                ref=1;
            }
        }
    }

    if(ref){
        stau += (beta-tau1);
    }

    size = c[3]->size;
    flag = c[3]->flag;
    n=0;
    for(j=0;j<c[3]->n;++j){
        s0 = c[3]->node[flag*size+j].state[0];
        s1 = c[3]->node[flag*size+j].state[1];
        if(s0!=s1){
            ++n;
        }
    }

    return lambda*stau/beta+(double)n/beta;
}

/**
 * @brief Perform qlm measurements on the system.
 *
 * This function computes various observables for a given system,
 * stores the results in a block of data, and writes the results to
 * a file once the block is full.
 *
 * @param c Pointer to an array of pointers to chains.
 * @param t Pointer to a table.
 * @param m Pointer to a model.
 * @param x The x dimension of the lattice.
 * @param y The y dimension of the lattice.
 * @param lambda The coupling constant used in the calculations.
 * @param fname The filename to which the results should be written.
 */
int qlm_block_data_id=0;
int qlm_block_data_size=1000;
double* qlm_block_data;
void qlm_measurement(chain** c, table* t, model* m, int x, int y, double lambda, char* fname){
    int xy = x*y;
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

    if(qlm_block_data==NULL){
        qlm_block_data = (double*)malloc(sizeof(double)*qlm_block_data_size*20);
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


    n=0;
    for(i=0;i<m->nsite;++i){
        size = c[i]->size;
        flag = c[i]->flag;

        for(j=0;j<c[i]->n;++j){
            s0 = c[i]->node[flag*size+j].state[0];
            s1 = c[i]->node[flag*size+j].state[1];
            if(s0!=s1){
                mtau[n]  = c[i]->node[flag*size+j].tau;
                msite[n] = i;
                ++n;
            }
        }
    }

    gsl_sort_index(msort,mtau,1,n);

    double Ma2=0;
    double Mb2=0;

#if 1
    int na=0;
    int nb=0;
    double stau_a=0.0;
    double stau_b=0.0;
    for(j=0;j<n;++j){
        i = msort[j];

        if(msite[i]<xy){
            Ma2 += (double)Ma*(double)Ma*(mtau[i]-stau_a);
            stau_a = mtau[i];
            Ma += -sigma[msite[i]]*2;
            na++;
        }
        else{
            Mb2 += (double)Mb*(double)Mb*(mtau[i]-stau_b);
            stau_b = mtau[i];
            Mb += -sigma[msite[i]]*2;
            nb++;
        }

        sigma[msite[i]] *= -1;


        //printf("%d %d %ld %.3f %d %d \n",j,i,msite[i],mtau[i],Ma,Mb);

    }
    Ma2 += (double)Ma*(double)Ma*(m->beta-stau_a);
    Mb2 += (double)Mb*(double)Mb*(m->beta-stau_a);

    Ma2 = Ma2/(m->beta);
    Mb2 = Mb2/(m->beta);
#else 
    Ma2 = (double)Ma*(double)Ma;
    Mb2 = (double)Mb*(double)Mb;
#endif

    chain *c_temp[4];
    double energy=0;
    for(i=0;i<m->nsite;++i){
        c_temp[0] = c[m->bond2site[i*NSPIN_MAX+0]];
        c_temp[1] = c[m->bond2site[i*NSPIN_MAX+1]];
        c_temp[2] = c[m->bond2site[i*NSPIN_MAX+2]];
        c_temp[3] = c[m->bond2site[i*NSPIN_MAX+3]];

        energy+=local_energy_density(c_temp,lambda, m->beta);
    }
    energy = energy/m->nsite;

    qlm_block_data[qlm_block_data_id+qlm_block_data_size*0] = Ma2;
    qlm_block_data[qlm_block_data_id+qlm_block_data_size*1] = Mb2;
    qlm_block_data[qlm_block_data_id+qlm_block_data_size*2] = fabs(Ma);
    qlm_block_data[qlm_block_data_id+qlm_block_data_size*3] = fabs(Mb);
    qlm_block_data[qlm_block_data_id+qlm_block_data_size*4] = (double)Ma*(double)Ma;
    qlm_block_data[qlm_block_data_id+qlm_block_data_size*5] = (double)Mb*(double)Mb;
    qlm_block_data[qlm_block_data_id+qlm_block_data_size*6] = (double)Ma*(double)Ma*(double)Ma*(double)Ma;
    qlm_block_data[qlm_block_data_id+qlm_block_data_size*7] = (double)Mb*(double)Mb*(double)Mb*(double)Mb;
    qlm_block_data[qlm_block_data_id+qlm_block_data_size*8] = (double)n;
    qlm_block_data[qlm_block_data_id+qlm_block_data_size*9] = energy;

    qlm_block_data_id++;

    if(qlm_block_data_id==qlm_block_data_size){
        int nobs=10;
        double obs[nobs];
        FILE* myfile = fopen(fname,"a");
        fprintf(myfile,"%d %d ",Ma,Mb);
        for(i=0;i<nobs;i++){
            obs[i]=0;
            for(j=0;j<qlm_block_data_size;j++){
                obs[i] += qlm_block_data[j+qlm_block_data_size*i];
            }
            obs[i] = obs[i]/qlm_block_data_size;

            fprintf(myfile,"%.16e ",obs[i]);
        }

        fprintf(myfile,"\n");
        fclose(myfile);

        qlm_block_data_id=0;
    }

    free(sigma);
}

/**
 * @brief Calculate the energy map for the lattice.
 *
 * This function calculates the energy map of the lattice by calling
 * local_energy_density for each site and accumulating the results.
 * After a specified number of sweeps, the energy map is written to
 * a file and memory is released.
 *
 * @param c Pointer to an array of pointers to chains.
 * @param m Pointer to a model.
 * @param x The x dimension of the lattice.
 * @param y The y dimension of the lattice.
 * @param lambda The coupling constant used in the calculations.
 * @param nsweep The number of sweeps to be performed before writing to file.
 * @param fname The filename to which the results should be written.
 */
int energy_map_id;
double* energy_map_data;
void qlm_energy_map(chain** c, model* m, int x, int y, double lambda, int nsweep, char* fname){
    int i;
    chain* c_temp[4];

    if(energy_map_data==NULL){
        energy_map_id=0;
        energy_map_data = (double*)malloc(sizeof(double)*m->nsite);
        for(i=0;i<m->nsite;++i) energy_map_data[i]=0;
    }

    energy_map_id++;

    for(i=0;i<m->nsite;++i){
        c_temp[0] = c[m->bond2site[i*NSPIN_MAX+0]];
        c_temp[1] = c[m->bond2site[i*NSPIN_MAX+1]];
        c_temp[2] = c[m->bond2site[i*NSPIN_MAX+2]];
        c_temp[3] = c[m->bond2site[i*NSPIN_MAX+3]];

        energy_map_data[i]+=local_energy_density(c_temp,lambda, m->beta);
    }

    if(energy_map_id==nsweep){
        FILE* myfile = fopen(fname,"w");
        for(i=0;i<m->nsite;++i){
            energy_map_data[i] /= (double)nsweep;
            fprintf(myfile,"%.10e ",energy_map_data[i]);
        }
        fprintf(myfile,"\n");
        fclose(myfile);

        free(energy_map_data);
    }
}

/**
 * @brief Check if the reference configurations are valid.
 *
 * This function checks if the reference configurations in a given table
 * satisfy the required conditions. If not, it prints an error message.
 * It is used for debugging purposes.
 *
 * @param t Pointer to a table.
 */
void qlm_check_ref_conf(table* t){
    int size = t->size;
    int n = t->n;
    int count=0;
    int nspin,type;
    int* state;

    for(int i=0;i<size;++i){
        if(t->list[i].key!=UINT64_MAX){
            type = nspin = t->list[i].type;
            if(type==5 || type==6){
                nspin = t->list[i].nspin;
                state = t->list[i].state;

                if(!((state[0]==state[1]) && (state[1]==state[2]))){
                    printf("violate the reference configuration!\n");
                }
                else if(!((state[0+nspin]==state[1+nspin]) && (state[1+nspin]==state[2+nspin]))){
                    printf("violate the reference configuration!\n");
                }
            }
            count++;
        }
    }

    if((n-count)!=0)
        printf("missing %d\n",n-count);
}

/**
 * @brief Main function for simulating a quantum lattice model (QLM) on a 2D triangular lattice.
 *
 * @param argc Number of input arguments.
 * @param argv Array of input arguments as strings.
 *        - argv[1]: Lattice size x.
 *        - argv[2]: Lattice size y.
 *        - argv[3]: Coupling constant -lambda.
 *        - argv[4]: Inverse temperature beta.
 *        - argv[5]: Distance for the initial state.
 *        - argv[6]: Number of thermalization steps (ntherm).
 *        - argv[7]: Number of sweeps for measurement (nsweep).
 *        - argv[8]: QLM block data size (optional, not used in the current code).
 *        - argv[9]: Seed for the random number generator.
 *        - argv[10]: Output data file name (optional, not used in the current code).
 *        - argv[11]: Output energy map file name (optional, not used in the current code).
 * @return 0 on successful completion of the simulation.
 *
 * The function simulates a QLM on a 2D triangular lattice, taking input arguments for lattice size,
 * coupling constants (lambda and beta), distance, number of thermalization steps (ntherm),
 * number of sweeps (nsweep), and random number generator seed. It generates a lattice model,
 * initializes a random state, and updates the state through a series of thermalization and
 * measurement steps.
 */
int main(int argc, char** argv){
    int x;
    int y;
    double lambda;
    double beta;
    int distance;
    int ntherm;
    int nsweep;
    int seed;
    char fname[128];
    char hname[128];
    char ename[128];
    char cname[128];

    if(argc==1){
        x = 8;
        y = 8;
        lambda = 1.0;
        beta = 10;
        distance = 0;
        ntherm = 10000;
        nsweep = 4000;
        seed = 1;
    }
    else{
        x = atoi(argv[1]);
        y = atoi(argv[2]);
        lambda = atof(argv[3]);
        beta = atof(argv[4]);
        distance = atoi(argv[5]);
        ntherm = atoi(argv[6]);
        nsweep = atoi(argv[7]);
        qlm_block_data_size = atoi(argv[8]);
        seed = atoi(argv[9]);
    }

    sprintf(fname,"data/qlm_x_%d_y_%d_beta_%.2f_lambda_%.4f_cdist_%d_seed_%d_.txt",x,y,beta,lambda,distance,seed);
    sprintf(hname,"data/qlm_energy_histogram_x_%d_y_%d_beta_%.2f_lambda_%.4f_cdist_%d_seed_%d_.txt",x,y,beta,lambda,distance,seed);
    //sprintf(fname,"data.txt");
    sprintf(ename,"map/qlm_energy_map_x_%d_y_%d_beta_%.2f_lambda_%.4f_cdist_%d_seed_%d_.txt",x,y,beta,lambda,distance,seed);
    sprintf(cname,"map/config_map_x_%d_y_%d_beta_%.2f_lambda_%.4f_cdist_%d_seed_%d_.txt",x,y,beta,lambda,distance,seed);

     //if(argc>10){
     //    strcpy(fname,argv[10]);
     //}
     //else if(argc>11){
     //    strcpy(ename,argv[11]);
     //}

    model* m = generate_QLM_2d_triangular(x,y,beta,lambda);

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    //state for generate initial state or charge counting
    int* state = (int*)malloc(sizeof(int)*x*y*2);

    chain* c[m->nsite];

    if(distance>0){
        initial_state_with_charge1(state,x,y,distance);

        for(int i=0;i<2*x*y;++i){
            c[i] = chain_alloc(2000);
            c[i]->state = state[i];
        }
    }
    else {
        double dis;
        dis = gsl_rng_uniform_pos(rng);
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
    }
    

    table* t = table_alloc(16);

    int quick_thermalize=0;
    int nq = 5;
    if( quick_thermalize){
        model mq[nq];
        for(int i=0;i<nq;++i){
            mq[i] = *m;
            for(int j=0;j<i;++j) (mq[i]).beta *= 0.5;
        }

        for(int i=0;i<nq;++i){
            for(int j=0;j<ntherm/nq;++j){
                update_A(c,t,&(mq[nq-i-1]),x,y,rng);
                update_B(c,t,&(mq[nq-i-1]),x,y,rng);
            }
        }
    }

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
                    printf("Vaiolate the periodic boundary condition!\n");
                    printf("# of node = %d\n",n);
                    printf("site=%d node_id=%d (s0,s1)=(%d,%d)\n",j,k,s0,s1);
                    chain_print_state(c[j]);
                    table_print_state(t);
                    printf("---------------------------------------\n");
                    exit(1);
                }

                s0 = c[j]->node[flag*size+k].state[0];
                s1 = c[j]->node[flag*size+k].state[1];
                if(s0!=s1){
                    ++ng;
                }
            }

        }

        for(int site_id=0;site_id<m->nsite;++site_id){
            state[site_id] = c[site_id]->state;
        }
        //print_charge(state,x,y);
        //counting_charge(state,x,y,distance);
        //qlm_check_ref_conf(t);
    }
    free(state);

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

        qlm_measurement(c,t,m,x,y,lambda,fname);
        qlm_energy_map(c,m,x,y,lambda,nsweep,ename);

        if((i+1)%2000==0) {
            // record the sanp shoot of configuration
            FILE* config_file = fopen(cname,"a");

            for(int i_site=0;i_site<(m->nsite);i_site++) {
                fprintf(config_file,"%d ",c[i_site]->state);
            }
            fprintf(config_file,"\n");

            fclose(config_file);

            // record the the current energy density
            FILE* energy_hist_file = fopen(hname,"a");

            chain *c_temp[4];
            double energy=0;
            for(i=0;i<m->nsite;++i){
                c_temp[0] = c[m->bond2site[i*NSPIN_MAX+0]];
                c_temp[1] = c[m->bond2site[i*NSPIN_MAX+1]];
                c_temp[2] = c[m->bond2site[i*NSPIN_MAX+2]];
                c_temp[3] = c[m->bond2site[i*NSPIN_MAX+3]];

                energy+=local_energy_density(c_temp,lambda, m->beta);
            }
            energy = energy/m->nsite;

            fprintf(energy_hist_file,"%.16e \n", energy);
            fclose(energy_hist_file);

        }

    }

    return 0;
}

