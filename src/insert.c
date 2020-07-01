#include "insert.h"

static double* insert_taus;
static int insert_ntau;
static int insert_taus_length;

static void generate_uniform_sampling(double weight, double beta, gsl_rng* rng){
    int i=0;
    double dis,tau;
    dis = gsl_rng_uniform_pos(rng);
    tau = -log(dis)/weight;

    while(tau<beta){
        if(i>=insert_taus_length){
            printf("generate_uniform_pos : run out of the storage!\n");
            exit(1);
        }
        insert_taus[i] = tau;
        tau = tau - log(dis)/weight;
        ++i;
    }

    insert_ntau = i;
}

void insert_horizontal_graphs(bond** b, site** s, int bond_id, int graph_id, double beta, gsl_rng* rng){
    assert((b[bond_id])->nspin==2);

    if(insert_taus==NULL){
        insert_taus_length = 10000;
        insert_taus = (double*)malloc(sizeof(double)*insert_taus_length);
        if(insert_taus==NULL){
            printf("insert_horizontal_graphs : memory allocate fail!\n");
            exit(1);
        }
    }

    int site_id[2];
    int sigma[2];
    site_id[0] = bond_get_site_id(b[bond_id],0);
    site_id[1] = bond_get_site_id(b[bond_id],1);
    site_node* sn[2];
    sn[0] = (s[site_id[0]])->first;
    sn[1] = (s[site_id[1]])->first;
    sigma[0] = (s[site_id[0]])->sigma;
    sigma[1] = (s[site_id[1]])->sigma;

    double weight = (b[bond_id])->weight[graph_id];
    generate_uniform_sampling(weight,beta,rng);

    double tau,tau1,tau2;
    for(int i=0;i<insert_ntau;++i){
        tau = insert_taus[i];
        if(sn[0]!=NULL)
            tau1 = (sn[0])->tau;
    }
}
