#include "update_graphs.h"

static void uniform_dist_sequential_generation(double** x, size_t* size, size_t* nx, double xmin, double xmax, double lambda, gsl_rng* rng){
    double tau = xmin;
    double dis;
    *nx=0;
    dis = gsl_rng_uniform_pos(rng);
    tau = tau-log(dis)/lambda;
    while(tau<xmax){
        if((*nx+1)>(*size)){
            *x = (double*)realloc(*x,2*(*size)*sizeof(double));
            *size = (*size)*2;
        }
        (*x)[*nx] = tau;

        (*nx)++;

        dis = gsl_rng_uniform_pos(rng);
        tau = tau-log(dis)/lambda;
    }
}

static void insert_graph_and_kinks(kinks** ks, bond** bd, int bond_id, int graph_id, double tau, int (*rule_nspin)(int*,int*)){
    int nspin = bond_get_nspin(bd[bond_id]);
    int sigma[nspin/2];
    int s[nspin];
    int site_id;
    for(int i=0;i<nspin/2;++i){
        site_id = bond_get_site_id(bd[bond_id],i);
        sigma[i] = kinks_sigma_from_tau(ks[site_id],tau);
        if(sigma[i]==INT_MAX) 
            return;
        
        s[i] = sigma[i];
        s[i+nspin/2] = sigma[i];
    }

    graph* g = bond_get_graph(bd[bond_id],graph_id);
    if(rule_nspin(g->rule,s)){
        int type_id = bond_check_available_id(bd[bond_id]);
        int kink_id[nspin/2];
        for(int i=0;i<nspin/2;++i){
            site_id = bond_get_site_id(bd[bond_id],i);
            kink_id[i] = kinks_insert(ks[site_id],bond_id,(nspin/2)*type_id+i,sigma[i],tau);
        }
        type_id = bond_insert_graph(bd[bond_id],graph_id,tau, kink_id);
    }
}

static void remove_graph_and_kinks(kinks** ks, bond** bd, int bond_id, int type_id){
    int nspin = bond_get_nspin(bd[bond_id]);
    int check_no_kink=1;
    int site_id,kink_id;

    for(int i=0;i<nspin/2;++i){
        site_id = bond_get_site_id(bd[bond_id],i);
        kink_id = bond_get_kink_id(bd[bond_id],type_id,i);
        if(!(kinks_check_no_kink(ks[site_id],kink_id)))
            check_no_kink = 0;
    }

    if(check_no_kink){
        for(int i=0;i<nspin/2;++i){
            site_id = bond_get_site_id(bd[bond_id],i);
            kink_id = bond_get_kink_id(bd[bond_id],type_id,i);
            kinks_remove(ks[site_id],kink_id);
        }

        bond_remove_graph(bd[bond_id],type_id);
    }
}

static void group_insert_graphs_and_kinks(kinks** ks, bond** bd, int bond_id, int graph_id, const double *tau, size_t ntau, int (*rule_nspin)(int*,int*)){
    size_t ntype = bond_get_ntype(bd[bond_id]);
    size_t size  = bond_get_size(bd[bond_id]);
    int nspin = bond_get_nspin(bd[bond_id]);
    int ngraph = bond_get_ngraph(bd[bond_id]);
    if(ntype+ntau>size){
        bond* bd_temp = bond_alloc(ntype+ntau,nspin,ngraph);
        bond_memcpy(bd_temp,bd[bond_id]);
        bond_free(bd[bond_id]);
        bd[bond_id] = bd_temp;
    }
    size_t nkink;
    for(int i=0;i<nspin/2;++i){
        int site_id = bond_get_site_id(bd[bond_id],i);
        size  = kinks_get_size(ks[site_id]);
        nkink = kinks_get_nkink(ks[site_id]);
        if(nkink+ntau>size){
            kinks* ks_temp = kinks_alloc(nkink+ntau);
            kinks_memcpy(ks_temp,ks[site_id]);
            kinks_free(ks[site_id]);
            ks[site_id] = ks_temp;
        }
    }

    for(int i=0;i<ntau;++i){
        insert_graph_and_kinks(ks,bd,bond_id,graph_id,tau[i],rule_nspin);
    }
}

void generate_graphs_with_uniform_dist(kinks** ks, bond** bd, int bond_id, double beta, double** taus, size_t* size, gsl_rng* rng, int (*rule_nspin)(int*,int*)){
    int ngraph = bond_get_ngraph(bd[bond_id]);
    size_t nt;
    double lambda;

    for(int graph_id=0; graph_id<ngraph;++graph_id){
        lambda = bond_get_weight(bd[bond_id],graph_id);
        uniform_dist_sequential_generation(taus,size,&nt,0,beta,lambda,rng);

        group_insert_graphs_and_kinks(ks,bd,bond_id,graph_id,*taus,nt,rule_nspin);
        
    }

    int nspin = bond_get_nspin(bd[bond_id]);
    for(int i=0;i<nspin/2;++i){
        int site_id = bond_get_site_id(bd[bond_id],i);
        kinks_sort_index_with_tau(ks[site_id]);
    }
}

void remove_all_graphs_with_no_kink(kinks** ks, bond** bd, int nbond){
    for(int bond_id=0;bond_id<nbond;++bond_id){
        size_t size  = bond_get_size(bd[bond_id]);
    
        for(int type_id=0;type_id<size;++type_id){
            int type = bond_get_type(bd[bond_id],type_id);
            if(type!=-1)
                remove_graph_and_kinks(ks,bd,bond_id,type_id);
        }
    }
}

#if 1

int main(int argc, char** argv){
    size_t size = 20;
    int nspin = 4;
    int ngraph=2;

    int site_id[2] = {0,1};
    graph* g[2] = {&GRAPH_HORI,&GRAPH_DIAG};
    double w[2] = {1.0,1.0};

    bond* bd1 = bond_alloc(size,nspin,ngraph);

    bond_set_site_id(bd1,nspin,site_id);
    bond_set_graph(bd1,ngraph,g,w);

    kinks* ks[2];
    ks[0] = kinks_alloc(size);
    ks[1] = kinks_alloc(size);

    kinks_set_sigma_i(ks[0],1);
    kinks_set_sigma_i(ks[1],1);


/*-------------------------------------------------------------------------*/
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,236298);

    size_t ns = 10;
    double beta = 100.0;
    double* taus = (double*)malloc(size*sizeof(double));

    for(int j=0;j<1000;++j){
        remove_all_graphs_with_no_kink(ks,&bd1,1);
        generate_graphs_with_uniform_dist(ks,&bd1,0,beta,&taus,&ns,rng,rule_4spin);
    }

/*-------------------------------------------------------------------------*/

    int ntype = bond_get_ntype(bd1);
    int nkink0 = kinks_get_nkink(ks[0]);
    int nkink1 = kinks_get_nkink(ks[1]);

    printf("ntype=%d | nkink0=%d | nkink1=%d\n",ntype,nkink0,nkink1);

    size = bond_get_size(bd1);
    for(size_t type_id=0; type_id<size; ++type_id){
        int  type  = bond_get_type(bd1,type_id);
        double tau = bond_get_tau(bd1,type_id);
        int kink_id0 = bond_get_kink_id(bd1,type_id,0);
        int kink_id1 = bond_get_kink_id(bd1,type_id,1);

        printf("%zu %d %e %d %d\n",type_id,type,tau,kink_id0,kink_id1);
    }

    size = kinks_get_size(ks[0]);
    for(size_t kink_id=0;kink_id<size;++kink_id){
        int active0  = kinks_get_active(ks[0],kink_id);
        int active1  = kinks_get_active(ks[1],kink_id);
        int bond_id0 = kinks_get_bond_id(ks[0],kink_id);
        int bond_id1 = kinks_get_bond_id(ks[1],kink_id);
        int type_id0 = kinks_get_type_id(ks[0],kink_id);
        int type_id1 = kinks_get_type_id(ks[1],kink_id);
        double tau0  = kinks_get_tau(ks[0],kink_id);
        double tau1  = kinks_get_tau(ks[1],kink_id);

        printf("%zu | %d %d %d %e | %d %d %d %e\n",kink_id,active0,bond_id0,type_id0,tau0,active1,bond_id1,type_id1,tau1);
    }

    bond_free(bd1);
    kinks_free(ks[0]);
    kinks_free(ks[1]);
}

#endif
