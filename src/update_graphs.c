#include "update_graphs.h"

insertion_plan* insertion_plan_alloc(int size, int max_nspin){
    insertion_plan* plan = (insertion_plan*)malloc(sizeof(insertion_plan));
    plan->size = size;
    plan->max_nspin = max_nspin;
    plan->ntau = 0;
    
    plan->taus = (double*)malloc(size*sizeof(double));
    plan->accept = (int*)malloc(size*sizeof(int));
    plan->sigma  = (int*)malloc(max_nspin/2*size*sizeof(int));

    for(int i=0;i<size;++i) plan->taus[i] = 0;

    return plan;
}

void insertion_plan_realloc(insertion_plan* plan, int size){
    assert(size>plan->size);
    plan->size = size;

    plan->taus = (double*)realloc(plan->taus,size*sizeof(double));
    plan->accept = (int*)realloc(plan->accept,size*sizeof(int));
    plan->sigma  = (int*)realloc(plan->sigma,plan->max_nspin/2*size*sizeof(int));

    for(int i=0;i<size;++i) plan->accept[i]=0;
}

void insertion_plan_free(insertion_plan* plan){
    free(plan->taus);
    free(plan->accept);
    free(plan->sigma);
    free(plan);
}

static void uniform_dist_sequential_generation(insertion_plan* plan, double xmin, double xmax, double lambda, gsl_rng* rng){
    double dis,x;
    plan->ntau = 0;

    dis = gsl_rng_uniform_pos(rng);
    x = xmin - log(dis)/lambda;
    while(x<xmax){
        if((plan->size)<=(plan->ntau)) 
            insertion_plan_realloc(plan,2*(plan->ntau+1));

        plan->taus[plan->ntau] = x;
        plan->ntau++;

        dis = gsl_rng_uniform_pos(rng);
        x = x - log(dis)/lambda;
    }

    assert((plan->ntau) < (plan->size));
}

void construct_insertion_plan(insertion_plan* plan, kinks** ks, bond** bd, double beta, int bond_id, int graph_id, gsl_rng* rng, int (*rule_nspin)(int*,int*)){

    double lambda = bond_get_weight(bd[bond_id],graph_id);
    uniform_dist_sequential_generation(plan,0,beta,lambda,rng);

    int nspin,spin_id,check;
    int site_id;
    int i;
    double tau;
    graph* g = bond_get_graph(bd[bond_id],graph_id);

    nspin = bond_get_nspin(bd[bond_id]);
    int sigma[nspin/2];
    int s[nspin];

    for(i=0;i<plan->ntau;++i){
        check=1;
        tau = plan->taus[i];
        for(spin_id=0;spin_id<nspin/2;++spin_id){
            site_id = bond_get_site_id(bd[bond_id],spin_id);
            sigma[spin_id] = kinks_sigma_from_tau(ks[site_id],tau);
            if(sigma[spin_id]==INT_MAX) 
                check=0;
        
            s[spin_id] = sigma[spin_id];
            s[spin_id+nspin/2] = sigma[spin_id];
        }

        if(rule_nspin(g->rule,s) && check){
            plan->accept[i] = 1;
            for(spin_id=0;spin_id<nspin/2;++spin_id){
                plan->sigma[(plan->max_nspin/2)*i+spin_id] = sigma[spin_id];
            }
        }
        else plan->accept[i] = 0;
    }
}

void insert_graph_and_kinks(kinks** ks, bond** bd, insertion_plan* plan, int bond_id, int graph_id){
    int i,tau_id,type_id,site_id;
    int ntau = plan->ntau;
    int nspin  = bond_get_nspin(bd[bond_id]);
    int ngraph = bond_get_ngraph(bd[bond_id]);
    int kink_id[nspin/2];
    double tau;

    int size_res;
    for(i=0;i<nspin/2;++i){
        site_id = bond_get_site_id(bd[bond_id],i);
        size_res = (ks[site_id])->nkink+ntau+100;
        if((ks[site_id])->size < size_res){
            kinks* ks_temp = kinks_alloc(size_res);
            kinks_memcpy(ks_temp,ks[site_id]);
            kinks_free(ks[site_id]);
            ks[site_id] = ks_temp;
        }
    }
    
    size_res = (bd[bond_id])->ntype+ntau+100;
    if((bd[bond_id])->size < size_res){
        bond* bd_temp = bond_alloc(size_res,nspin,ngraph);
        bond_memcpy(bd_temp,bd[bond_id]);
        bond_free(bd[bond_id]);
        bd[bond_id] = bd_temp;
    }

    int bsize = bond_get_size(bd[bond_id]);
    int* label = (int*)malloc(sizeof(int)*bsize);

    i=0;
    type_id=0;
    while(i<ntau){
        if(bond_get_type(bd[bond_id],type_id)==-1){
            label[i] = type_id;
            ++i;
        }
        ++type_id;
    }

    for(tau_id=0;tau_id<ntau;++tau_id){
        if(plan->accept[tau_id]){
            tau = plan->taus[tau_id];
            //type_id = bond_check_available_id(bd[bond_id]);
            type_id = label[tau_id];
            for(i=0;i<nspin/2;++i){
                site_id = bond_get_site_id(bd[bond_id],i);
                kink_id[i] = kinks_insert(ks[site_id],bond_id,(nspin/2)*type_id+i,plan->sigma[plan->max_nspin/2*tau_id+i],tau);
            }

            bond_insert_graph(bd[bond_id],graph_id,tau,kink_id,type_id);
        }
    }

    free(label);
    //for(i=0;i<nspin/2;++i){
    //    site_id = bond_get_site_id(bd[bond_id],i);
    //    kinks_sort_index_with_tau(ks[site_id]);
    //}
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

void remove_all_graphs_with_no_kink(kinks** ks, bond** bd, int nsite, int nbond){
    for(int bond_id=0;bond_id<nbond;++bond_id){
        int size  = bond_get_size(bd[bond_id]);
    
        for(int type_id=0;type_id<size;++type_id){
            int type = bond_get_type(bd[bond_id],type_id);
            if(type!=-1)
                remove_graph_and_kinks(ks,bd,bond_id,type_id);
        }
    }
    for(int site_id=0;site_id<nsite;++site_id)
        kinks_sort_index_with_tau(ks[site_id]);
}

static insertion_plan* user_plan;
void update_graph_user_friendly(kinks** ks, bond** bd, int nsite, int nbond, int max_nspin, double beta, gsl_rng* rng){
    if(user_plan==NULL){
        int size = 100;
        user_plan = insertion_plan_alloc(size,max_nspin);
    }

    remove_all_graphs_with_no_kink(ks,bd,nsite,nbond);

    assert((user_plan->max_nspin)==max_nspin);

    int bond_id,graph_id;
    int ngraph,nspin;
    int (*rule_nspin)(int*,int*);
    for(bond_id=0;bond_id<nbond;++bond_id){
        ngraph = bond_get_ngraph(bd[bond_id]);
        nspin  = bond_get_nspin(bd[bond_id]);
        
        //add new rule if nspin!=4
        if(nspin==4) rule_nspin = rule_4spin;
        else{
            printf("not yet to implement rule_nspin!\n");
            exit(-1);
        }
        for(graph_id=0;graph_id<ngraph;++graph_id){
            construct_insertion_plan(user_plan,ks,bd,beta,bond_id,graph_id,rng,rule_nspin);
            insert_graph_and_kinks(ks,bd,user_plan,bond_id,graph_id);

        }
    }
    for(int site_id=0;site_id<nsite;++site_id)
        kinks_sort_index_with_tau(ks[site_id]);
}
