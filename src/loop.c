#include "loop.h"

// global variables
static link* outer_link;
static link* inner_link;

static int get_link_id(int nbond, int max_nspin, int bond_id, int type_id, int spin_id){
    return max_nspin*(type_id*nbond+bond_id)+spin_id;
}

link* link_alloc(int nbond, int max_type_id, int max_nspin){
    link* lk = (link*)malloc(sizeof(link));
    lk->nbond = nbond;
    lk->max_type_id = max_type_id;
    lk->max_nspin = max_nspin;
    lk->size = nbond*max_type_id*max_nspin;

    lk->index = (int*)malloc((lk->size)*sizeof(int));

    for(int i=0;i<lk->size;++i) lk->index[i]=-1;

    return lk;
}

void link_realloc(link* lk, int max_type_id){
    lk->max_type_id = max_type_id;
    lk->size = (lk->nbond)*max_type_id*(lk->max_nspin);

    lk->index = (int*)realloc(lk->index,(lk->size)*sizeof(int));
    for(int i=0;i<lk->size;++i) lk->index[i]=-1;
}

void link_free(link* lk){
    free(lk->index);
    free(lk);
}

int link_check_init(link* lk){
    int check=1;
    for(int i=0;i<lk->size;++i){
        if(lk->index[i]!=-1){
            check=0;
            break;
        }
    }

    return check;
}

void link_measure_size(bond** bd, int nbond, int* max_type_id, int* max_nspin){
    *max_nspin=0;
    *max_type_id=0;
    int nspin;
    int intype;
    for(int bond_id=0;bond_id<nbond;++bond_id){
        nspin = bond_get_nspin(bd[bond_id]);
        intype = bond_get_size(bd[bond_id]);

        if(nspin>*max_nspin) *max_nspin = nspin;

        if(intype>*max_type_id) *max_type_id = intype;
    }
}

static void link_set_linking(link* lk, int l1, int l2){
    assert(l1<lk->size);
    assert(l2<lk->size);

    lk->index[l1] = l2;
    lk->index[l2] = l1;
}

static void link_set_one_way_linking(link* lk, int l1, int l2){
    assert(l1<lk->size);
    assert(l2<lk->size);

    lk->index[l1] = l2;
}

static int link_get_index(link* lk, int p){
    assert(p<lk->size);

    return lk->index[p];
}

static void link_erase_index(link* lk, int p){
    assert(p<lk->size);

    lk->index[p] = -1;
}

static void link_insert_flip_flag(link* lk, int p){
    assert(p<lk->size);

    lk->index[p] = -2;
}

void loop_construct_outer_link(link* lk, kinks** ks, bond** bd, int nsite, int nbond){
    int rank_id,kink_id,site_id,bond_id,type_id,spin_id,nkink;
    int max_nspin = lk->max_nspin;
    int nspin,fl,ll1,ll2;
    
    assert(nbond==lk->nbond);
    assert(link_check_init(lk));

    for(site_id=0;site_id<nsite;++site_id){
        nkink = kinks_get_nkink(ks[site_id]);
        if(nkink>0){
            kink_id = kinks_get_sort(ks[site_id],0);
        
            assert(kinks_get_active(ks[site_id],kink_id));

            bond_id = kinks_get_bond_id(ks[site_id],kink_id);
            type_id = kinks_get_type_id(ks[site_id],kink_id);
            nspin = bond_get_nspin(bd[bond_id]);
            spin_id = type_id%(nspin/2);
            type_id = type_id/(nspin/2);
            ll2 = get_link_id(nbond,max_nspin,bond_id,type_id,spin_id);
            fl = ll2;
            for(rank_id=1;rank_id<nkink;++rank_id){
                ll1 = ll2+nspin/2;
                kink_id = kinks_get_sort(ks[site_id],rank_id);

                assert(kinks_get_active(ks[site_id],kink_id));

                bond_id = kinks_get_bond_id(ks[site_id],kink_id);
                type_id = kinks_get_type_id(ks[site_id],kink_id);
                nspin = bond_get_nspin(bd[bond_id]);
                spin_id = type_id%(nspin/2);
                type_id = type_id/(nspin/2);

                ll2 = get_link_id(nbond,max_nspin,bond_id,type_id,spin_id);
                link_set_linking(lk,ll1,ll2);
            }
            ll1 = ll2+nspin/2;
            link_set_linking(lk,fl,ll1);
        }
    }
}


void loop_construct_inner_link(link* lk, bond** bd, int nbond){
    int max_nspin = lk->max_nspin;
    int l1,l2,size,bond_id,type_id,spin_id;
    int nspin,type;
    graph* g;
    
    assert(nbond==lk->nbond);
    assert(link_check_init(lk));

    for(bond_id=0;bond_id<nbond;++bond_id){
        nspin = bond_get_nspin(bd[bond_id]);
        size = bond_get_size(bd[bond_id]);
        for(type_id=0;type_id<size;++type_id){
            type = bond_get_type(bd[bond_id],type_id);
            if(type!=-1){
                g = bond_get_graph(bd[bond_id],type);
                for(spin_id=0;spin_id<nspin;++spin_id){
                    l1 = get_link_id(nbond,max_nspin,bond_id,type_id,spin_id);
                    l2 = l1 + g->link[spin_id] - spin_id;
                    link_set_one_way_linking(lk,l1,l2);
                }
            }
        }
    }
}

void loop_cluster_identify(link* outer_lk, link* inner_lk, gsl_rng* rng){
    assert((outer_lk->size)==(inner_lk->size));

    int scan_id,now_id,next_id;
    int size = outer_lk->size;
    double dis;
    
    for(scan_id=0;scan_id<size;++scan_id){
        if(link_get_index(outer_lk,scan_id)>=0){
            dis = gsl_rng_uniform_pos(rng);
            if(dis<0.5){
                now_id = scan_id;
                while(link_get_index(outer_lk,now_id)>=0){
                    next_id = link_get_index(outer_lk,now_id);
                    link_erase_index(outer_lk,now_id);
                    link_erase_index(outer_lk,next_id);

                    now_id = link_get_index(inner_lk,next_id);
                }
            }
            else{
                now_id = scan_id;
                while(link_get_index(outer_lk,now_id)>=0){
                    next_id = link_get_index(outer_lk,now_id);
                    link_insert_flip_flag(outer_lk,now_id);
                    link_insert_flip_flag(outer_lk,next_id);

                    now_id = link_get_index(inner_lk,next_id);
                }
            }
        }
    }

    for(scan_id=0;scan_id<size;++scan_id) 
        link_erase_index(inner_lk,scan_id);
}

void loop_cluster_flip(kinks** ks, bond** bd, link* outer_lk, int nbond){
    int max_nspin = outer_lk->max_nspin;
    int size,site_id,kink_id,bond_id,type_id,spin_id,link_id;
    int nspin,flag;
    
    assert(nbond==outer_lk->nbond);

    for(bond_id=0;bond_id<nbond;++bond_id){
        nspin = bond_get_nspin(bd[bond_id]);
        size = bond_get_size(bd[bond_id]);
        for(type_id=0;type_id<size;++type_id){
            if(bond_get_type(bd[bond_id],type_id)!=-1){
                for(spin_id=0;spin_id<nspin/2;++spin_id){
                    link_id = get_link_id(nbond,max_nspin,bond_id,type_id,spin_id);
                    flag = link_get_index(outer_lk,link_id);
                    if(flag==-2){
                        site_id = bond_get_site_id(bd[bond_id],spin_id);
                        kink_id = bond_get_kink_id(bd[bond_id],type_id,spin_id);
                        kinks_flip_sigma_b(ks[site_id],kink_id);
                        link_erase_index(outer_lk,link_id);
                    }
                }
                for(spin_id=0;spin_id<nspin/2;++spin_id){
                    link_id = get_link_id(nbond,max_nspin,bond_id,type_id,spin_id+nspin/2);
                    flag = link_get_index(outer_lk,link_id);
                    if(flag==-2){
                        site_id = bond_get_site_id(bd[bond_id],spin_id);
                        kink_id = bond_get_kink_id(bd[bond_id],type_id,spin_id);
                        kinks_flip_sigma_a(ks[site_id],kink_id);
                        link_erase_index(outer_lk,link_id);
                    }
                }
            }
        }
    }
}

void loop_cluster_sigma_i(kinks** ks, int nsite, gsl_rng* rng){
    int site_id,kink_id;
    double dis;
    int sigma;

    for(site_id=0;site_id<nsite;++site_id){
        if(kinks_get_nkink(ks[site_id])==0){
            dis = gsl_rng_uniform_pos(rng);
            if(dis<0.5) kinks_set_sigma_i(ks[site_id],1);
            else kinks_set_sigma_i(ks[site_id],-1);
        }
        else{
            kink_id = kinks_get_sort(ks[site_id],0);
            sigma = kinks_get_sigma_b(ks[site_id],kink_id);
            kinks_set_sigma_i(ks[site_id],sigma);
        }
    }
}

int loop_cluster_check_periodic(kinks** ks, int nsite){
    int site_id,kink_id,nkink;
    int check=1,sigma_i,sigma_f;

    for(site_id=0;site_id<nsite;++site_id){
        nkink = kinks_get_nkink(ks[site_id]);
        if(nkink!=0){
            kink_id = kinks_get_sort(ks[site_id],0);
            sigma_i = kinks_get_sigma_b(ks[site_id],kink_id);
            kink_id = kinks_get_sort(ks[site_id],nkink-1);
            sigma_f = kinks_get_sigma_a(ks[site_id],kink_id);

            if(sigma_i!=sigma_f){
                check=0;
                break;
            }
        }
    }

    return check;
}

void loop_cluster_update_user_friendly(kinks** ks, bond** bd, int nsite, int nbond, gsl_rng* rng){
    int max_type_id;
    int max_nspin;

    link_measure_size(bd,nbond,&max_type_id,&max_nspin);

    if(outer_link==NULL){ 
        outer_link = link_alloc(nbond,max_type_id,max_nspin);
        inner_link = link_alloc(nbond,max_type_id,max_nspin);
    }
    else if(outer_link->max_type_id<max_type_id){
        link_realloc(outer_link,max_type_id);
        link_realloc(inner_link,max_type_id);
    }

    assert(outer_link!=NULL);
    assert(inner_link!=NULL);
    assert((outer_link->size)==(inner_link->size));
    assert((outer_link->nbond)==(inner_link->nbond));
    assert((outer_link->max_type_id)==(inner_link->max_type_id));
    assert((outer_link->max_nspin)==(inner_link->max_nspin));

    //start the loop cluster update
    assert(loop_cluster_check_periodic(ks,nsite));

    loop_construct_outer_link(outer_link,ks,bd,nsite,nbond);
    loop_construct_inner_link(inner_link,bd,nbond);
    loop_cluster_identify(outer_link,inner_link,rng);
    loop_cluster_flip(ks,bd,outer_link,nbond);
    loop_cluster_sigma_i(ks,nsite,rng);

    assert(loop_cluster_check_periodic(ks,nsite));
}
