#include "configuration.h"

kinks* kinks_alloc(int size){
    kinks* ks = (kinks*)malloc(sizeof(kinks));
    ks->size     = size;
    ks->nkink    = 0;
    ks->sigma_i  = 0;     
    ks->active   = (int*)malloc(size*sizeof(int));
    ks->sigma_b  = (int*)malloc(size*sizeof(int));
    ks->sigma_a  = (int*)malloc(size*sizeof(int));
    ks->bond_id  = (int*)malloc(size*sizeof(int));
    ks->type_id = (int*)malloc(size*sizeof(int));
    ks->sort     = (size_t*)malloc(size*sizeof(size_t));
    ks->tau = (double*)malloc(size*sizeof(double));

    for(int i=0;i<size;++i){
        ks->active[i]   = 0;
        ks->sigma_b[i]  = 0;
        ks->sigma_a[i]  = 0;
        ks->bond_id[i]  = -1;
        ks->type_id[i] = -1;
        ks->sort[i]     = i;
        ks->tau[i] = DBL_MAX;;
    }

    return ks;
}

void kinks_free(kinks* ks){
    free(ks->active);
    free(ks->sigma_b);
    free(ks->sigma_a);
    free(ks->bond_id);
    free(ks->type_id);
    free(ks->sort);
    free(ks->tau);
    free(ks);
}

void kinks_memcpy(kinks* dest, const kinks* src){
    assert(dest->size>=src->size);

    dest->nkink   = src->nkink;
    dest->sigma_i = src->sigma_i;

    for(int i=0;i<src->size;++i){
        dest->active[i]   = src->active[i];
        dest->sigma_b[i]  = src->sigma_b[i];
        dest->sigma_a[i]  = src->sigma_a[i];
        dest->bond_id[i]  = src->bond_id[i];
        dest->type_id[i] = src->type_id[i];
        dest->tau[i] = src->tau[i];
        dest->sort[i] = src->sort[i];
    }
}

void kinks_set_sigma_i(kinks* ks, int sigma_i){
    ks->sigma_i = sigma_i;
}

void kinks_flip_sigma_b(kinks* ks, int kink_id){
    assert(ks->active[kink_id]);

    ks->sigma_b[kink_id] *= -1;
}

void kinks_flip_sigma_a(kinks* ks, int kink_id){
    assert(ks->active[kink_id]);

    ks->sigma_a[kink_id] *= -1;
}

int kinks_sigma_from_tau(const kinks* ks, double tau){
    int sigma=0;
    if(ks->nkink==0) sigma = ks->sigma_i;
    else{
        int start = 0;
        int end = ks->size-1;
        int mid;
        int check=1;

        if(tau<ks->tau[ks->sort[start]]){
            sigma = ks->sigma_b[ks->sort[start]];
            check = 0;
        }
        else if(tau>ks->tau[ks->sort[end]]){
            sigma = ks->sigma_a[ks->sort[end]];
            check = 0;
        }
        
        while(check){
            mid = start + (end-start)/2;
            if(ks->tau[ks->sort[mid]]<tau) start=mid;
            else if(ks->tau[ks->sort[mid]]>tau) end=mid;
            else{
                sigma = INT_MAX;
                check = 0;
            }

            if((end-start)==1){
                sigma = ks->sigma_a[ks->sort[start]];
                check = 0;
            }
        }
    }

    return sigma;
}

int kinks_check_available_id(kinks* ks){
    int i, kink_id=-1;
    for(i=0;i<ks->size;++i){
        if(ks->active[i]==0){
            kink_id = i;
            break;
        }
    }

    assert(kink_id!=-1);

    return kink_id;
}

int kinks_check_no_kink(kinks* ks, int kink_id){
    assert(ks->active[kink_id]);
    if(ks->sigma_a[kink_id]==ks->sigma_b[kink_id]) return 1;
    else return 0;
}

int kinks_get_size(kinks* ks){
    return ks->size;
}

int kinks_get_nkink(kinks* ks){
    return ks->nkink;
}

int kinks_get_sigma_i(kinks* ks){
    return ks->sigma_i;
}

int kinks_get_active(kinks* ks, int kink_id){
    assert(kink_id<ks->size);

    return ks->active[kink_id];
}

int kinks_get_sigma_b(kinks* ks, int kink_id){
    assert(kink_id<ks->size);

    return ks->sigma_b[kink_id];
}

int kinks_get_sigma_a(kinks* ks, int kink_id){
    assert(kink_id<ks->size);

    return ks->sigma_a[kink_id];
}

int kinks_get_bond_id(kinks* ks, int kink_id){
    assert(kink_id<ks->size);

    return ks->bond_id[kink_id];
}

int kinks_get_type_id(kinks* ks, int kink_id){
    assert(kink_id<ks->size);

    return ks->type_id[kink_id];
}

double kinks_get_tau(kinks* ks, int kink_id){
    assert(kink_id<ks->size);

    return ks->tau[kink_id];
}

int kinks_get_sort(kinks* ks, int rank){
    assert(rank<ks->nkink);

    return ks->sort[rank];
}


int kinks_insert(kinks* ks, int bond_id, int type_id, int sigma, double tau){
    int i, kink_id=-1;
    for(i=ks->nkink;i<ks->size;++i){
        if(ks->active[ks->sort[i]]==0){
            kink_id = ks->sort[i];
            break;
        }
    }

    assert(kink_id!=-1);

    ks->active[kink_id]   = 1;
    ks->sigma_b[kink_id]  = sigma;
    ks->sigma_a[kink_id]  = sigma;
    ks->bond_id[kink_id]  = bond_id;
    ks->type_id[kink_id] = type_id;
    ks->tau[kink_id] = tau;
    ks->nkink++;

    return kink_id;
}

void kinks_remove(kinks* ks, int kink_id){
    assert(ks->sigma_b[kink_id]==ks->sigma_a[kink_id]);

    ks->active[kink_id]   = 0;
    ks->sigma_b[kink_id]  = 0;
    ks->sigma_a[kink_id]  = 0;
    ks->bond_id[kink_id]  = -1;
    ks->type_id[kink_id] = -1;
    ks->tau[kink_id] = DBL_MAX;
    ks->nkink--;
}

void kinks_sort_index_with_tau(kinks* ks){
   gsl_sort_index(ks->sort,ks->tau,1,ks->size); 
}

void kinks_print_state(kinks* ks){
    kinks_sort_index_with_tau(ks);
    int size,rank_id,kink_id,bond_id,type_id,nkink;
    int sigma_a,sigma_b,sigma_i;
    double tau;
    
    size = kinks_get_size(ks);
    nkink = kinks_get_nkink(ks);
    sigma_i = kinks_get_sigma_i(ks);

    printf("##################################################\n");
    printf("State of this kinks...\n");
    printf("size = %d,  nkink = %d, sigma_i = %d \n",size,nkink,sigma_i);
    for(rank_id=0;rank_id<nkink;++rank_id){
        kink_id = kinks_get_sort(ks,rank_id);
        sigma_b = kinks_get_sigma_b(ks,kink_id);
        sigma_a = kinks_get_sigma_a(ks,kink_id);
        bond_id = kinks_get_bond_id(ks,kink_id);
        type_id = kinks_get_type_id(ks,kink_id);
        tau = kinks_get_tau(ks,kink_id);

        printf("%.3f (%d,%d) bond_id=%d type_id=%d kink_id=%d\n",tau,sigma_b,sigma_a,bond_id,type_id,kink_id);
    }
    printf("##################################################\n");
}
