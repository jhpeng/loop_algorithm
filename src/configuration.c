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
    ks->graph_id = (int*)malloc(size*sizeof(int));
    ks->sort     = (size_t*)malloc(size*sizeof(size_t));
    ks->tau = (double*)malloc(size*sizeof(double));

    for(int i=0;i<size;++i){
        ks->active[i]   = 0;
        ks->sigma_b[i]  = 0;
        ks->sigma_a[i]  = 0;
        ks->bond_id[i]  = -1;
        ks->graph_id[i] = -1;
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
    free(ks->graph_id);
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
        dest->graph_id[i] = src->graph_id[i];
        dest->tau[i] = src->tau[i];
        dest->sort[i] = src->sort[i];
    }
}

void kinks_set_sigma_i(kinks* ks, int sigma_i){
    ks->sigma_i = sigma_i;
}

int kinks_sigma_from_tau(const kinks* ks, double tau){
    int sigma=0;
    if(ks->nkink==0) sigma = ks->sigma_i;
    else{
        int start = 0;
        int end = ks->size-1;
        int mid;
        int check=1;
        
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

int kinks_insert(kinks* ks, int bond_id, int graph_id, int sigma, double tau){
    int i, kink_id=-1;
    for(i=0;i<ks->size;++i){
        if(ks->active[i]==0){
            kink_id = i;
            break;
        }
    }

    assert(kink_id!=-1);

    ks->active[kink_id]   = 1;
    ks->sigma_b[kink_id]  = sigma;
    ks->sigma_a[kink_id]  = sigma;
    ks->bond_id[kink_id]  = bond_id;
    ks->graph_id[kink_id] = graph_id;
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
    ks->graph_id[kink_id] = -1;
    ks->tau[kink_id] = DBL_MAX;
    ks->nkink--;
}

void kinks_sort_index_with_tau(kinks* ks){
   gsl_sort_index(ks->sort,ks->tau,1,ks->size); 
}

#if 1
#include <gsl/gsl_rng.h>

int main(int argc, char** argv){
    int size=100;
    kinks* ks = kinks_alloc(size);
    kinks* ks_cpy = kinks_alloc(size);
    kinks_set_sigma_i(ks,1);
    
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,0);

    for(int i=0;i<50;++i){
        double tau = 10.0*gsl_rng_uniform_pos(rng);
        int kink_id = kinks_insert(ks,0,i,1,tau);
        printf("%d\n",kink_id);
    }
    kinks_remove(ks,10);
    kinks_remove(ks,20);

    double tau = 10.0*gsl_rng_uniform_pos(rng);
    int kink_id = kinks_insert(ks,0,50,1,tau);
    printf("%d\n",kink_id);

    kinks_sort_index_with_tau(ks);

    for(int i=0;i<100;++i){
        int active = ks->active[i];
        int sigma_a = ks->sigma_a[i];
        int sigma_b = ks->sigma_b[i];
        int bond_id = ks->bond_id[i];
        int graph_id = ks->graph_id[i];
        double tau = ks->tau[i];
        size_t sort = ks->sort[i];
        printf("%d, %d, %d, %d, %d, %.8e, %zu\n",active,sigma_a,sigma_b,bond_id,graph_id,tau,sort);
    }

    kinks_memcpy(ks_cpy,ks);

    kinks_free(ks_cpy);
    kinks_free(ks);
}
#endif
