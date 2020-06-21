#include "graph.h"

int rule_4spin(int* rule, int s0, int s1, int s2, int s3){
    s0 = (s0+1)/2;
    s1 = (s1+1)/2;
    s2 = (s2+1)/2;
    s3 = (s3+1)/2;
    return rule[s0+2*s1+4*s2+8*s3];
}

bond* bond_alloc(size_t size, int nspin, int ntype){
    bond* bd = (bond*)malloc(sizeof(bond));
    bd->size = size;
    bd->nspin = nspin;
    bd->ntype = ntype;
    bd->ngraph = 0;
    bd->site_id = (int*)malloc(nspin/2*sizeof(int));
    bd->type = (graph**)malloc(ntype*sizeof(graph*));
    bd->weight = (double*)malloc(ntype*sizeof(double));
    bd->graphs = (int*)malloc(size*sizeof(int));
    bd->tau = (double*)malloc(size*sizeof(double));
    bd->kink_id = (int*)malloc(nspin/2*size*sizeof(int));

    for(size_t i=0;i<size;++i){
        bd->graphs[i] = -1;
        bd->tau[i] = DBL_MAX;
        for(int j=0;j<nspin/2;++j) bd->kink_id[i*nspin/2+j]=-1;
    }

    return bd;
}

void bond_free(bond* bd){
    free(bd->site_id);
    free(bd->type);
    free(bd->weight);
    free(bd->graphs);
    free(bd->tau);
    free(bd->kink_id);
    free(bd);
}

void bond_memcpy(bond* dest, const bond* src){
    if(dest->size>=src->size){
        if(dest->nspin!=src->nspin){
            printf("bond_memcpy : memory copy fail! nspin should be the same!\n");
            exit(-1);
        }
        else if(dest->ntype!=src->ntype){
            printf("bond_memcpy : memory copy fail! ntype should be the same!\n");
            exit(-1);
        }

        int nspin = src->nspin;
        int ntype = src->ntype;
        size_t size = src->size;
        
        for(int i=0;i<nspin/2;++i) dest->site_id[i] = src->site_id[i];
        for(int i=0;i<ntype;++i){
            dest->type[i] = src->type[i];
            dest->weight[i] = src->weight[i];
        }

        dest->ngraph = src->ngraph;
    
        for(size_t i=0;i<size;++i){
            dest->graphs[i] = src->graphs[i];
            dest->tau[i] = src->tau[i];
            for(int j=0;j<nspin/2;++j) 
                dest->kink_id[i*nspin/2+j]=src->kink_id[i*nspin/2+j];
        }
    }
    else{
        printf("bond_memcpy : memory copy fail! Size of dest should be larger than size of src!\n");
        exit(-1);
    }
}

void bond_set_site_id(bond* bd, int nspin, const int* site_id){
    if(bd->nspin==nspin){
        for(int i=0;i<nspin/2;++i) bd->site_id[i] = site_id[i];
    }
    else{
        printf("bond_set_site_id : Set site id fail! nspin is not equal.\n");
        exit(-1);
    }
}

void bond_set_graph_type(bond* bd, int ntype, graph** type, const double* weight){
    if(bd->ntype==ntype){
        for(int i=0;i<ntype;++i){
            bd->type[i] = type[i];
            bd->weight[i] = weight[i];
        }
    }
    else{
        printf("bond_set_graph_type : Set graph type id fail! ntype is not equal.\n");
        exit(-1);
    }
}

int bond_insert_graph(bond* bd, int type, double tau, int* kink_id){
    int graph_id=-1;
    for(size_t i=0;i<bd->size;++i){
        if(bd->graphs[i]==-1){
            graph_id = i;
            break;
        }
    }

    if(graph_id!=-1){
        bd->graphs[graph_id] = type;
        bd->tau[graph_id] = tau;
        for(int i=0;i<bd->nspin/2;++i)
            bd->kink_id[bd->nspin/2*graph_id+i] = kink_id[i];

        bd->ngraph++;
    }
    else{
        printf("bond_insert_graph : Insert graph fail! Run out of storage.\n");
        exit(-1);
    }
    
    return graph_id;
}

void bond_remove_graph(bond* bd, int graph_id){
    if(bd->graphs[graph_id]!=-1){
        bd->graphs[graph_id] = -1;
        bd->tau[graph_id] = DBL_MAX;
        bd->kink_id[graph_id] = -1;
        bd->ngraph--;
    }
    else{
        printf("bond_remove_graph : Remove graph fail! Remove wrong graph_id.\n");
        exit(-1);
    }
}

#if 1
int main(int argc, char** argv){
    size_t size = 20;
    int nspin = 4;
    int ntype = 2;

    int site_id[2] = {0,1};
    graph* g[2] = {&GRAPH_DIAG,&GRAPH_HORI};
    double w[2] = {1.0,2.0};

    //printf("%d\n",rule_4spin(GRAPH_RULE_HORI,-1,1,1,-1));
    //printf("%d\n",rule_4spin(GRAPH_RULE_HORI,1,-1,1,-1));
    //printf("%d\n",rule_4spin(GRAPH_RULE_HORI,-1,1,-1,1));
    //printf("%d\n",rule_4spin(GRAPH_RULE_HORI,1,-1,-1,1));

    bond* bd1 = bond_alloc(size,nspin,ntype);
    bond* bd2 = bond_alloc(2048,nspin,ntype);

    bond_set_site_id(bd1,nspin,site_id);
    bond_set_graph_type(bd1,ntype,g,w);

    int kink_id[2] = {3,5};
    int graph_id = bond_insert_graph(bd1,0,3.14,kink_id);

    for(int i=0;i<size;++i){
        printf("%d %e \n",bd1->graphs[i],bd1->tau[i]);
    }

    bond_remove_graph(bd1,graph_id);

    for(int i=0;i<size;++i){
        printf("%d %e \n",bd1->graphs[i],bd1->tau[i]);
    }

    bond_memcpy(bd2,bd1);

    bond_free(bd1);
    bond_free(bd2);
}
#endif
