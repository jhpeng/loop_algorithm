#include "graph.h"

/*============= define 4-spin graph here ==============*/
//
//             ^
//             |       2|        |3
//             |        |********|
//             |        |--------|
//             |        |--------|
//    temporal |        |********|
//             |       0|        |1
//
//                    -------------->
//                     spatial
/*=====================================================*/
//diagonal graph  
int GRAPH_LINK_DIAG[4]  = {3,2,1,0};
int GRAPH_RULE_DIAG[16] = {1,0,0,0,
                           0,0,1,0,
                           0,1,0,0,
                           0,0,0,1};

graph GRAPH_DIAG = {4,GRAPH_LINK_DIAG,GRAPH_RULE_DIAG};

//horizontal graph  
int GRAPH_LINK_HORI[4]  = {1,0,3,2};
int GRAPH_RULE_HORI[16] = {0,0,0,0,
                           0,1,1,0,
                           0,1,1,0,
                           0,0,0,0};

graph GRAPH_HORI = {4,GRAPH_LINK_HORI,GRAPH_RULE_HORI};

int rule_4spin(int* rule, int* s){
    s[0] = (s[0]+1)/2;
    s[1] = (s[1]+1);
    s[2] = (s[2]+1);
    s[3] = (s[3]+1);
    return rule[s[0]+s[1]+2*s[2]+4*s[3]];
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
    assert(dest->size>=src->size);
    assert(dest->nspin==src->nspin);
    assert(dest->ntype==src->ntype);

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

void bond_set_site_id(bond* bd, int nspin, const int* site_id){
    assert(bd->nspin==nspin);

    for(int i=0;i<nspin/2;++i) bd->site_id[i] = site_id[i];
}

void bond_set_graph_type(bond* bd, int ntype, graph** type, const double* weight){
    assert(bd->ntype==ntype);

    for(int i=0;i<ntype;++i){
        bd->type[i] = type[i];
        bd->weight[i] = weight[i];
    }
}

int bond_check_available_id(bond* bd){
    int graph_id=-1;
    for(size_t i=0;i<bd->size;++i){
        if(bd->graphs[i]==-1){
            graph_id = i;
            break;
        }
    }

    assert(graph_id!=-1);

    return graph_id;
}

int bond_insert_graph(bond* bd, int type, double tau, int* kink_id){
    int graph_id=-1;
    for(size_t i=0;i<bd->size;++i){
        if(bd->graphs[i]==-1){
            graph_id = i;
            break;
        }
    }

    assert(graph_id!=-1);

    bd->graphs[graph_id] = type;
    bd->tau[graph_id] = tau;
    for(int i=0;i<bd->nspin/2;++i)
        bd->kink_id[bd->nspin/2*graph_id+i] = kink_id[i];

    bd->ngraph++;
    
    return graph_id;
}

void bond_remove_graph(bond* bd, int graph_id){
    assert(bd->graphs[graph_id]!=-1);
        
    bd->graphs[graph_id] = -1;
    bd->tau[graph_id] = DBL_MAX;
    bd->kink_id[graph_id] = -1;
    bd->ngraph--;
}

#if 0
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
