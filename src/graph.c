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

