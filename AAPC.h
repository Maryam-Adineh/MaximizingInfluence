#ifndef AAPC_H
#define AAPC_H

#include "Graph.h"
#include "limit.h"

class AAPC {
private:
    int n;
    int *outputList;
    char file[STR_LEN];
    bool * activated;
    double ** p1;
    double ** p2;
public:
    void run(Graph *graph, double ratio, int step);
    void run_celf_optimized(Graph *graph, double ratio, int step);
    double calculate_probabilities(Graph *g, double ratio, int step);
    int *getSet();
};


#endif //AAPC_H
