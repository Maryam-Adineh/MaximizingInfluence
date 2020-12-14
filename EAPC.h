#ifndef EAPC_H
#define EAPC_H

#include <deque>
#include "limit.h"
#include "Graph.h"

class EAPC {
private:
    int n;
    int *outputList;
    char file[STR_LEN];
    bool * activated;
    deque<int> neighborhood;
    int *levels;
    double ** prob;
    double * finalProb;
public:
    void run(Graph *g, double ratio, double epsilon);
    void run_celf_optimized(Graph *g, double ratio, double epsilon);
    void calculateProbs(Graph *g, int id, double rate, double level);
    void calculateProbsFromSeed(Graph *g, int id, double rate, double level);
    double marginal_gain_BFS(Graph *g, int i, double rate, double level);
    int* getSet();
};


#endif //EAPC_H
