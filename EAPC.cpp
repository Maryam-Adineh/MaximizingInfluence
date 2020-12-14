#include <cstring>
#include <math.h>
#include "EAPC.h"


void EAPC::run(Graph *g, double ratio, double epsilon) {
    sprintf(file, "EAPC_%s_p=%.2f_e==%.4f.txt", g->getDataset().c_str(), ratio, epsilon);
    FILE *out = fopen(file, "w");
    n = g->getMaxNode();
    outputList = new int[setSize];
    activated = new bool[n]();
    levels = new int[n];
    memset(levels, -1, sizeof(int) * n);
    finalProb = new double[n]();
    prob = new double *[n];
    for (int i = 0; i < n; i++)
        prob[i] = new double[2]();

    double level = ceil(log(epsilon) / log(ratio));

    for (int k = 0; k < setSize; k++) {
        double max = -1000;
        int seed;
        for (int id = 0; id < n; id++)
            if (!activated[id]) {
                calculateProbs(g, id, ratio, level);

                double totalProb = 0;
                double last_totalProb = 0;
                while (!neighborhood.empty()) {
                    int node = neighborhood.front();
                    totalProb += prob[node][1];
                    last_totalProb += finalProb[node];
                    levels[node] = -1;
                    neighborhood.pop_front();
                }
                if (totalProb - last_totalProb > max) {
                    max = totalProb - last_totalProb;
                    seed = id;
                }
            }
        activated[seed] = true;
        outputList[k] = seed;
        fprintf(out, "%d\n", seed);

        calculateProbsFromSeed(g, seed, ratio, level);

        while (!neighborhood.empty()) {
            int node = neighborhood.front();
            levels[node] = -1;
            neighborhood.pop_front();
        }
    }

    for (int i = 0; i < n; i++) {
        delete[] prob[i];
    }

    delete[] finalProb;
    delete[] prob;
    delete[] levels;
    delete[] activated;

    fclose(out);
}

void EAPC::calculateProbs(Graph *g, int id, double rate, double level) {
    neighborhood.push_back(id);
    levels[id] = 0;
    prob[id][0] = finalProb[id];
    prob[id][1] = 1;
    int qcounter = 0;
    while (qcounter < neighborhood.size()) {
        int currentNode = neighborhood.at(qcounter);
        if (levels[currentNode] < level) {
            int numOfNeighbors = g->getNeighbor(currentNode);
            for (int i = 0; i < numOfNeighbors; i++) {
                Edge e = g->getEdge(currentNode, i);
                if (!activated[e.v] && levels[e.v] == -1) {
                    neighborhood.push_back(e.v);
                    levels[e.v] = levels[currentNode] + 1;
                    prob[e.v][0] = finalProb[e.v];
                    prob[e.v][1] = finalProb[e.v];
                }
                prob[e.v][1] =
                        1 - ((1 - prob[e.v][1]) * (1 - (prob[currentNode][1] * rate)) /
                             (1 - (prob[currentNode][0] * rate)));
            }
        }
        qcounter++;
    }
}

void EAPC::calculateProbsFromSeed(Graph *g, int id, double rate, double level) {
    neighborhood.push_back(id);
    levels[id] = 0;
    prob[id][0] = finalProb[id];
    prob[id][1] = 1;
    int qcounter = 0;
    while (qcounter < neighborhood.size()) {
        int currentNode = neighborhood.at(qcounter);
        if (levels[currentNode] < level) {
            int numOfNeighbors = g->getNeighbor(currentNode);
            for (int i = 0; i < numOfNeighbors; i++) {
                Edge e = g->getEdge(currentNode, i);
                if (!activated[e.v] && levels[e.v] == -1) {
                    neighborhood.push_back(e.v);
                    levels[e.v] = levels[currentNode] + 1;
                    prob[e.v][0] = finalProb[e.v];
                    prob[e.v][1] = finalProb[e.v];
                }
                prob[e.v][1] =
                        1 - ((1 - prob[e.v][1]) * (1 - (prob[currentNode][1] * rate)) /
                             (1 - (prob[currentNode][0] * rate)));
            }
        }
        finalProb[currentNode] = prob[currentNode][1];
        qcounter++;
    }
}

void EAPC::run_celf_optimized(Graph *g, double ratio, double epsilon) {
    sprintf(file, "EAPC_%s_p=%.2f_e=%.3f.txt", g->getDataset().c_str(), ratio, epsilon);
    FILE *out = fopen(file, "w");
    n = g->getMaxNode();

    prob = new double *[n];
    for (int i = 0; i < n; i++) {
        prob[i] = new double[2]();
    }

    vector<double> improve(n + 1);
    vector<int> lastupdate(n);
    vector<int> heap(n);
    activated = new bool[n]();
    outputList = new int[setSize];
    finalProb = new double[n]();
    levels = new int[n];
    memset(levels, -1, sizeof(int) * n);

    double level = ceil(log(epsilon) / log(ratio));

    for (int i = 0; i < n; i++) {
        heap[i] = i;
        lastupdate[i] = -1;
        improve[i] = (double) (n + 1);
    }

    for (int i = 0; i < setSize; i++) {
        while (lastupdate[heap[0]] != i) {
            lastupdate[heap[0]] = i;
            activated[heap[0]] = true;
            improve[heap[0]] = marginal_gain_BFS(g, heap[0], ratio, level);
            activated[heap[0]] = false;

            int x = 0;
            while (x * 2 + 2 <= n - i) {
                int newx = x * 2 + 1;
                if ((newx + 1 < n - i) && (improve[heap[newx]] < improve[heap[newx + 1]]))
                    newx++;
                if (improve[heap[x]] < improve[heap[newx]]) {
                    int t = heap[x];
                    heap[x] = heap[newx];
                    heap[newx] = t;
                    x = newx;
                } else
                    break;
            }
        }

        activated[heap[0]] = true;
        outputList[i] = heap[0];
        fprintf(out, "%d\n", outputList[i]);

        calculateProbsFromSeed(g, outputList[i], ratio, level);

        while (!neighborhood.empty()) {
            int node = neighborhood.front();
            levels[node] = -1;
            neighborhood.pop_front();
        }

        heap[0] = heap[n - i - 1];
        int x = 0;

        while (x * 2 + 2 <= n - i) {
            int newx = x * 2 + 1;
            if ((newx + 1 < n - i) && (improve[heap[newx]] < improve[heap[newx + 1]]))
                newx++;
            if (improve[heap[x]] < improve[heap[newx]]) {
                int t = heap[x];
                heap[x] = heap[newx];
                heap[newx] = t;
                x = newx;
            } else
                break;
        }
    }

    fclose(out);

    for (int i = 0; i < n; i++)
        delete[] prob[i];

    delete[] prob;
    delete[] activated;
    delete[] levels;
    delete[] finalProb;
}

double EAPC::marginal_gain_BFS(Graph *g, int i, double rate, double level) {
    double totalProb = 0;
    double last_totalProb = 0;

    calculateProbs(g, i, rate, level);

    while (!neighborhood.empty()) {
        int node = neighborhood.front();
        totalProb += prob[node][1];
        last_totalProb += finalProb[node];
        levels[node] = -1;
        neighborhood.pop_front();
    }

    return totalProb - last_totalProb;
}


int* EAPC::getSet() {
    return outputList;
}