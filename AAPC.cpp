#include "AAPC.h"

using namespace std;

void AAPC::run(Graph *g, double ratio, int step) {
    sprintf(file, "AAPC_%s_p=%.2f_step=%d.txt", g->getDataset().c_str(), ratio, step);
    FILE *output = fopen(file, "w");
    n = g->getMaxNode();

    p1 = new double *[n];
    for (int i = 0; i < n; i++)
        p1[i] = new double[step]();

    p2 = new double *[n];
    for (int i = 0; i < n; i++)
        p2[i] = new double[step]();

    activated = new bool[n]();
    outputList = new int[setSize];

    for (int i = 0; i < n; i++) {
        p1[i][0] = 0;
        p2[i][0] = 0;
    }

    for (int s = 0; s < setSize; s++) {
        int seed;
        double maxInf = 0;
        double totalProb;
        for (int id = 0; id < n; id++)
            if (!activated[id]) {
                activated[id] = true;
                for (int i = 1; i < step; i++) {
                    for (int j = 0; j < n; j++)
                        if (activated[j])
                            p1[j][i] = 0;
                        else
                            p1[j][i] = 1;
                    for (int j = 0; j < n; j++) {
                        int numOfNeighbors = g->getNeighbor(j);
                        for (int k = 0; k < numOfNeighbors; k++) {
                            Edge e = g->getEdge(j, k);
                            p1[e.v][i] *= (1 - p1[j][i - 1] * ratio);
                        }
                    }
                    for (int j = 0; j < n; j++)
                        p1[j][i] = (1 - p1[j][i]) * (1 - p2[j][i - 1]);
                    for (int j = 0; j < n; j++)
                        p2[j][i] = 1 - (1 - p2[j][i - 1]) * (1 - p1[j][i]);
                }

                activated[id] = false;
                totalProb = 0;
                for (int j = 0; j < n; j++)
                    totalProb += p2[j][step - 1];
                if (totalProb > maxInf) {
                    maxInf = totalProb;
                    seed = id;
                }
            }
        activated[seed] = true;
        p1[seed][0] = 1;
        p2[seed][0] = 1;
        outputList[s] = seed;
        fprintf(output, "%d\n", seed);
    }

    delete[] activated;
    for (int i = 0; i < n; i++)
        delete[]p1[i];
    delete[] p1;

    for (int i = 0; i < n; i++)
        delete[]p2[i];
    delete[] p2;

    fclose(output);
}

void AAPC::run_celf_optimized(Graph *g, double ratio, int step) {
    sprintf(file, "AAPC_%s_p=%.2f_Step=%d.txt", g->getDataset().c_str(), ratio, step);
    FILE *out = fopen(file, "w");
    n = g->getMaxNode();

    double old = 0.0;
    activated = new bool[n]();

    p1 = new double *[n];
    for (int i = 0; i < n; i++)
        p1[i] = new double[step]();

    p2 = new double *[n]();
    for (int i = 0; i < n; i++)
        p2[i] = new double[step];

    vector<double> improve(n + 1);
    vector<int> lastupdate(n);
    vector<int> heap(n);
    outputList = new int[setSize];

    for (int i = 0; i < n; i++) {
        p1[i][0] = 0;
        p2[i][0] = 0;
    }

    for (int i = 0; i < n; i++) {
        heap[i] = i;
        lastupdate[i] = -1;
        improve[i] = (double) (n + 1);
    }

    for (int i = 0; i < setSize; i++) {
        while (lastupdate[heap[0]] != i) {
            lastupdate[heap[0]] = i;
            activated[heap[0]] = true;
            improve[heap[0]] = calculate_probabilities(g, ratio, step) - old;
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
        p1[outputList[i]][0] = 1;
        p2[outputList[i]][0] = 1;
        old += improve[heap[0]];

        fprintf(out, "%d\n", outputList[i]);

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

    delete[] activated;

    for (int i = 0; i < n; i++)
        delete[] p1[i];
    delete[]p1;

    for (int i = 0; i < n; i++)
        delete[]p2[i];
    delete[] p2;

    fclose(out);
}

double AAPC::calculate_probabilities(Graph *g, double ratio, int step) {
    for (int i = 1; i < step; i++) {
        for (int j = 0; j < n; j++)
            if (activated[j])
                p1[j][i] = 0;
            else
                p1[j][i] = 1;

        for (int j = 0; j < n; j++) {
            int numOfNeighbors = g->getNeighbor(j);
            for (int k = 0; k < numOfNeighbors; k++) {
                Edge e = g->getEdge(j, k);
                p1[e.v][i] *= (1 - p1[j][i - 1] * ratio);
            }
        }

        for (int j = 0; j < n; j++)
            p1[j][i] = (1 - p1[j][i]) * (1 - p2[j][i - 1]);

        for (int j = 0; j < n; j++)
            p2[j][i] = 1 - (1 - p2[j][i - 1]) * (1 - p1[j][i]);
    }

    double totalProb = 0;
    for (int j = 0; j < n; j++)
        totalProb += p2[j][step - 1];

    return totalProb;
}

int *AAPC::getSet() {
    return outputList;
}

