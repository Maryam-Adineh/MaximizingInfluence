#include <iostream>
#include <cstring>
#include <cmath>
#include "limit.h"
#include "Graph.h"
#include "IndependCascade.h"
#include "AAPC.h"
#include "EAPC.h"

using namespace std;
double timer;
clock_t start;

int setSize = 0;

void toSimulate(char *file, int *sett, double (*Run)(int set[], int size, int iter),
        int iteration) {
    FILE *out = fopen(file, "w");
    int *set = new int[setSize];
    for (int t = 0; t < setSize; t++) {
        set[t] = sett[t];
        fprintf(out, "%02d \t %10lg\n", t + 1, Run(set, t + 1, iteration));
    }
    fclose(out);
    delete set;
}

int main(int argc, char **argv) {
    Graph graph;
    string dataset;
    char timeFile[STR_LEN], infFile[STR_LEN];
    FILE *timeFile_p;
    double rate, epsilon;
    int step, iteration;
    bool direction = true;
    if (argc != 7) {
        cout << "USAGE:\n"
             << " ./InfluenceMaximization <AAPC|EAPC> <dataset> <d|u> <p> <k> <R>";
        return 1;
    }

    if (!strncmp(argv[3], "u", 1))
        direction = false;
    rate = atof(argv[4]);
    setSize = atoi(argv[5]);
    iteration = atoi(argv[6]);

    if (!graph.build(argv[2], direction)) {
        cerr << "Failed to build graph" << endl;
        return 1;
    }

    time_t time1 = time(NULL);
    cout << "*************************************************\n";
    cout << "current time: " << ctime(&time1);
    cout << "dataset: " << graph.getDataset().c_str() << "\trate: " << rate << endl;
    cout << "*************************************************\n";
    IndependCascade ic(&graph, rate);
    graph.removeMultiplicity();

    if (!strncmp(argv[1],"AAPC",strlen(argv[1]))) {
        AAPC aapc;
        step = 4;
        start = clock();
        aapc.run_celf_optimized(&graph, rate, step);
        timer = (double) (clock() - start) / CLOCKS_PER_SEC;
        sprintf(timeFile, "AAPC_time_%s_p=%.2f_step=%d.txt",
                graph.getDataset().c_str(), rate, step);
        timeFile_p = fopen(timeFile, "w");
        fprintf(timeFile_p, "%f", timer);
        fclose(timeFile_p);
        sprintf(infFile, "AAPC_IC_%s_p=%.2f_step=%d.txt",
                graph.getDataset().c_str(), rate, step);
        cout << "Calculating influence spread is started\n";
        toSimulate(infFile, aapc.getSet(), IndependCascade::run, iteration);
    } else if (!strncmp(argv[1], "EAPC", strlen(argv[1]))) {
        EAPC eapc;
        epsilon = 0.01;
        start = clock();
        eapc.run(&graph, rate, epsilon);
        timer = (double) (clock() - start) / CLOCKS_PER_SEC;
        sprintf(timeFile, "EAPC_time_%s_p=%.2f_e=%.4f.txt",
                graph.getDataset().c_str(), rate, epsilon);
        timeFile_p = fopen(timeFile, "w");
        fprintf(timeFile_p, "%f", timer);
        fclose(timeFile_p);
        sprintf(infFile, "EAPC_IC_%s_p=%.2f_e=%.4f.txt",
                graph.getDataset().c_str(), rate, epsilon);
        cout << "Calculating influence spread is started\n";
        toSimulate(infFile, eapc.getSet(), IndependCascade::run, iteration);
    } else {
        cout << "Error: Unknown algorithm\n";
        return 1;
    }

    cout << "Finished successfully\n";
    return 0;
}