#include "REPS.h"
#include <iostream>

using namespace std;

REPS::REPS() {
    method_name = "REPS";
}

REPS::~REPS() {}

void REPS::build_paths(Graph _graph, vector<SDpair> _requests) {
    cerr << "[REPS] Gurobi not available, REPS disabled." << endl;
    paths.clear();
    graph = _graph;
    requests = _requests;
}
