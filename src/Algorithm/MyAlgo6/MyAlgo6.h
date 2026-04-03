#ifndef __MYALGO6_H
#define __MYALGO6_H

#include    "../AlgorithmBase/AlgorithmBase.h"
#include    "../../Network/Graph/Graph.h"
#include    "../../config.h"

using namespace std;

class MyAlgo6 : public AlgorithmBase {
    vector<vector<vector<pair<int, int>>>> merge_shape;
    vector<vector<pair<int, int>>> recursion_build(int length);
    Shape_vector build_merge_shape(vector<int> path);

public:
    MyAlgo6(const Graph& graph, const vector<pair<int, int>>& requests, const map<SDpair, vector<Path>>& paths);
    void run();
};

#endif