#ifndef __ALGORITHMBASE_H
#define __ALGORITHMBASE_H

#include "../../Network/Graph/Graph.h"
#include "../../Network/PathMethod/PathMethodBase/PathMethod.h"
#include "../../config.h"
using namespace std;


class AlgorithmBase {
protected:
    string algorithm_name;
    Graph graph;
    map<string, double> res;
    vector<double> cdf;
    vector<SDpair> requests;
    int time_limit, memory_total, request_cnt;
    map<SDpair, vector<Path>> paths;
    void update_res();
    double A, B, n, T, tao;
    double bar(double F);
    double Fswap(double Fa, double Fb);
    double t2F(double t);
    double F2t(double F);
    double pass_tao(double F);
    vector<Path> get_paths(int src, int dst);
public:
    AlgorithmBase(const Graph& graph, const vector<SDpair>& requests, const map<SDpair, vector<Path>>& paths);
    map<string, double> get_res();
    double get_res(string str);
    vector<double> get_cdf();
    string get_name();
    virtual ~AlgorithmBase();
    virtual void run() = 0;
};

#endif