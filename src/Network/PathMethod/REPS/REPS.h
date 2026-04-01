#include "../PathMethodBase/PathMethod.h"
// #include "gurobi_c++.h"  // disabled: Gurobi not available on this system

using namespace std;

class REPS : public PathMethod {
    void PFT_LP(vector<double> &t_plum, vector<map<pair<int, int>, double>> &f_plum);
    pair<Path, double> dijkstra(int src, int dst, map<pair<int, int>, double>&f_plum_i);
    void update_f_plum(Path path, double width, map<pair<int, int>, double>&f_plum_i);
public:
    REPS();
    ~REPS();
    void build_paths(Graph graph, vector<SDpair> requests);
};