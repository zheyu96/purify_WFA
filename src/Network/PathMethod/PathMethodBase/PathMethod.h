#ifndef __PATHMETHOD_H
#define __PATHMETHOD_H

#include "../../Graph/Graph.h"
#include "../../../config.h"

using namespace std;
using Path = vector<int>;
using SDpair = pair<int, int>;


class PathMethod {
protected:
    map<SDpair, vector<Path>> paths;
    string method_name;
    Graph graph;
    vector<SDpair> requests;
    void reserve_path(Path path);
    void reserve_path(Path path, int amount);
public:
    PathMethod();
    ~PathMethod();
    string get_name();
    const vector<Path>& get_paths(int src, int dst);
    const map<SDpair, vector<Path>>& get_paths() const;
    virtual void build_paths(Graph _graph, vector<SDpair> _requests) = 0;
};

#endif