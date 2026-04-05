#include "PathMethod.h"

PathMethod::PathMethod() {}
PathMethod::~PathMethod() {}

const vector<Path>& PathMethod::get_paths(int src, int dst) {
    return paths[{src, dst}];
}

const map<SDpair, vector<Path>>& PathMethod::get_paths() const {
    return paths;
}

string PathMethod::get_name() {
    return method_name;
}

void PathMethod::reserve_path(Path path) {
    graph.reserve_path(path);
    paths[{path[0], path.back()}].push_back(path);
    if(paths[{path[0], path.back()}].size() >= 1) {
        for(int i = 0; i < (int)requests.size(); i++) {
            if(requests[i] == make_pair(path[0], path.back())) {
                requests.erase(requests.begin() + i);
                return;
            }
        }
    }
}

void PathMethod::reserve_path(Path path, int amount) {
    assert(graph.check_path_resource(path, amount));
    graph.reserve_path(path, amount);
    paths[{path[0], path.back()}].push_back(path);
    if(paths[{path[0], path.back()}].size() >= 1) {
        for(int i = 0; i < (int)requests.size(); i++) {
            if(requests[i] == make_pair(path[0], path.back())) {
                requests.erase(requests.begin() + i);
                return;
            }
        }
    }
}