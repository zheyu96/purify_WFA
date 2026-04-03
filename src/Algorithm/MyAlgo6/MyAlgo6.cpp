#include "MyAlgo6.h"

MyAlgo6::MyAlgo6(const Graph& _graph, const vector<SDpair>& _requests, const map<SDpair, vector<Path>>& _paths):
    AlgorithmBase(_graph, _requests, _paths) {
    algorithm_name = "ASAP";
    merge_shape.resize(graph.get_num_nodes() + 1);
    for(SDpair sdpair : requests) {
        sort(paths[sdpair].begin(), paths[sdpair].end(), [&](Path &a, Path &b) {
            return graph.path_Pr(a) < graph.path_Pr(b);
        });
    }
    // sort(requests.begin(), requests.end(), [&](SDpair &a, SDpair &b) {
    //     double Pr_a = (paths[a].empty() ? 1 : graph.path_Pr(paths[a][0]));
    //     double Pr_b = (paths[b].empty() ? 1 : graph.path_Pr(paths[b][0]));
    //     return Pr_a < Pr_b;
    // });
}

vector<vector<pair<int, int>>> MyAlgo6::recursion_build(int length) {
    if(!merge_shape[length].empty()) return merge_shape[length];
    if(length == 1) {
        return merge_shape[length] = {};
    }
    if(length == 2) {
        return merge_shape[length] = {{{0, 1}}, {{0, 1}}};
    }

    int right_size = (length + 1) / 2, left_size = length - right_size + 1;
    vector<vector<pair<int, int>>> left = recursion_build(left_size);
    vector<vector<pair<int, int>>> right = recursion_build(right_size);

    int offest = left.back().front().second - right.front().front().second;
    for(int i = 0; i < (int)right.size(); i++) {
        for(int j = 0; j < (int)right[i].size(); j++) {
            right[i][j].first += offest;
            right[i][j].second += offest;
        }
    }

    for(int i = 0; i < (int)left.size(); i++) {
        merge_shape[length].push_back(left[i]);
    }

    merge_shape[length].back().push_back(right.front().front());

    for(int i = 1; i < (int)right.size(); i++) {
        merge_shape[length].push_back(right[i]);
    }

    merge_shape[length].front().front().second++;
    merge_shape[length].back().front().second++;

    // cerr << "len = " << length << endl;
    // for(int i = 0; i < (int)merge_shape[length].size(); i++) {
    //     cerr << "-----" << endl;
    //     for(int j = 0; j < (int)merge_shape[length][i].size(); j++) {
    //         cerr << "{" << merge_shape[length][i][j].first << ", " << merge_shape[length][i][j].second << "}" << endl;
    //     }
    // }
    return merge_shape[length];
}
Shape_vector MyAlgo6::build_merge_shape(vector<int> path) {
    vector<vector<pair<int, int>>> result = recursion_build(path.size());
    Shape_vector shape;
    for(int i = 0; i < (int)path.size(); i++) {
        shape.push_back({path[i], result[i]});
    }
    return shape;
}
void MyAlgo6::run() {
    for(int i = 0; i < (int)requests.size(); i++) {
        int src = requests[i].first;
        int dst = requests[i].second;
        vector<Path> paths = get_paths(src, dst);
        for(Path path : paths) {
            Shape_vector shape = build_merge_shape(path);
            bool cant = false;
            for(int i = 0; i < (int)shape.size(); i++) {
                for(int j = 0; j < (int)shape[i].second.size(); j++) {
                    if(shape[i].second[j].second >= graph.get_time_limit()) {
                        cant = true;
                    }
                }
            }
            if(!cant && graph.check_resource_ASAP(shape, false)) {
                graph.reserve_shape_ASAP(shape);
            }
            break;
        }
    }

    update_res();
    cerr << "[" << algorithm_name << "] end" << endl;
}