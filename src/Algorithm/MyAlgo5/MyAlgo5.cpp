#include "MyAlgo5.h"

MyAlgo5::MyAlgo5(const Graph& _graph, const vector<SDpair>& _requests, const map<SDpair, vector<Path>>& _paths):
    AlgorithmBase(_graph, _requests, _paths) {
    algorithm_name = "Linear";
    linear_shape.resize(graph.get_num_nodes() + 1);
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

vector<vector<pair<int, int>>> MyAlgo5::recursion_build(int length) {
    if(!linear_shape[length].empty()) return linear_shape[length];
    if(length == 1) {
        return linear_shape[length] = {};
    }
    if(length == 2) {
        return linear_shape[length] = {{{0, 1}}, {{0, 1}}};
    }

    vector<vector<pair<int, int>>> left = recursion_build(length - 1);
    vector<vector<pair<int, int>>> right = recursion_build(2);

    int offest = left.back().front().second - right.front().front().second;
    for(int i = 0; i < (int)right.size(); i++) {
        for(int j = 0; j < (int)right[i].size(); j++) {
            right[i][j].first += offest;
            right[i][j].second += offest;
        }
    }

    for(int i = 0; i < (int)left.size(); i++) {
        linear_shape[length].push_back(left[i]);
    }

    linear_shape[length].back().push_back(right.front().front());

    for(int i = 1; i < (int)right.size(); i++) {
        linear_shape[length].push_back(right[i]);
    }

    linear_shape[length].front().front().second++;
    linear_shape[length].back().front().second++;

    // cerr << "len = " << length << endl;
    // for(int i = 0; i < (int)linear_shape[length].size(); i++) {
    //     cerr << "-----" << endl;
    //     for(int j = 0; j < (int)linear_shape[length][i].size(); j++) {
    //         cerr << "{" << linear_shape[length][i][j].first << ", " << linear_shape[length][i][j].second << "}" << endl;
    //     }
    // }
    return linear_shape[length];
}
Shape_vector MyAlgo5::build_linear_shape(Path path) {
    vector<vector<pair<int, int>>> result = recursion_build(path.size());
    Shape_vector shape;
    for(int i = 0; i < (int)path.size(); i++) {
        shape.push_back({path[i], result[i]});
    }
    return shape;
}
void MyAlgo5::run() {
    vector<Shape_vector> shapes;
    for(int i = 0; i < (int)requests.size(); i++) {
        int src = requests[i].first;
        int dst = requests[i].second;

        vector<Path> paths = get_paths(src, dst);
        for(Path path : paths) {
            Shape_vector shape = build_linear_shape(path);
            bool cant = false;
            for(int i = 0; i < (int)shape.size(); i++) {
                for(int j = 0; j < (int)shape[i].second.size(); j++) {
                    if(shape[i].second[j].second >= graph.get_time_limit()) {
                        cant = true;
                    }
                }
            }

            if(!cant && graph.check_resource(shape, false)) {
                shapes.push_back(shape);
            }
            break;
        }
    }
    
    vector<pair<double, Shape_vector>> fidelity_shapes;

    for(Shape_vector shape : shapes) {
        double fidelity = Shape(shape).get_fidelity(A, B, n, T, tao, graph.get_F_init());
        fidelity_shapes.emplace_back(fidelity, shape);
    }

    sort(fidelity_shapes.begin(), fidelity_shapes.end());

    set<SDpair> used;
    for(auto P : fidelity_shapes) {
        Shape_vector shape = P.second;
        int src = shape[0].first, dst = shape.back().first;
        if(used.count({src, dst})) continue;
        used.insert({src, dst});
        if(graph.check_resource(shape, true)) {
            graph.reserve_shape(shape);
        }
    }
    update_res();
    cerr << "[" << algorithm_name << "] end" << endl;
}