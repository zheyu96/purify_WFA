#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <random>
#include <algorithm>
#include <chrono> 
#include <omp.h>   
#include <cassert>
#include <ctime>

#include "./config.h"
#include "Network/Graph/Graph.h"
#include "Algorithm/AlgorithmBase/AlgorithmBase.h"
#include "Algorithm/WernerAlgo2_time/WernerAlgo2_time.h" 
#include "Algorithm/WernerAlgo_time/WernerAlgo_time.h"
#include "Network/PathMethod/Greedy/Greedy.h"

using namespace std;

// --- 輔助函式：產生隨機請求 ---
SDpair generate_new_request(int num_of_node){
    random_device rd;
    default_random_engine generator = default_random_engine(rd());
    uniform_int_distribution<int> unif(0, num_of_node-1);
    int node1 = unif(generator), node2 = unif(generator);
    while(node1 == node2) node2 = unif(generator);
    return make_pair(node1, node2);
}

// --- 輔助函式：根據 Fidelity 門檻產生請求 ---
vector<SDpair> generate_requests_fid(Graph graph, int requests_cnt, double fid_th, double hop_th) {
    int n = graph.get_num_nodes();
    vector<pair<SDpair, double>> cand[22];
    random_device rd;
    default_random_engine generator = default_random_engine(rd());
    uniform_int_distribution<int> unif(0, 1e9);
    
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if(i == j) continue;
            double fid = graph.get_ini_fid(i, j);
            if(fid > fid_th && graph.distance(i, j) >= hop_th) {
                int index = fid / 0.05;
                if(index < 0) index = 0;
                if(index > 20) index = 20;
                cand[index].emplace_back(make_pair(make_pair(i, j), (double)graph.distance(i, j)));
            }
        }
    }

    for(int i = 0; i < 22; i++) {
        shuffle(cand[i].begin(), cand[i].end(), generator);
    }

    vector<SDpair> requests;
    int pos[22] = {0};
    int idx = 0;
    while(requests.size() < requests_cnt) {
        int cnt = unif(generator) % 5 + 4;
        cnt = min(cnt, (int)(requests_cnt - requests.size()));
        int try_idx = 0;
        while(cand[21 - idx].empty() && try_idx < 22) {
            idx = (idx + 1) % 22;
            try_idx++;
        }
        if(try_idx >= 22) break; 

        for(int i = 0; i < cnt && pos[21-idx] < cand[21-idx].size(); i++) {
            requests.push_back(cand[21 - idx][pos[21 - idx]].first);
            pos[21 - idx]++;
        }
        idx = (idx + 1) % 22;
    }
    return requests;
}

int main() {
    string file_path = "../data/";
    system("mkdir -p ../data/ans");
    system("mkdir -p ../data/log");
    system("mkdir -p ../data/input");

    // 1. 設定預設參數
    map<string, double> default_setting;
    default_setting["num_nodes"] = 100;
    default_setting["request_cnt"] = 50;
    default_setting["time_limit"] = 13;
    default_setting["avg_memory"] = 10;
    default_setting["tao"] = 0.002;
    default_setting["fidelity_threshold"] = 0.7;
    default_setting["epsilon"] = 0.75; 
    default_setting["Zmin"] = 0.02702867239;
    default_setting["bucket_eps"] = 0.01;
    default_setting["time_eta"] = 0.001;
    
    // 2. 設定實驗變量 (X 軸)
    map<string, vector<double>> change_parameter;
    change_parameter["epsilon"] = {0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95}; 
    change_parameter["bucket_eps"] = {0.001, 0.01, 0.05, 0.1, 0.5};

    vector<string> X_names = {"epsilon", "bucket_eps"}; 
    vector<string> Y_names = {"fidelity_gain", "succ_request_cnt", "runtime"};
    int round = 5; // 測試時可設小，正式實驗建議 50

    // 3. 預先產生固定請求 (確保對照組公平性)
    vector<vector<SDpair>> default_requests(round);
    cerr << "Pre-generating requests..." << endl;
    for(int r = 0; r < round; r++) {
        string filename = file_path + "input/round_" + to_string(r) + ".input";
        ifstream f(filename);
        if(!f.good()){
            string command = "python3 graph_generator.py " + filename + " " + to_string((int)default_setting["num_nodes"]);
            system(command.c_str());
        }
        Graph temp_g(filename, 20, 0.9, 10, 0.5, 0.98, 0.7, 0.25, 0.75, 2, 10, 0.002, 0.027, 0.01, 0.001);
        default_requests[r] = generate_requests_fid(temp_g, (int)default_setting["request_cnt"], 0.7, 3);
        assert(!default_requests[r].empty());
    }

    // 4. 開始實驗
    for (string X_name : X_names) {
        for(double change_value : change_parameter[X_name]) {
            cerr << "Testing " << X_name << " = " << change_value << endl;
            // result[round][algo_name][metric_name]
            vector<map<string, map<string, double>>> result(round);

            #pragma omp parallel for
            for(int r = 0; r < round; r++) {
                string filename = file_path + "input/round_" + to_string(r) + ".input";
                map<string, double> input_param = default_setting;
                input_param[X_name] = change_value;

                Graph graph(filename, (int)input_param["time_limit"], 0.9, (int)input_param["avg_memory"], 
                            0.5, 0.98, input_param["fidelity_threshold"], 0.25, 0.75, 2, 10, 
                            input_param["tao"], input_param["Zmin"], 
                            input_param["bucket_eps"], input_param["time_eta"]);

                Greedy path_method;
                path_method.build_paths(graph, default_requests[r]);
                const auto& paths = path_method.get_paths();

                // --- 執行 WernerAlgo2_time (ZFA2) ---
                AlgorithmBase* algo2 = new WernerAlgo2_time(graph, default_requests[r], paths, 
                                                           input_param["epsilon"], input_param["bucket_eps"]);
                auto s2 = chrono::high_resolution_clock::now();
                algo2->run();
                auto e2 = chrono::high_resolution_clock::now();
                
                string name2 = algo2->get_name(); 
                #pragma omp critical
                {
                    result[r][name2]["fidelity_gain"] = algo2->get_res("fidelity_gain");
                    result[r][name2]["succ_request_cnt"] = algo2->get_res("succ_request_cnt");
                    result[r][name2]["runtime"] = chrono::duration<double>(e2 - s2).count();
                }
                delete algo2;

                // --- 執行 WernerAlgo_time (ZFA) ---
                AlgorithmBase* algo1 = new WernerAlgo_time(graph, default_requests[r], paths, input_param["epsilon"]);
                auto s1 = chrono::high_resolution_clock::now();
                algo1->run();
                auto e1 = chrono::high_resolution_clock::now();

                string name1 = algo1->get_name(); 
                #pragma omp critical
                {
                    result[r][name1]["fidelity_gain"] = algo1->get_res("fidelity_gain");
                    result[r][name1]["succ_request_cnt"] = algo1->get_res("succ_request_cnt");
                    result[r][name1]["runtime"] = chrono::duration<double>(e1 - s1).count();
                }
                delete algo1;
            }

            // 5. 寫入結果
            // 假設 get_name() 回傳的是 "ZFA_time" 和 "ZFA2_time"
            vector<string> target_algos = {"ZFA_time", "ZFA2_time"};
            
            for(string a_name : target_algos) {
                for(string Y_name : Y_names) {
                    string ans_file = file_path + "ans/Greedy_" + X_name + "_" + Y_name + "_" + a_name + ".ans";
                    ofstream ofs(ans_file, ios::app);
                    
                    double sum = 0;
                    int count = 0;
                    for(int r = 0; r < round; r++) {
                        if(result[r].count(a_name)) {
                            sum += result[r][a_name][Y_name];
                            count++;
                        }
                    }
                    
                    if(count > 0) {
                        ofs << change_value << " " << sum / count << endl;
                    }
                    ofs.close();
                }
            }
        }
    }
    return 0;
}