
#include "./config.h"
#include "Network/Graph/Graph.h"
#include "Algorithm/AlgorithmBase/AlgorithmBase.h"
#include "Algorithm/MyAlgo1/MyAlgo1.h"
#include "Algorithm/MyAlgo2/MyAlgo2.h"
#include "Algorithm/MyAlgo3/MyAlgo3.h"
#include "Algorithm/MyAlgo4/MyAlgo4.h"
#include "Algorithm/MyAlgo5/MyAlgo5.h"
#include "Algorithm/MyAlgo6/MyAlgo6.h"
#include "Algorithm/WernerAlgo/WernerAlgo.h"
#include "Algorithm/WernerAlgo2/WernerAlgo2.h"
#include "Algorithm/WernerAlgo3/WernerAlgo3.h"
#include "Algorithm/WernerAlgo_UB/WernerAlgo_UB.h"
#include "Network/PathMethod/PathMethodBase/PathMethod.h"
#include "Network/PathMethod/Greedy/Greedy.h"
#include "Network/PathMethod/QCAST/QCAST.h"
#include "Network/PathMethod/REPS/REPS.h"

using namespace std;

SDpair generate_new_request(int num_of_node){
    random_device rd;
    default_random_engine generator = default_random_engine(rd());
    uniform_int_distribution<int> unif(0, num_of_node-1);
    int node1 = unif(generator), node2 = unif(generator);
    while(node1 == node2) node2 = unif(generator);

    return make_pair(node1, node2);
}

vector<SDpair> generate_requests(Graph &graph, int requests_cnt, int length_lower, int length_upper) {
    int n = graph.get_num_nodes();
    vector<SDpair> cand;
    random_device rd;
    default_random_engine generator = default_random_engine(rd());
    uniform_int_distribution<int> unif(0, 1e9);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if(i == j) continue;
            int dist = graph.distance(i, j);
            if(dist >= length_lower && dist <= length_upper) {
                cand.emplace_back(i, j);
            }
        }
    }

    random_shuffle(cand.begin(), cand.end());

    vector<SDpair> requests;
    for(SDpair sdpair : cand) {
        int cnt = unif(generator) % 4 + 3;
        while(cnt--) requests.push_back(sdpair);
    }

    while((int)requests.size() < requests_cnt) {
        requests.emplace_back(generate_new_request(n));
    }

    while((int)requests.size() > requests_cnt) {
        requests.pop_back();
    }

    return requests;
}
vector<SDpair> generate_requests_fid(Graph &graph, int requests_cnt,double fid_th,double hop_th, double fid_upper = 1) {
    int n = graph.get_num_nodes();
    vector<pair<SDpair,double>> cand[22];
    random_device rd;
    default_random_engine generator = default_random_engine(rd());
    uniform_int_distribution<int> unif(0, 1e9);
    int sd_cnt=0;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if(i == j) continue;
            double fid = graph.get_ini_fid(i,j);
            //cerr<<"fid of "<<i<<" "<<j<<" : "<<fid<<endl;
            assert(fid>=0.0&&fid<=1.0);
            if(fid > fid_th && fid <= fid_upper && graph.distance(i,j)>=hop_th) {
                int index = fid/0.05;
                //index-=5;
                if(index < 0) continue;
                if(index > 20) index = 20;
                int d=graph.distance(i, j),f0=fid,prob=pow(0.1,d)*pow(0.9,max(d-1,0));
                //double score = f0+prob*100-0.1*d;
                cand[index].emplace_back(std::make_pair(std::make_pair(i, j), graph.distance(i, j)));
                if(graph.distance(i,j)>=1)sd_cnt++;
            }
        }
    }
     cerr << "\033[1;32m"<< "[SD ini pairs] : "<<sd_cnt<< "\033[0m"<< endl;
    /*for(int i=21;i>=0;i--){
        if(!cand[i].empty()){
            random_shuffle(cand[i].begin(), cand[i].end());
        }
    } */
    /* for(int i=21;i>=0;i--){
        sort(cand[i].begin(),cand[i].end(),[](const pair<SDpair,double>& L,const pair<SDpair,double>& R){
            return L.second > R.second;
        }) ;
    }  */
    for(int i=0;i<22;i++){
        random_shuffle(cand[i].begin(), cand[i].end());
    }
    // 檢查是否有任何候選
    bool any_cand = false;
    for (int i = 0; i < 22; i++) if (!cand[i].empty()) any_cand = true;
    if (!any_cand) {
        cerr << "[generate_requests_fid] WARNING: no candidates found (fid_th=" << fid_th << ", hop_th=" << hop_th << ")" << endl;
        return {};
    }

    vector<SDpair> requests;
    int pos[22];
    for(int i=0;i<22;i++) pos[i]=0;
    int idx=0;
    while((int)requests.size()<requests_cnt){
        int cnt=unif(generator) % 5 +4;
        cnt=min(cnt,(int)(requests_cnt-(int)requests.size()));
        // 找下一個非空桶（有保護）
        int tries = 0;
        while(cand[21-idx].empty()){
            idx++;
            if(idx>=22) idx=0;
            if(++tries > 22) break;  // 防止無限迴圈
        }
        if(tries > 22) break;
        if(!cand[21-idx].empty()){
            for(int i=0;i<cnt;i++){
                requests.push_back(cand[21-idx][pos[21-idx]].first);
            }
            pos[21-idx]++;
            pos[21-idx]%=cand[21-idx].size();
        }
        idx=(idx+1)%22;
    }
    if ((int)requests.size() < requests_cnt)
        cerr << "[generate_requests_fid] only generated " << requests.size() << "/" << requests_cnt << " requests" << endl;
    return requests;
}
// 生成「purification 能帶來優勢」的 request
// 三個優先級（由高到低）：
//   A) 嚴格甜蜜點：w_no < w_th AND w_pur >= w_th（不做 purify 過不了，做了就能過）
//   B) 邊緣受益者：w_no 只比 w_th 高一點（< w_th * margin_ratio），purify 後 w_pur 遠高於 w_th
//      → 這些 pair 不做 purify 勉強過，但 fidelity 很低，做 purify 後 Pr*w 大幅提升
//   C) 必須 purify 的長 hop（補充用）
vector<SDpair> generate_requests_purify_needed(Graph &graph, int requests_cnt, int min_hop = 2) {
    int n = graph.get_num_nodes();
    double fid_th = graph.get_fidelity_threshold();
    double w_th = (4.0 * fid_th - 1.0) / 3.0;

    // BFS 找最短路徑
    auto bfs_path = [&](int src, int dst) -> vector<int> {
        vector<int> parent(n, -1);
        vector<bool> vis(n, false);
        queue<int> que;
        vis[src] = true;
        que.push(src);
        while (!que.empty()) {
            int u = que.front(); que.pop();
            if (u == dst) break;
            for (int v : graph.adj_list[u]) {
                if (!vis[v]) {
                    vis[v] = true;
                    parent[v] = u;
                    que.push(v);
                }
            }
        }
        if (!vis[dst]) return {};
        vector<int> path;
        for (int v = dst; v != -1; v = parent[v]) path.push_back(v);
        reverse(path.begin(), path.end());
        return path;
    };

    // Eq.8: r rounds pumping purification
    auto purified_werner = [](double w_e, int rounds) -> double {
        double w_cur = w_e;
        for (int r = 0; r < rounds; r++) {
            w_cur = (3.0*w_cur*w_e + 3.0*w_cur + 3.0*w_e - 1.0)
                  / (9.0*w_cur*w_e - 3.0*w_cur - 3.0*w_e + 5.0);
        }
        return w_cur;
    };

    // 考慮時間衰退：W-domain 中每個 edge 加上 decoherence
    double eta = graph.get_tao() / graph.get_time_limit();  // per-slot η in code's convention
    double kappa = graph.get_n();  // κ

    // 邊緣受益者的判定：w_no 超過 w_th 但不超過 margin_ratio 倍
    // 這些 pair 不做 purify 勉強過，但 fidelity 低，做 purify 後 Pr*w 提升顯著
    const double margin_ratio = 1.08;  // w_no < w_th * 1.08 就算「邊緣」

    // 優先級：A=嚴格甜蜜點(score 2.0+), B=邊緣受益者(score 1.0+)
    vector<pair<double, SDpair>> candidates;

    // 診斷用
    struct HopDiag { int total=0, pass_no=0, marginal=0, sweet=0, fail_both=0; double sum_w=0; };
    map<int, HopDiag> diag;

    double tao_val = graph.get_tao();
    double T_param = graph.get_T();
    double n_param = graph.get_n();
    double w_decay_per_slot = exp(-pow(tao_val / T_param, n_param));

    const int max_purify_rounds = 3;  // 最多考慮 3 rounds purification

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            vector<int> path = bfs_path(i, j);
            if (path.empty()) continue;
            int h = (int)path.size() - 1;
            if (h < min_hop) continue;

            // 收集每條 edge 的 w_e
            vector<double> w_edges(h);
            double w_no_purify = 1.0;
            for (int k = 0; k < h; k++) {
                w_edges[k] = graph.get_link_werner(path[k], path[k+1]);
                w_no_purify *= w_edges[k];
            }

            int slots_no = 1 + (int)ceil(log2(max(h, 1)));
            double w_no_decayed = w_no_purify * pow(w_decay_per_slot, slots_no);

            diag[h].total++;
            diag[h].sum_w += w_no_decayed;

            // 嘗試 1, 2, 3 rounds purification，找最少 rounds 就能過 threshold 的
            int best_rounds = -1;
            double best_w_pur = 0;
            for (int rr = 1; rr <= max_purify_rounds; rr++) {
                double w_pur = 1.0;
                for (int k = 0; k < h; k++)
                    w_pur *= purified_werner(w_edges[k], rr);
                int slots_pur = (1 + rr) + (int)ceil(log2(max(h, 1)));
                double w_pur_decayed = w_pur * pow(w_decay_per_slot, slots_pur);
                if (w_pur_decayed >= w_th) {
                    best_rounds = rr;
                    best_w_pur = w_pur_decayed;
                    break;  // 用最少 rounds 就行
                }
            }

            if (w_no_decayed < w_th && best_rounds > 0) {
                // A: 嚴格甜蜜點 — 不做 purify 過不了，做 best_rounds 輪能過
                diag[h].sweet++;
                // score: rounds 少的優先（資源省），同 rounds 內 margin 大的優先
                double score = 3.0 - best_rounds * 0.3 + best_w_pur / w_th * 0.1;
                candidates.push_back({score, {i, j}});
                candidates.push_back({score, {j, i}});
            } else if (w_no_decayed >= w_th && w_no_decayed < w_th * margin_ratio
                       && best_rounds > 0) {
                // B: 邊緣受益者
                diag[h].marginal++;
                double boost = best_w_pur / w_no_decayed;
                double score = 1.0 + boost / 10.0;
                candidates.push_back({score, {i, j}});
                candidates.push_back({score, {j, i}});
            } else if (w_no_decayed >= w_th * margin_ratio) {
                diag[h].pass_no++;
            } else {
                diag[h].fail_both++;
            }
        }
    }

    // 診斷
    cerr << "\033[1;33m" << "[purify_needed] diagnostics (w_th=" << w_th
         << ", margin=" << margin_ratio << ", max_rounds=" << max_purify_rounds << "):" << "\033[0m" << endl;
    for (auto &[hop, stats] : diag) {
        cerr << "  hop=" << hop
             << " | pairs=" << stats.total
             << " | comfy_pass=" << stats.pass_no
             << " | marginal=" << stats.marginal
             << " | strict_sweet=" << stats.sweet
             << " | fail_both=" << stats.fail_both
             << " | avg_w_no=" << (stats.total > 0 ? stats.sum_w / stats.total : 0)
             << endl;
    }
    cerr << "\033[1;33m" << "[purify_needed] total candidates=" << candidates.size()
         << " (strict + marginal)" << "\033[0m" << endl;

    if (candidates.empty()) {
        cerr << "\033[1;31m" << "[purify_needed] WARNING: no pairs found! "
             << "Consider lowering min_fidelity or raising fidelity_threshold or min_hop."
             << "\033[0m" << endl;
        return {};
    }

    // 按 hop 數分桶，每桶內 shuffle，然後 round-robin 均勻抽取
    map<int, vector<SDpair>> hop_buckets;
    for (auto &[score, sd] : candidates) {
        vector<int> p = bfs_path(sd.first, sd.second);
        int h = p.empty() ? 0 : (int)p.size() - 1;
        hop_buckets[h].push_back(sd);
    }

    random_device rd;
    default_random_engine gen(rd());

    // 每桶 shuffle
    vector<pair<int, vector<SDpair>>> buckets_vec;
    for (auto &[h, vec] : hop_buckets) {
        shuffle(vec.begin(), vec.end(), gen);
        buckets_vec.push_back({h, vec});
    }

    cerr << "\033[1;33m" << "[purify_needed] hop distribution:";
    for (auto &[h, vec] : buckets_vec)
        cerr << " hop" << h << "=" << vec.size();
    cerr << "\033[0m" << endl;

    // Round-robin：輪流從每個 hop 桶取，每次取 2~4 個同 SD pair
    vector<SDpair> requests;
    vector<int> pos(buckets_vec.size(), 0);
    uniform_int_distribution<int> rep_dist(2, 4);
    int bucket_idx = 0;
    while ((int)requests.size() < requests_cnt) {
        // 找到一個還有 pair 的桶
        bool found = false;
        for (int try_cnt = 0; try_cnt < (int)buckets_vec.size(); try_cnt++) {
            int bi = (bucket_idx + try_cnt) % (int)buckets_vec.size();
            if (pos[bi] < (int)buckets_vec[bi].second.size()) {
                int rep = min(rep_dist(gen), requests_cnt - (int)requests.size());
                for (int r = 0; r < rep; r++)
                    requests.push_back(buckets_vec[bi].second[pos[bi]]);
                pos[bi]++;
                bucket_idx = (bi + 1) % (int)buckets_vec.size();
                found = true;
                break;
            }
        }
        if (!found) {
            // 所有桶都用完，循環重來
            for (int i = 0; i < (int)pos.size(); i++) pos[i] = 0;
            for (auto &[h, vec] : buckets_vec) shuffle(vec.begin(), vec.end(), gen);
        }
    }
    requests.resize(requests_cnt);

    shuffle(requests.begin(), requests.end(), gen);
    return requests;
}

int main(){
    string file_path = "../data/";

    map<string, double> default_setting;
    default_setting["num_nodes"] = 100;
    default_setting["request_cnt"] = 50;
    default_setting["entangle_lambda"] = 0.045;
    default_setting["time_limit"] = 20;
    default_setting["avg_memory"] = 20; // 16
    default_setting["tao"] = 0.002;
    default_setting["path_length"] = 4;
    // === Purification 甜蜜點參數 ===
    // min_fidelity 降低使 w_e 落在 ~0.73-0.87，讓 3-4 hop 就需要 purify
    // fidelity_threshold 提高使 w_th 更嚴格，擴大甜蜜點範圍
    // 原始值: min_fidelity=0.89, fidelity_threshold=0.7
    default_setting["min_fidelity"] = 0.78;
    default_setting["max_fidelity"] = 0.93;
    default_setting["swap_prob"] = 0.9;
    default_setting["fidelity_threshold"] = 0.75;
    default_setting["entangle_time"] = 0.00025;
    default_setting["entangle_prob"] = 0.01;
    default_setting["Zmin"]=0.02702867239;
    default_setting["bucket_eps"]=0.01;
    default_setting["time_eta"]=0.001;
    default_setting["hop_count"]=3;
    default_setting["delta_P"]=0.000000001;
    map<string, vector<double>> change_parameter;
    change_parameter["request_cnt"] = {10,20,30,40,50,60,70,80};
    change_parameter["num_nodes"] = {40, 70, 100, 130, 160};
    change_parameter["min_fidelity"] = {0.6, 0.7, 0.8, 0.9, 0.95};
    change_parameter["avg_memory"] = {2,4,6, 8, 10,12,14};
    // change_parameter["tao"] = {0.3, 0.4, 0.5, 0.6, 0.7};
    change_parameter["tao"] = {0.0015, 0.00175, 0.002,0.00225,0.0025};
    change_parameter["path_length"] = {3, 6, 9, 12, 15};
    change_parameter["swap_prob"] = {0.6, 0.7, 0.8, 0.9,0.95};
    change_parameter["fidelity_threshold"] = {0.55,0.6,0.65,0.7, 0.75, 0.8,0.85,0.9};
    change_parameter["time_limit"] = {3,5,7, 9, 11, 13, 15};
    change_parameter["entangle_lambda"] = {0.0125, 0.025, 0.035, 0.045, 0.055, 0.065};
    change_parameter["entangle_time"] = {0.0001, 0.00025, 0.0004, 0.00055, 0.0007,0.00085,0.001};
    change_parameter["entangle_prob"] = {0.0001, 0.001, 0.01, 0.1, 1};
    change_parameter["hop_count"] = {1,2,3,4,5,6};
    //change_parameter["Zmin"]={0.028,0.150,0.272,0.394,0.518};
    change_parameter["bucket_eps"]={0.00001,0.0001,0.001,0.01,0.1};
    change_parameter["time_eta"]={0.00001,0.0001,0.001,0.01,0.1};
    int round = 1;
    vector<vector<SDpair>> default_requests(round);
    #pragma omp parallel for
    for(int r = 0; r < round; r++) {
        int num_nodes = default_setting["num_nodes"];
        int avg_memory = default_setting["avg_memory"];
        // int request_cnt = default_setting["request_cnt"];
        int time_limit = default_setting["time_limit"];
        double min_fidelity = default_setting["min_fidelity"];
        double max_fidelity = default_setting["max_fidelity"];
        double Zmin=default_setting["Zmin"];
        double bucket_eps=default_setting["bucket_eps"];
        double time_eta=default_setting["time_eta"];
        double swap_prob = default_setting["swap_prob"];
        double fidelity_threshold = default_setting["fidelity_threshold"];
        int length_upper = default_setting["path_length"] + 1;
        int length_lower = default_setting["path_length"] - 1;
        map<string, double> input_parameter = default_setting;
        vector<map<string, map<string, double>>> result(round);
        // double entangle_lambda = input_parameter["entangle_lambda"];
        // double entangle_time = input_parameter["entangle_time"];
        double entangle_prob = input_parameter["entangle_prob"];
        string filename = file_path + "input/round_" + to_string(r) + ".input";
        string command = "python3 graph_generator.py ";
        double A = 0.25, B = 0.75, tao = default_setting["tao"], T = 10, n = 2;
        // derandom
        string parameter = to_string(num_nodes);
        cerr << (command + filename + " " + parameter) << endl;
        if(system((command + filename + " " + parameter).c_str()) != 0){
            cerr<<"error:\tsystem proccess python error"<<endl;
            exit(1);
        }
        Graph graph(filename, time_limit, swap_prob, avg_memory, min_fidelity, max_fidelity, fidelity_threshold, A, B, n, T, tao,Zmin,bucket_eps,time_eta,input_parameter["delta_P"]);
        // === 混合生成 request ===
        // 目標：purify 演算法明顯優於 non-purify，但 non-purify 也能接一些（答案不為 0）
        // 比例：~60% purify 甜蜜點（只有 purify 能接）+ ~40% baseline（所有人都能接）
        //
        // 容量估算：100 nodes × avg_mem=20 × time_limit=20 = 40k memory-timeslots
        // 每個 purify shape ~30-45 mem-slots，baseline ~16-25 → 實際可服務 ~30-60 個
        // 生成量要大於可服務量（讓演算法有選擇空間），但不要太多（浪費 LP 時間）
        int total_cnt = 80;  // 略高於可服務量，讓 LP 有選擇空間
        double purify_ratio = 0.6;
        int purify_cnt = (int)(total_cnt * purify_ratio);    // ~48 purify
        int baseline_cnt = total_cnt - purify_cnt;            // ~32 baseline

        auto purify_reqs = generate_requests_purify_needed(graph, purify_cnt, 2);
        // baseline: fidelity 夠高，不需要 purify 就能通過（所有演算法都能服務）
        auto baseline_reqs = generate_requests_fid(graph, baseline_cnt, fidelity_threshold + 0.01, 2);
        // 如果 baseline 不夠，放寬條件
        if ((int)baseline_reqs.size() < baseline_cnt) {
            cerr << "[requests] baseline only got " << baseline_reqs.size() << "/" << baseline_cnt
                 << ", relaxing fid_th to 0.6" << endl;
            baseline_reqs = generate_requests_fid(graph, baseline_cnt, 0.6, 2);
        }

        if (purify_reqs.empty() && baseline_reqs.empty()) {
            cerr << "[fallback] no candidates at all, using generate_requests_fid with loose params" << endl;
            default_requests[r] = generate_requests_fid(graph, total_cnt, 0.5, 2);
        } else {
            // 交錯排列：purify_per_cycle 個 purify + 1 個 baseline，確保任何前綴都有兩種
            // 比例: cycle size = purify_per_cycle + 1, purify 佔 purify_per_cycle/(purify_per_cycle+1)
            default_requests[r].clear();
            int pi = 0, bi = 0;
            // purify_ratio=0.6 → 3 purify + 2 baseline per cycle
            int pur_per_cycle = 3, base_per_cycle = 2;
            while ((int)default_requests[r].size() < total_cnt) {
                for (int k = 0; k < pur_per_cycle && (int)default_requests[r].size() < total_cnt; k++) {
                    if (pi < (int)purify_reqs.size())
                        default_requests[r].push_back(purify_reqs[pi++]);
                    else if (bi < (int)baseline_reqs.size())
                        default_requests[r].push_back(baseline_reqs[bi++]);
                }
                for (int k = 0; k < base_per_cycle && (int)default_requests[r].size() < total_cnt; k++) {
                    if (bi < (int)baseline_reqs.size())
                        default_requests[r].push_back(baseline_reqs[bi++]);
                    else if (pi < (int)purify_reqs.size())
                        default_requests[r].push_back(purify_reqs[pi++]);
                }
            }
            default_requests[r].resize(total_cnt);
        }

        // === 印出最終 request 的詳細統計 ===
        {
            map<int, int> hop_dist;
            for (auto &sd : default_requests[r]) {
                int d = graph.distance(sd.first, sd.second);
                hop_dist[d]++;
            }
            cerr << "\033[1;36m"
                 << "========== Request Generation Done ==========" << endl
                 << "  total=" << default_requests[r].size()
                 << " | purify_sweet=" << purify_reqs.size()
                 << " | baseline=" << baseline_reqs.size() << endl
                 << "  hop distribution: ";
            for (auto &[h, cnt] : hop_dist)
                cerr << h << "hop=" << cnt << " ";
            cerr << endl
                 << "  design: ~" << (int)(purify_ratio*100) << "% need purify (ZFA2 advantage)"
                 << ", ~" << (int)((1-purify_ratio)*100) << "% all-algo baseline (ZFA/MyAlgo1 also score)" << endl
                 << "================================================"
                 << "\033[0m" << endl;
        }
        assert(!default_requests[r].empty());
    }




    // vector<string> X_names = {"time_limit", "request_cnt", "num_nodes", "avg_memory", "tao"};
    //vector<string> X_names = {"request_cnt"};
    vector<string> X_names = { "request_cnt", "time_limit", "tao",  "fidelity_threshold" , "avg_memory","hop_count" };
    //vector<string> X_names = {"Zmin","bucket_eps","time_eta"};
    vector<string> Y_names = {"fidelity_gain", "succ_request_cnt"};
    vector<string> algo_names = {"ZFA2","ZFA","MyAlgo1", "MyAlgo2", "MyAlgo3"};
    // init result


    vector<PathMethod*> path_methods;
    path_methods.emplace_back(new Greedy());
    /* path_methods.emplace_back(new QCAST());
    path_methods.emplace_back(new REPS()); */
    for(PathMethod *path_method : path_methods) {

        for(string X_name : X_names) {
            for(string Y_name : Y_names){
                if(path_method->get_name() != "Greedy" && X_name != "request_cnt")
                    continue; 
                string filename = "ans/" + path_method->get_name() + "_" + X_name + "_" + Y_name + ".ans";
                fstream file( file_path + filename, ios::out );
            }
        }

        for(string X_name : X_names) {
            if(path_method->get_name() != "Greedy" && X_name != "request_cnt")
                continue; 
                
            map<string, double> input_parameter = default_setting;

            for(double change_value : change_parameter[X_name]) {
                vector<map<string, map<string, double>>> result(round);
                input_parameter[X_name] = change_value;

                // int num_nodes = input_parameter["num_nodes"];
                int avg_memory = input_parameter["avg_memory"];
                int request_cnt = input_parameter["request_cnt"];
                int time_limit = input_parameter["time_limit"];
                double min_fidelity = input_parameter["min_fidelity"];
                double max_fidelity = input_parameter["max_fidelity"];
                double Zmin = input_parameter["Zmin"];
                double bucket_eps=input_parameter["bucket_eps"];
                double time_eta=input_parameter["time_eta"];
                // double entangle_lambda = input_parameter["entangle_lambda"];
                // double entangle_time = input_parameter["entangle_time"];
                double entangle_prob = input_parameter["entangle_prob"];
                double swap_prob = input_parameter["swap_prob"];
                double fidelity_threshold = input_parameter["fidelity_threshold"];
                int hop_count = input_parameter["hop_count"];
                // int length_upper, length_lower;
                // if(input_parameter["path_length"] == -1) {
                //     length_upper = num_nodes;
                //     length_lower = 6;
                // } else {
                //     length_upper = input_parameter["path_length"] + 1;
                //     length_lower = input_parameter["path_length"] - 1;
                // }

                int sum_has_path = 0;
                //#pragma omp parallel for
                for(int r = 0; r < round; r++) {
                    string filename = file_path + "input/round_" + to_string(r) + ".input";
                    ofstream ofs;
                    ofs.open(file_path + "log/" + path_method->get_name() + "_" + X_name + "_in_" + to_string(change_value) + "_Round_" + to_string(r) + ".log");

                    time_t now = time(0);
                    char* dt = ctime(&now);
                    cerr  << "時間 " << dt << endl << endl;
                    ofs << "時間 " << dt << endl << endl;




                    double A = 0.25, B = 0.75, tao = input_parameter["tao"], T = 10, n = 2;
                    Graph graph(filename, time_limit, swap_prob, avg_memory, min_fidelity, max_fidelity, fidelity_threshold, A, B, n, T, tao,Zmin,bucket_eps,time_eta,input_parameter["delta_P"]);

                    ofs << "--------------- in round " << r << " -------------" <<endl;
                    vector<pair<int, int>> requests;
                    if(hop_count==3){
                        int idx=0;
                        for(int i = 0; i < request_cnt; i++) {
                            /* while(graph.get_ini_fid(default_requests[r][idx].first,default_requests[r][idx].second)<fidelity_threshold){
                                idx=(idx+1)%default_requests[r].size();
                            } */
                            requests.emplace_back(default_requests[r][idx]);
                            idx=(idx+1)%default_requests[r].size();
                        }
                    }
                    else{
                        requests=generate_requests_fid(graph,request_cnt,0,hop_count);
                    }
                    Graph path_graph = graph;
                    path_graph.increase_resources(10);
                    PathMethod *new_path_method;
                    if(path_method->get_name() == "Greedy") new_path_method = new Greedy();
                    else if(path_method->get_name() == "QCAST") new_path_method = new QCAST();
                    else if(path_method->get_name() == "REPS") new_path_method = new REPS();
                    else {
                        cerr << "unknown path method" << endl;
                        assert(false);
                    }

                    new_path_method->build_paths(path_graph, requests);
                    cout << "found path" << endl;
                    map<SDpair, vector<Path>> paths = new_path_method->get_paths();
                    map<SDpair, set<Path>> paths_st;
                    for(auto [sdpair, pathss] : paths) {
                        for(Path path : pathss) {
                            paths_st[sdpair].insert(path);
                        }
                    }

                    paths.clear();
                    for(auto [sdpair, pathss] : paths_st) {
                        for(Path path : pathss) {
                            paths[sdpair].push_back(path);
                        }
                    }

                    int path_len = 0, path_cnt = 0, mx_path_len = 0;

                    int has_path = 0;
                    for(SDpair sdpair : requests) {
                        int mi_path_len = INF;
                        has_path += !paths[sdpair].empty();
                        for(Path path : paths[sdpair]) {
                            mi_path_len = min(mi_path_len, (int)path.size());
                            for(int i = 1; i < (int)path.size(); i++) {
                                assert(graph.adj_set[path[i]].count(path[i - 1]));
                            }
                        }
                        if(mi_path_len != INF) {
                            mx_path_len = max(mx_path_len, mi_path_len);
                            path_cnt++;
                            path_len += mi_path_len;
                        }
                    }

                    sum_has_path += has_path;
                    cerr << "Path method: " << path_method->get_name() << "\n";
                    cerr << "Request cnt: " << request_cnt << "\n";
                    cerr << "Has Path cnt: " << has_path << "\n";
                    cerr << "Avg path length = " << path_len / (double)path_cnt << "\n";
                    cerr << "Max path length = " << mx_path_len << "\n";
                    vector<AlgorithmBase*> algorithms;
                    //algorithms.emplace_back(new WernerAlgo_UB(graph,requests,paths));
                    algorithms.emplace_back(new WernerAlgo3(graph,requests,paths));  // ZFA_UB (LP upper bound with purify)
                    algorithms.emplace_back(new WernerAlgo2(graph,requests,paths));
                    algorithms.emplace_back(new WernerAlgo(graph,requests,paths));
                    if(X_name!="Zmin"&&X_name!="bucket_eps"&&X_name!="time_eta"){
                        algorithms.emplace_back(new MyAlgo1(graph, requests, paths));
                        algorithms.emplace_back(new MyAlgo2(graph, requests, paths));
                        algorithms.emplace_back(new MyAlgo3(graph, requests, paths));
                        /* algorithms.emplace_back(new MyAlgo4(graph, requests, paths));
                        algorithms.emplace_back(new MyAlgo5(graph, requests, paths));
                        algorithms.emplace_back(new MyAlgo6(graph, requests, paths)); */
                    }


                    //#pragma omp parallel for
                    for(int i = 0; i < (int)algorithms.size(); i++) {
                        algorithms[i]->run();
                    }



                    for(int i = 0; i < (int)algorithms.size(); i++) {
                        for(string Y_name : Y_names) {
                            result[r][algorithms[i]->get_name()][Y_name] = algorithms[i]->get_res(Y_name);
                        }
                    }

                    now = time(0);
                    dt = ctime(&now);
                    cerr << "時間 " << dt << endl << endl;
                    ofs << "時間 " << dt << endl << endl;
                    ofs.close();

                    for(auto &algo : algorithms){
                        delete algo;
                    }
                    algorithms.clear();

                }

                map<string, map<string, double>> sum_res;
                // for(string algo_name : algo_names){
                //     for(int r = 0; r < round; r++){
                //         result[r][algo_name]["waiting_time"] /= result[T][algo_name]["total_request"];
                //         result[r][algo_name]["encode_ratio"] = result[T][algo_name]["encode_cnt"] / (result[T][algo_name]["encode_cnt"] + result[T][algo_name]["unencode_cnt"]);
                //         result[r][algo_name]["succ-finished_ratio"] = result[T][algo_name]["throughputs"] / result[T][algo_name]["finished_throughputs"];
                //         result[r][algo_name]["fail-finished_ratio"] = 1 - result[T][algo_name]["succ-finished_ratio"];
                //         result[r][algo_name]["path_length"] = result[T][algo_name]["path_length"] / result[T][algo_name]["finished_throughputs"];
                //         result[r][algo_name]["divide_cnt"] = result[T][algo_name]["divide_cnt"] / result[T][algo_name]["finished_throughputs"];
                //         result[r][algo_name]["use_memory_ratio"] = result[T][algo_name]["use_memory"] / result[T][algo_name]["total_memory"];
                //         result[r][algo_name]["use_channel_ratio"] = result[T][algo_name]["use_channel"] / result[T][algo_name]["total_channel"];
                //     }
                // }

                for(string Y_name : Y_names) {
                    string filename = "ans/" + path_method->get_name() + "_" + X_name + "_" + Y_name + ".ans";
                    ofstream ofs;
                    ofs.open(file_path + filename, ios::app);
                    ofs << change_value << ' ';

                    for(string algo_name : algo_names){
                        for(int r = 0; r < round; r++){
                            sum_res[algo_name][Y_name] += result[r][algo_name][Y_name];
                        }
                        ofs << sum_res[algo_name][Y_name] / round << ' ';
                    }
                    ofs << endl;
                    ofs.close();
                }
            }
        }
    }
    return 0;
}