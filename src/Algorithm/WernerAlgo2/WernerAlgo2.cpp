#include "WernerAlgo2.h"
#include <fstream> 
#include <iostream>
#include <cmath>

using namespace std;

WernerAlgo2::WernerAlgo2(Graph graph,vector<pair<int,int>> requests,map<SDpair, vector<Path>> paths): AlgorithmBase(graph, requests, paths)
{
    algorithm_name = "ZFA2";
}

void WernerAlgo2::variable_initialize() {
    // 與 MyAlgo1 類似：初始化 dual 與目標
    int m = (int)requests.size()
          + graph.get_num_nodes() * graph.get_time_limit();

    double delta = (1 + epsilon) * (1.0 / pow((1 + epsilon) * m, 1.0 / epsilon));
    obj = m * delta;

    alpha.assign(requests.size(), delta);
    x.clear();
    x.resize(requests.size());
    int V = graph.get_num_nodes();
    int T = graph.get_time_limit();
    dpp.eps_bucket = graph.get_bucket_eps();
    double F_th=graph.get_fidelity_threshold();
    double w_th=(4.0*F_th-1.0)/3.0;
    dpp.Zhat = sqrt(-log(w_th))+1e-9;
    dpp.Zmin = graph.get_Zmin();
    dpp.T    = time_limit-1;
    dpp.tau_max=min(time_limit-1,5);
    dpp.eta  = graph.get_tao()/graph.get_time_limit();
    dpp.deltaP=graph.get_delta_P();
    beta.assign(V, vector<double>(T, INF));

    for (int v = 0; v < V; ++v) {
        for (int t = 0; t < T; ++t) {
            int cap = graph.get_node_memory_at(v, t);
            beta[v][t] = (cap == 0) ? INF : (delta / cap);
        }
    }
}

Shape_vector WernerAlgo2::separation_oracle(){
    // 針對每個 request 找最小成本 shape，選最好的回傳
    double most_violate=1e9;
    bool flag=0;
    Shape_vector todo_shape;
    vector<int> best_purify_rounds;
    cerr << "[ZFA2:oracle] requests.size()=" << requests.size() << endl;
    for(int i=0;i<requests.size();i++){
        int src=requests[i].first,dst=requests[i].second;
        vector<Path> paths=get_paths(src,dst);
        cerr << "[ZFA2:oracle] req " << i << " (" << src << "->" << dst << ") paths=" << paths.size() << endl;
        for(int p=0;p<paths.size();p++){
            // Example initialization: adjust T and n as needed for your context
            int T=dpp.T+5;
            int n=paths[p].size()+5;
            cerr << "[ZFA2:oracle] path " << p << " len=" << paths[p].size() << " DP_table T=" << T << " n=" << n << endl;
            DP_table.clear();
            DP_table.resize(T);
            for(int i=0;i<DP_table.size();i++){
                DP_table[i].resize(n);
                for(int j=0;j<DP_table[i].size();j++)
                    DP_table[i][j].resize(n);
            }
            for(int t=1;t<=dpp.T;t++){
                cerr << "[ZFA2:oracle] run_dp_in_t t=" << t << "/" << (int)dpp.T << endl;
                run_dp_in_t(paths[p],dpp,t);
                cerr << "[ZFA2:oracle] eval_best_J t=" << t << endl;
                auto cur_val=eval_best_J(0,paths[p].size()-1,t,alpha[i]);
                if(cur_val.first<most_violate){
                    most_violate=cur_val.first;
                    vector<int> cur_rounds;
                    cerr << "[ZFA2:oracle] backtrack_shape t=" << t << endl;
                    todo_shape=backtrack_shape(cur_val.second,paths[p],cur_rounds);
                    best_purify_rounds=cur_rounds;
                }
            }
        }
    }
    if(!todo_shape.empty()){
        shape_purify_map[todo_shape]=best_purify_rounds;
    }
    return todo_shape;
}
WernerAlgo2::ZLabel WernerAlgo2::gen_leaf_label(int s,int e,int st,int tlen) {
    double Bleaf=0.0;
    if(st-tlen<0) return ZLabel();
    // bounds check
    if(s >= (int)beta.size() || e >= (int)beta.size()) {
        cerr << "[ZFA2:leaf] ERROR s=" << s << " or e=" << e << " >= beta.size()=" << beta.size() << endl;
        return ZLabel();
    }
    for(int i=0;i<=tlen;i++){
        if(st-i < 0 || st-i >= (int)beta[s].size()) {
            cerr << "[ZFA2:leaf] ERROR st-i=" << st-i << " out of beta[s].size()=" << beta[s].size() << endl;
            return ZLabel();
        }
        double bt=beta[s][st-i]+beta[e][st-i];
        Bleaf+=bt*Purify_in_vt[tlen-1][i];
    }
    double w_ini=graph.get_link_werner(s,e);
    double w_cur=w_ini,p_cur=graph.get_entangle_succ_prob(s,e);
    for(int i=1;i<=tlen-1;i++){
        p_cur*=(9.0L*w_cur*w_ini-3.0L*w_ini-3.0L*w_cur+5)/8.0L;
        w_cur=(3*w_cur*w_ini+3*w_ini+3*w_cur-1.0L)/(9*w_cur*w_ini-3*w_ini-3*w_cur+5.0L);
        p_cur*=graph.get_entangle_succ_prob(s,e);
    }
    double Zleaf=sqrt(-log(w_cur));
    double Pleaf=log(p_cur);
    if(Zleaf>dpp.Zhat) return ZLabel();
    return ZLabel(Bleaf,Zleaf,Pleaf,Op::LEAF,tlen-1,s,e,st,-1);
} 
void WernerAlgo2::run_dp_in_t(const Path& path, const DPParam& dpp,int t) {
    const int T = graph.get_time_limit();
    const int n = (int)path.size();
    cerr << "[ZFA2:dp] run_dp_in_t t=" << t << " T=" << T << " n=" << n
         << " DP_table.size()=" << DP_table.size() << endl;
    if(t >= (int)DP_table.size()) { cerr << "[ZFA2:dp] ERROR t >= DP_table.size()" << endl; return; }

    // -------- t = 1..T-1 外圈時間迴圈 --------
    for(int a=0;a<n-1;a++)
        for(int b=a+1;b<n;b++){
            if(a >= (int)DP_table[t].size() || b >= (int)DP_table[t][a].size()) {
                cerr << "[ZFA2:dp] ERROR a=" << a << " b=" << b << " out of bounds" << endl;
                continue;
            }
            int s=path[a],e=path[b];
            vector<ZLabel> cand;
            //leaf
            if(a+1==b){
                for(int i=0;i<=purify_time;i++){
                    if(t-i-1<=0) continue;
                    ZLabel L=gen_leaf_label(s,e,t,i+1);
                    if(L.Z<=dpp.Zhat){
                        L.ent_time={t-i-1,t};
                        cand.push_back(L);
                    }
                }
            }
            //continue
            auto pre=DP_table[t-1][a][b];
            for(int p_id=0;p_id<pre.size();p_id++){
                double Zp=pre[p_id].Z+dpp.eta;
                if(Zp<=dpp.Zhat){
                    double Bp=pre[p_id].B+beta[s][t]+beta[e][t];
                    double Pp=pre[p_id].P;
                    ZLabel L(Bp,Zp,Pp,Op::CONT,-1,a,b,t,-1,p_id);
                    cand.push_back(L);
                }
            }
            //merge
            for(int k=a+1;k<b;k++){
                auto L1=DP_table[t-1][a][k],L2=DP_table[t-1][k][b];
                if(L1.size()==0||L2.size()==0) continue;
                for(int lid=0;lid<L1.size();lid++)
                    for(int rid=0;rid<L2.size();rid++){
                        auto left_seg=L1[lid],right_seg=L2[rid];
                        double Zp=sqrt((left_seg.Z+dpp.eta)*(left_seg.Z+dpp.eta)+
                                        (right_seg.Z+dpp.eta)*(right_seg.Z+dpp.eta));
                        double swap_prob=log(graph.get_node_swap_prob(path[k]));
                        double Pp=left_seg.P+right_seg.P+swap_prob;
                        if(Zp<=dpp.Zhat){
                            double Bp=left_seg.B+right_seg.B+beta[s][t]+beta[e][t];
                            ZLabel L(Bp,Zp,Pp,Op::MERGE,-1,a,b,t,k,-1,lid,rid);
                            cand.push_back(L);
                        }
                    }
            }
            vector<ZLabel> non_leaf;
            for (auto& L : cand) {
                if (L.op != Op::LEAF)
                    non_leaf.push_back(L);
            }
            bucket_by_ZP(non_leaf);// trimming non-leaf
            cand.erase(
                remove_if(cand.begin(), cand.end(),
                        [](const ZLabel& L){ return L.op != Op::LEAF; }),
                cand.end());
            cand.insert(cand.end(), non_leaf.begin(), non_leaf.end());
            DP_table[t][a][b] = cand;
        }
}

void WernerAlgo2::pareto_prune_byZ(vector<ZLabel>& cand) {
    if (cand.empty()) return;
    sort(cand.begin(), cand.end(), [](const ZLabel& x, const ZLabel& y){
        if(x.Z!=y.Z) return x.Z < y.Z;
        return x.B<y.B;
    });
    vector<ZLabel> kept;
    double bestB = INF;
    for (auto& L : cand) {
        if (L.B + 1e-12 < bestB) {
            kept.push_back(L);
            bestB = L.B;
        }
    }
    cand.swap(kept);
}

void WernerAlgo2::bucket_by_ZP(vector<ZLabel>& cand) {
    if (cand.empty()) return;
    double q=1+dpp.eps_bucket;
    double invLogQ=1.0/log(q);
    double deltaP=dpp.deltaP;
    map<pair<double,double>,ZLabel> buckets;
    for(auto L:cand){
        double kW;
        if(L.Z<=dpp.Zmin) kW=0.0;
        else{
            kW=floor(log(L.Z/dpp.Zmin)*invLogQ+1e-12);
            if(kW<0) kW=0.0;
        } 
        double kP=floor(-L.P/deltaP);
        auto key=make_pair(kW,kP);
        if(buckets.count(key)==0||L.B<buckets[key].B)
            buckets[key]=L;
    }
    vector<ZLabel> bucketed;
    for(auto L:buckets)
        bucketed.push_back(L.second);
    //pareto_prune_byZ(bucketed);
    sort(bucketed.begin(), bucketed.end(), [](const ZLabel& x, const ZLabel& y){
        return x.Z < y.Z;
    });
    cand.swap(bucketed);
}

Shape_vector WernerAlgo2::backtrack_shape(ZLabel leaf,const vector<int>& path, vector<int>& out_purify_rounds){
    int left_id=path[leaf.a],right_id=path[leaf.b];
    if(leaf.op==Op::LEAF){
        Shape_vector result;
        if (leaf.ent_time.size() < 2) return Shape_vector{};
        assert(leaf.ent_time.size() == 2);
        for(int i=0;i<leaf.purify_type+1;i++){
            int start_time=leaf.ent_time[0]+i,end_time=leaf.ent_time[0]+i;
            for(int j=0;j<Purify_in_vt[leaf.purify_type][i];j++){
                result.push_back({left_id,{{start_time,end_time}}});
                result.push_back({right_id,{{start_time,end_time}}});
            }
        }
        // LEAF: 一條 link，purify rounds = purify_type (即 tlen-1)
        out_purify_rounds = { leaf.purify_type };
        return result;
    }
    if(leaf.op==Op::CONT){
        assert(leaf.parent_id>=0&&leaf.parent_id<DP_table[leaf.t-1][leaf.a][leaf.b].size());
        ZLabel pre_label=DP_table[leaf.t-1][leaf.a][leaf.b][leaf.parent_id];
        // CONT: idle 不改變 purify rounds，直接透傳
        Shape_vector last_time=backtrack_shape(pre_label,path,out_purify_rounds);
        if (last_time.empty()) return Shape_vector{};
        auto & prel=last_time.front().second[0],&prer=last_time.back().second[0];
        assert(last_time.front().first==path[leaf.a]);
        assert(last_time.back().first==path[leaf.b]);
        assert(prel.second==leaf.t-1);
        assert(prer.second==leaf.t-1);
        prel.second++;
        prer.second++;
        return last_time;
    }
    if(leaf.op==Op::MERGE){
        Shape_vector left_result,right_result,result;
        assert(leaf.k>=0);
        int k_id=path[leaf.k];
        // MERGE: 左右各自遞迴取 purify rounds，再串接
        vector<int> left_rounds, right_rounds;
        ZLabel left_leaf=DP_table[leaf.t-1][leaf.a][leaf.k][leaf.left_id];
        left_result=backtrack_shape(left_leaf,path,left_rounds);
        ZLabel right_leaf=DP_table[leaf.t-1][leaf.k][leaf.b][leaf.right_id];
        right_result=backtrack_shape(right_leaf,path,right_rounds);
        if(DEBUG) {
            assert(left_result.front().first == path[leaf.a]);
            assert(left_result.front().second[0].second == leaf.t - 1);
            assert(left_result.front().second.size() == 1);
            assert(left_result.back().first == k_id);
            assert(right_result.front().first == k_id);
            assert(right_result.back().first == path[leaf.b]);
            assert(right_result.back().second[0].second == leaf.t - 1);
            assert(left_result.back().second.size() == 1);
        }

        for(int i = 0; i < (int)left_result.size(); i++) {
            result.push_back(left_result[i]);
        }
        result.back().second.push_back(right_result.front().second.front());
        for(int i = 1; i < (int)right_result.size(); i++) {
            result.push_back(right_result[i]);
        }

        result.front().second[0].second++;
        result.back().second[0].second++;
        // 串接左右的 purify rounds
        out_purify_rounds = left_rounds;
        out_purify_rounds.insert(out_purify_rounds.end(), right_rounds.begin(), right_rounds.end());
        return result;
    }
    out_purify_rounds.clear();
    return Shape_vector{};
}
int WernerAlgo2::split_dis(int s,int d,WernerAlgo2::ZLabel& L){
    if(L.op!=WernerAlgo2::Op::MERGE||L.k<0) return 1e18/4;
    int mid=(s+d)/2;
    return abs(mid-L.k);
}
pair<double,WernerAlgo2::ZLabel> WernerAlgo2::eval_best_J(int s, int d, int t, double alp){
    double bestJ=1e18;
    int bestdis=1e18/4;
    int flag=0;
    ZLabel tmp={};
    for(auto L:DP_table[t][s][d]){
        double J=(alp+L.B)*exp(L.Z*L.Z-L.P);
        int dis=split_dis(s,d,L);
        if(J+EPS<bestJ||(fabs(J-bestJ)<=EPS&&dis<bestdis)){
            bestJ=J;
            tmp=L;
            bestdis=dis;
            flag=1;
        }
    }
    if(flag) return {bestJ,tmp};
    else return {INF,tmp};
}

void WernerAlgo2::run() {
    cerr << "[ZFA2] run() start, requests.size()=" << requests.size() << endl;
    int round = 1;
    while (round-- && !requests.empty()) {
        variable_initialize();
        cerr << "[ZFA2] variable_initialize done, dpp.T=" << dpp.T << " dpp.Zhat=" << dpp.Zhat << endl;
        //cerr << "\033[1;31m"<< "[WernerAlgo's parameter] : "<< dpp.Zmin<<" "<<dpp.eps_bucket<<" "<<dpp.eta<< "\033[0m"<< endl;
        int it=0;
        double eps=1e-4;
        while (obj+eps < 1.0) {
            cerr << "[ZFA2] oracle iteration " << it++ << ", obj=" << (double)obj << endl;
            //if(++it>200) break;
            Shape_vector shape=separation_oracle();
            cerr << "[ZFA2] oracle done, shape.size()=" << shape.size() << endl;
            if (shape.empty()) break;
            // 先用MyAlgo1的框架刻出來
            cerr << "[ZFA2:run] shape details:" << endl;
            for(int i=0;i<(int)shape.size();i++){
                cerr << "  shape[" << i << "] node=" << shape[i].first << " intervals=" << shape[i].second.size() << ":";
                for(auto& p : shape[i].second) cerr << " [" << p.first << "," << p.second << "]";
                cerr << endl;
            }
            double q = 1.0;
            for(int i=0;i<shape.size();i++){
                map<int,int> need_amount;
                for(pair<int,int> usedtime:shape[i].second){
                    int start=usedtime.first,end=usedtime.second;
                    for(int t=start;t<=end;t++)
                        need_amount[t]++;
                }
                for(pair<int,int>P:need_amount){
                    int t=P.first;
                    double theta=P.second;
                    cerr << "[ZFA2:run] get_node_memory_at node=" << shape[i].first << " t=" << t << " V=" << graph.get_num_nodes() << " T=" << graph.get_time_limit() << endl;
                    q=min(q,graph.get_node_memory_at(shape[i].first,t)/theta);
                }
            }
            if(q<=1e-10) break;
            int req_idx=-1;
            for(int i=0;i<requests.size();i++){
                int ln=shape.front().first,rn=shape.back().first;
                if(requests[i]==make_pair(ln,rn)){
                    if(req_idx==-1||alpha[req_idx]>alpha[i]){
                        req_idx=i;
                    }
                }
            }
            if(req_idx==-1) break;
            x[req_idx][shape]+=q;
            double ori=alpha[req_idx];
            alpha[req_idx]=alpha[req_idx]*(1+epsilon*q);
            obj+=(alpha[req_idx]-ori);
            for(int i=0;i<shape.size();i++){
                map<int,int> need_amount;
                for(pair<int,int> usedtime:shape[i].second){
                    int start=usedtime.first,end=usedtime.second;
                    for(int t=start;t<=end;t++)
                        need_amount[t]++;
                }

                for(pair<int, int> P : need_amount) {
                    int t = P.first;
                    int node_id = shape[i].first;
                    double theta = P.second;
                    double original = beta[node_id][t];
                    if(graph.get_node_memory_at(node_id, t) == 0) {
                        beta[node_id][t] = INF;
                    } else {
                        beta[node_id][t] = beta[node_id][t] * (1 + epsilon * (q / (graph.get_node_memory_at(node_id, t) / theta)));
                    }
                    obj += (beta[node_id][t] - original) * theta;
                }
            }
        }
        vector<pair<double, Shape_vector>> shapes;

        for(int i = 0; i < (int)requests.size(); i++) {
            for(auto P : x[i]) {
                shapes.push_back({P.second, P.first});
            }
        }

        sort(shapes.begin(), shapes.end(), [](pair<double, Shape_vector> left, pair<double, Shape_vector> right) {
            return left.first > right.first;
        });

        vector<bool> used(requests.size(), false);
        vector<int> finished;
        
        // [新增] 統計資料結構：Key=Duration, Value=List of {W_raw, W_new}
        map<int, vector<pair<double, double>>> purification_stats;

        for(pair<double, Shape_vector> P : shapes) {
            // 用 purify_rounds 構造 Shape（若有的話）
            vector<int> pr;
            if(shape_purify_map.count(P.second))
                pr = shape_purify_map[P.second];
            Shape shape = pr.empty() ? Shape(P.second) : Shape(P.second, pr);
            int request_index = -1;
            for(int i = 0; i < (int)requests.size(); i++) {
                if(used[i] == false && requests[i] == make_pair(shape.get_node_mem_range().front().first, shape.get_node_mem_range().back().first)) {
                    request_index = i;
                }
            }

            if(request_index == -1 || used[request_index]) continue;
            
            // 呼叫 check_resource 時開啟純化 (Enable Purification)
            if(graph.check_resource(shape, true, true)) {
                used[request_index] = true;
                graph.reserve_shape(shape, true);
                finished.push_back(request_index);

                // [新增] 收集純化資訊
                Shape_vector sv = shape.get_node_mem_range();
                vector<int> pur_rounds = shape.get_link_purify_rounds();
                // 遍歷每一段 Link
                for(size_t i = 0; i < sv.size() - 1; ++i) {
                    int rounds = (i < pur_rounds.size()) ? pur_rounds[i] : 0;
                    if (rounds > 0) {
                        int u = sv[i].first;
                        int v = sv[i+1].first;
                        double w_e = graph.get_link_werner(u, v);
                        double w_cur = w_e;
                        for (int r = 0; r < rounds; r++) {
                            double num = 3.0L * w_cur * w_e + 3.0L * w_cur + 3.0L * w_e - 1.0L;
                            double den = 9.0L * w_cur * w_e - 3.0L * w_cur - 3.0L * w_e + 5.0L;
                            w_cur = num / den;
                        }
                        purification_stats[rounds].push_back({w_e, w_cur});
                    }
                }
            }
        }

        sort(finished.rbegin(), finished.rend());
        for(auto fin : finished) {
            requests.erase(requests.begin() + fin);
        }

        // [新增] 將統計結果寫入檔案
        // 檔案路徑與 main.cpp 中的 log 路徑一致 (假設在 ../data/log/)
        // 使用 append 模式，避免覆蓋之前的回合記錄
        string log_file_path = "../data/log/ZFA2_Purification_Stats.txt";
        ofstream log_file(log_file_path, ios::app);
        
        if (log_file.is_open()) {
            if (!purification_stats.empty()) {
                log_file << "--- Round Log ---" << endl;
                for(auto const& [len, vec] : purification_stats) {
                    log_file << "Rounds " << len << " (Count: " << vec.size() << "):" << endl;
                    for(auto const& pair : vec) {
                        double w_raw = pair.first;
                        double w_new = pair.second;
                        // 轉換為 Fidelity 方便閱讀: F = (3W+1)/4
                        double f_raw = (3.0 * w_raw + 1.0) / 4.0;
                        double f_new = (3.0 * w_new + 1.0) / 4.0;
                        
                        log_file << "  W: " << w_raw << " -> " << w_new 
                                 << " | F: " << f_raw << " -> " << f_new << endl;
                    }
                }
                log_file << "-----------------" << endl;
            }
            log_file.close();
        } else {
            cerr << "[Warning] Unable to open log file: " << log_file_path << endl;
        }

        // 同時印到 cerr 方便即時觀察
        cerr << "\n=== [WernerAlgo2 Purification Stats] ===" << endl;
        if (purification_stats.empty()) {
            cerr << "No links were purified (all rounds = 0)." << endl;
        } else {
            for(auto const& [len, vec] : purification_stats) {
                cerr << "Rounds " << len << " (Count: " << vec.size() << "):" << endl;
                // 為了避免洗版，終端機只印數量，詳細內容看檔案
            }
        }
        cerr << "========================================\n" << endl;
    }
    update_res();
    cerr << "[" << algorithm_name << "] end" << endl;
}