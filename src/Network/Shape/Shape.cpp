#include "Shape.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <cassert>
#include "../../config.h"
using namespace std;

Shape::Shape(Shape_vector _node_mem_range):
    node_mem_range(_node_mem_range) {}

Shape::Shape(Shape_vector _node_mem_range, vector<int> _link_purify_rounds):
    node_mem_range(_node_mem_range), link_purify_rounds(_link_purify_rounds) {}

Shape::Shape() {}

Shape_vector Shape::get_node_mem_range() {
    return node_mem_range;
}

vector<int> Shape::get_link_purify_rounds() {
    return link_purify_rounds;
}

// [修改 1] 函式定義必須包含 bool enable_purification (注意：這裡不寫預設值 = false，預設值只寫在 .h)
double Shape::get_fidelity(double _A, double _B, double _n, double _T, double _tao, 
                           map<pair<int, int> , double> F_init, bool enable_purification) {
    A = _A, B = _B, n = _n, T = _T, tao = _tao;
    
    // [修改 2] 將傳入的參數存入成員變數，這樣 recursion_get_fidelity 才能讀到
    this->purification_enabled = enable_purification;

    // cerr << A << " " << B << " " << n << " " << T << " " << tao << endl;
    check_valid();
    return recursion_get_fidelity(0, (int)node_mem_range.size() - 1, F_init);
}

double Shape::recursion_get_fidelity(int left, int right, map<pair<int, int> , double> &F_init) {
    // Base Case: 葉節點 (Link)
    if(left == right - 1) {
        double raw_f = F_init[{node_mem_range[left].first, node_mem_range[right].first}];

        // 若有提供 link_purify_rounds 且該 link 有 purification，使用正確輪數 (論文 Eq.7-8)
        int rounds = 0;
        if (this->purification_enabled && left < (int)link_purify_rounds.size()) {
            rounds = link_purify_rounds[left];
        }

        if (rounds > 0) {
            // Fidelity -> Werner parameter: w = (4F - 1) / 3
            double w_e = (4.0 * raw_f - 1.0) / 3.0;
            double w_cur = w_e;

            // Pumping purification (Eq.7-8)
            for (int r = 0; r < rounds; r++) {
                double num = 3.0 * w_cur * w_e + 3.0 * w_cur + 3.0 * w_e - 1.0;
                double den = 9.0 * w_cur * w_e - 3.0 * w_cur - 3.0 * w_e + 5.0;
                w_cur = num / den;
            }

            // Werner parameter -> Fidelity: F = (3w + 1) / 4
            double purified_f = (3.0 * w_cur + 1.0) / 4.0;
            return pass_tao(purified_f);
        }

        // 無 purification，直接回傳原始 Fidelity 並衰減
        return pass_tao(raw_f);
    }

    // --- 以下為原本的 Merge 遞迴邏輯 (保持不變) ---
    int latest = left + 1;
    for(int i = left + 1; i < right; i++) {
        if(node_mem_range[i].second[0].second > node_mem_range[latest].second[0].second) {
            latest = i;
        }
    }

    double Fa = recursion_get_fidelity(left, latest, F_init);
    double Fb = recursion_get_fidelity(latest, right, F_init);

    int now_swap_time = node_mem_range[latest].second[0].second;
    int next_swap_time = min(node_mem_range[left].second.back().second, node_mem_range[right].second.front().second);
    double pass_time = tao * (next_swap_time - now_swap_time);
    
    // cerr << "l r late = " << left << " " << right << " " << latest << endl;
    // cerr << "Fa Fb result = " << Fa << " " << Fb << " " << Fswap(t2F(F2t(Fa) + pass_time), t2F(F2t(Fb) + pass_time)) << endl;
    return t2F(F2t(Fswap(pass_tao(Fa), pass_tao(Fb))) + (pass_time - tao));
}

double Shape::bar(double F) {
    return (1.0 - F);
}
double Shape::Fswap(double Fa, double Fb) {
    if(Fa <= A + EPS || Fb <= A + EPS) return 0;
    return Fa * Fb + (1.0 / 3.0) * bar(Fa) * bar(Fb);
}
double Shape::t2F(double t) {
    if(t >= 1e5) return 0;
    return A + B * exp(-pow(t / T, n));
}
double Shape::F2t(double F) {
    if(F <= A + EPS) return 1e9;
    return T * pow(-log((F - A) / B), 1.0 / n);
}

double Shape::pass_tao(double F) {
    return t2F(F2t(F) + tao);
}

void Shape::check_valid() {
    if(node_mem_range[0].second.size() != 1) {
        cerr << "the usage of memory at src should be 1, but " << node_mem_range[0].second.size() << endl;
        exit(1);
    }
    if(node_mem_range.back().second.size() != 1) {
        cerr << "the usage of memory at dst should be 1, but " << node_mem_range.back().second.size() << endl;
        exit(1);
    }
    for(int i = 1; i < (int)node_mem_range.size() - 1; i++) {
        int node_id = node_mem_range[i].first;
        vector<pair<int, int>> memory_range = node_mem_range[i].second;
        if(memory_range.size() != 2) {
            cerr << "node " << node_id << endl;
            cerr << "the usage of memory at med node should be 2, but " << memory_range.size() << endl;
            exit(1);
        }
        
        // memory_range[[0]and memory_range[[1]
        if(memory_range.front().second != memory_range[1].second) {
            cerr << "two range should finish same time (after swap)" << endl;
            exit(1);
        }
        assert(memory_range[0].first >= 0);
        assert(memory_range[1].first >= 0);
        assert(memory_range[0].first <= memory_range[0].second);
        assert(memory_range[1].first <= memory_range[1].second);
    }

    recursion_check(0, node_mem_range.size() - 1);
}

void Shape::recursion_check(int left, int right) {
    if(left == right - 1) {
        int left_start = node_mem_range[left].second.back().first;
        int left_end = node_mem_range[left].second.back().second;
        int right_start = node_mem_range[right].second.front().first;
        int right_end = node_mem_range[right].second.front().second;
        assert(left_start == right_start);
        assert(left_end - 1 >= left_start);
        assert(right_end - 1 >= right_start);
        return;
    }

    int latest = left + 1;
    for(int i = left + 1; i < right; i++) {
        if(node_mem_range[i].second[0].second > node_mem_range[latest].second[0].second) {
            latest = i;
        }
    }

    int latest_finish_time = node_mem_range[latest].second[0].second;
    int left_finish_time = node_mem_range[left].second.back().second;
    int left_start_time = node_mem_range[left].second.back().first;
    int right_finish_time = node_mem_range[right].second.front().second;
    int right_start_time = node_mem_range[right].second.front().first;
    if(latest_finish_time >= left_finish_time || latest_finish_time >= right_finish_time) {
        cerr << "swap time error, swap after finish" << endl;
        cerr << "swap time = " << latest_finish_time << endl;
        cerr << "left fin = " << left_finish_time << endl;
        cerr << "right fin = " << right_finish_time << endl;
        exit(1);
    }

    if(latest_finish_time < left_start_time || latest_finish_time < right_start_time) {
        cerr << "swap time error, swap before entangle" << endl;
        exit(1);
    }

    recursion_check(left, latest);
    recursion_check(latest, right);
}

void Shape::print() {
    cerr << "length = " << node_mem_range.size() << endl;
    for(int i = 0; i < (int)node_mem_range.size(); i++) {
        int id = node_mem_range[i].first;
        cerr << "--- id = " << id << " ---" << endl;
        for(auto P : node_mem_range[i].second) {
            cerr << "{" << P.first << ", " << P.second << "}" << endl;
        }
    }
}

bool Shape::operator< (Shape &ls) {
    return node_mem_range < ls.node_mem_range;
}

bool Shape::operator== (Shape &ls) {
    return node_mem_range == ls.node_mem_range;
}