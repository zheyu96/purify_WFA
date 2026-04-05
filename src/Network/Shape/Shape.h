#ifndef __SHAPE_H
#define __SHAPE_H

#include <vector>
#include <utility>
#include <map>
#include <cmath>
#include <iostream>
#include <cassert>

using namespace std;

// 定義 Shape_vector
typedef vector<pair<int, vector<pair<int, int>>>> Shape_vector;

class Shape {
private:
    Shape_vector node_mem_range;
    double A, B, n, T, tao;

    // [新增] 必須要在這裡宣告這個變數，Shape.cpp 才能使用它
    bool purification_enabled = false;

    // 每條 link 的 purification 輪數（大小 = node數-1），空代表不啟用
    vector<int> link_purify_rounds;

    double recursion_get_fidelity(int left, int right, const map<pair<int, int> , double> &F_init);
    // void check_valid(); // <--- 原本在這裡 (private)，導致錯誤
    void recursion_check(int left, int right);
    
    // 輔助計算函式
    double bar(double F);
    double Fswap(double Fa, double Fb);
    double t2F(double t);
    double F2t(double F);
    double pass_tao(double F);

public:
    Shape(Shape_vector _node_mem_range);
    Shape(Shape_vector _node_mem_range, vector<int> _link_purify_rounds);
    Shape();

    const Shape_vector& get_node_mem_range() const;
    const vector<int>& get_link_purify_rounds() const;

    // [修改] 這裡的宣告必須包含新的參數 bool enable_purification
    double get_fidelity(double _A, double _B, double _n, double _T, double _tao,
                        const map<pair<int, int> , double>& F_init, bool enable_purification = false);

    // [修正] 將 check_valid 移動到 public 區域
    void check_valid(); 

    void print();
    bool operator< (Shape &ls);
    bool operator== (Shape &ls);
};

#endif