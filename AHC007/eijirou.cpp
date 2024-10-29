//#pragma GCC target("avx2")
//#pragma GCC optimize("O3")
//#pragma GCC optimize("unroll-loops")
//
//#include <bits/stdc++.h>
//#include <atcoder/all>
//using namespace std;
//using namespace atcoder;
//
//double uniform(double l, double r) {
//  // l以上r以下の実数をランダムに生成する.
//  return l + (r - l) * (double)rand() / RAND_MAX;
//}
//
//struct Solver {
//  int N;
//  int M;
//  vector<int> x;
//  vector<int> y;
//  vector<int> u;
//  vector<int> v;
//  int repeat;
//  double min_cost;
//  double max_cost;
//  vector<int> distances;
//  vector<vector<pair<int, int>>> random_distances;
//  dsu uf;
//
//  int get_distance(int i, int j) {
//    // Euclid distance
//    return round(sqrt((x.at(i) - x.at(j)) * (x.at(i) - x.at(j)) + (y.at(i) - y.at(j)) * (y.at(i) - y.at(j))));
//  }
//
//  void init() {
//    N = 400;
//    M = 1995;
//    x = vector<int>(N);
//    y = vector<int>(N);
//    u = vector<int>(M);
//    v = vector<int>(M);
//    repeat = 160; // モンテカルロ法のシミュレーション回数
//    min_cost = 1.1; // ユークリッド距離にかける最小値
//    max_cost = 2.9; // ユークリッド距離にかける最大値
//    distances = vector<int>(M);
//    random_distances = vector<vector<pair<int, int>>>(repeat, vector<pair<int, int>>(M));
//    uf = dsu(N);
//  }
//
//  void input() {
//    for (int i = 0; i < N; i++) {
//      cin >> x.at(i) >> y.at(i);
//    }
//    for (int i = 0; i < M; i++) {
//      cin >> u.at(i) >> v.at(i);
//    }
//  }
//
//  void prepare() {
//    for (int i = 0; i < M; i++) {
//      distances.at(i) = get_distance(u.at(i), v.at(i));
//    }
//    for (int i = 0; i < repeat; i++) {
//      for (int j = 0; j < M; j++) {
//        random_distances.at(i).at(j) = make_pair(round(uniform(min_cost, max_cost) * distances.at(j)), j);
//      }
//      sort(random_distances.at(i).begin(), random_distances.at(i).end());
//    }
//  }
//
//  int kruskal(dsu uf, int start, int case_number) {
//    for (pair<int, int> p : random_distances.at(case_number)) {
//      int cost = p.first;
//      int idx = p.second;
//      if (start < idx) {
//        uf.merge(u.at(idx), v.at(idx));
//        if (uf.same(u.at(start), v.at(start))) {
//          return cost;
//        }
//      }
//    }
//    return -1;
//  }
//
//  bool monte_carlo(int start) {
//    int sum_costs = 0;
//    for (int i = 0; i < repeat; i++) {
//      sum_costs += kruskal(uf, start, i);
//      if (i == 0 && sum_costs == -1) {
//        return true;
//      }
//    }
//    return repeat * distances.at(start) < sum_costs;
//  }
//
//  void solve() {
//    for (int i = 0; i < M; i++) {
//      cin >> distances.at(i);
//      if (!uf.same(u.at(i), v.at(i)) && monte_carlo(i)) {
//        cout << 1 << endl << flush;
//        uf.merge(u.at(i), v.at(i));
//      }
//      else {
//        cout << 0 << endl << flush;
//      }
//    }
//  }
//};
//
//int main() {
//  Solver solver;
//  solver.init();
//  solver.input();
//  solver.prepare();
//  solver.solve();
//}