#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cctype>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <deque>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <utility>
#include <vector>

using namespace std;

typedef long long int ll;

// タイマー
namespace
{
  std::chrono::steady_clock::time_point start_time_clock;

  void start_timer()
  {
    start_time_clock = std::chrono::steady_clock::now();
  }

  double get_elapsed_time()
  {
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time_clock;
    return elapsed.count();
  }
}

// 乱数
namespace
{
  static uint32_t rand_xorshift()
  {
    static uint32_t x = 123456789;
    static uint32_t y = 362436069;
    static uint32_t z = 521288629;
    static uint32_t w = 88675123;
    uint32_t t = x ^ (x << 11);
    x = y;
    y = z;
    z = w;
    w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
    return w;
  }

  static double rand_01()
  {
    return (rand_xorshift() + 0.5) * (1.0 / UINT_MAX);
  }

  static double rand_range(double l, double r)
  {
    return l + (r - l) * rand_01();
  }

  static uint32_t rand_range(uint32_t l, uint32_t r)
  {
    return l + rand_xorshift() % (r - l + 1); // [l, r]
  }

  void shuffle_array(int* arr, int n)
  {
    for (int i = n - 1; i >= 0; i--) {
      int j = rand_xorshift() % (i + 1);
      int swa = arr[i];
      arr[i] = arr[j];
      arr[j] = swa;
    }
  }
}

inline void normalize(vector<double> vec)
{
  double sum = 0.0;
  for (int i = 0; i < vec.size(); i++) {
    sum += vec[i];
  }
  for (int i = 0; i < vec.size(); i++) {
    vec[i] /= sum;
  }
}

inline void normalize_matrix(vector<vector<double>>& mat, int dir, int index)
{
  if (dir == 0) {
    // 行方向の正規化
    double sum = 0.0;
    for (int j = 0; j < mat[index].size(); j++) {
      sum += mat[index][j];
    }
    for (int j = 0; j < mat[index].size(); j++) {
      mat[index][j] /= sum;
    }
  }
  else if (dir == 1) {
    // 列方向の正規化
    double sum = 0.0;
    for (int i = 0; i < mat.size(); i++) {
      sum += mat[i][index];
    }
    for (int i = 0; i < mat.size(); i++) {
      mat[i][index] /= sum;
    }
  }
}

const double TIME_LIMIT = 3.8;
int exec_mode;

constexpr int MAX_L = 50;
constexpr int MAX_N = 100;
constexpr int MAX_Q = 10000;

int L, N, S;
int Y[MAX_N], X[MAX_N];

class Layout
{
public:
  int L;
  vector<vector<int>> p;
  ll cost;
  int mid_y;
  int mid_x;

  Layout(int l) {
    L = l;
    p.resize(l, vector<int>(l, 0));
    cost = 0;
    mid_y = L / 2;
    mid_x = L / 2;
  }

  void calc_cost() {
    cost = 0;
    for (int i = 0; i < L; i++) {
      for (int j = 0; j < L; j++) {
        cost += (p[i][j] - p[(i + 1) % L][j]) * (p[i][j] - p[(i + 1) % L][j]);
        cost += (p[i][j] - p[i][(j + 1) % L]) * (p[i][j] - p[i][(j + 1) % L]);
      }
    }
  }

  void layout_1() {
    for (int i = 0; i < L; i++) {
      for (int j = 0; j < L; j++) {
        if (i == mid_y && j == mid_x) {
          p[i][j] = 1000;
        }
        else {
          p[i][j] = max(0, 500 - (abs(mid_y - i) + abs(mid_x - j)) * 1000 / L);
        }
      }
    }
    calc_cost();
  }
};

class LocalCase
{
public:
  int query_count;
  ll cost = 0;
  vector<int> A;
  vector<int> f;

  void init(int n) {
    query_count = 0;
    cost = 0;
    A.resize(n, 0);
    f.resize(MAX_Q, 0);
  }

  int query(int i, int y, int x, const Layout& layout) {
    if (query_count >= MAX_Q) {
      return -1;
    }

    int qy = (Y[A[i]] + y) % L;
    int qx = (X[A[i]] + x) % L;
    cost += 100 * (10 + abs(y) + abs(x));
    query_count++;
    return layout.p[qy][qx] + f[query_count - 1];
  }
};
LocalCase local_case;

class Estimation
{
public:
  vector<int> E;
  vector<vector<double>> prob;
  vector<double> tmp_prob;

  Estimation(int n) {
    E.resize(n, 0);
    prob.resize(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        prob[i][j] = 1.0 / n;
      }
    }
    tmp_prob.resize(n, 0.0);
  }

  double calc_ruiseki(double num) {
    return 1.0 / 2.0 * (1.0 + erf(num / sqrt(2.0 * S * S)));
  }

  double calc_kakuritsu(int kitai, int output) {
    double res = 0.0;
    int diff = output - kitai;
    if (output == 0) {
      res = calc_ruiseki(diff + 0.5) - calc_ruiseki(diff - S * 10);
    }
    else if (output == 1000) {
      res = calc_ruiseki(diff + S * 10) - calc_ruiseki(diff - 0.5);
    }
    else {
      res = calc_ruiseki(diff + 0.5) - calc_ruiseki(diff - 0.5);
    }
    return res;
  }

  void update_prob(int index_in, int y, int x, int m, const Layout& layout) {
    for (int i = 0; i < N; i++) {
      int kitai = layout.p[(y + Y[i]) % L][(x + X[i]) % L];
      tmp_prob[i] = calc_kakuritsu(kitai, m);
    }
    normalize(tmp_prob);
    for (int j = 0; j < N; j++) {
      double other_prob = (1.0 - tmp_prob[j]) / (N - 1);
      for (int i = 0; i < N; i++) {
        if (i == index_in) {
          prob[i][j] *= tmp_prob[j];
        }
        else {
          prob[i][j] *= other_prob;
        }
      }
    }
    // 行方向に正規化
    for (int i = 0; i < N; i++) {
      normalize_matrix(prob, 0, i);
    }
    // 列方向に正規化
    for (int j = 0; j < N; j++) {
      normalize_matrix(prob, 1, j);
    }
    // 行方向に正規化
    for (int i = 0; i < N; i++) {
      normalize_matrix(prob, 0, i);
    }
  }

  double get_minimum_prob() {
    double res = 1.0;
    for (int i = 0; i < N; i++) {
      double max_prob = 0.0;
      for (int j = 0; j < N; j++) {
        max_prob = max(max_prob, prob[i][j]);
      }
      res = min(res, max_prob);
    }
    return res;
  }

  void calculate_estimation() {
    for (int i = 0; i < N; i++) {
      double max_prob = 0.0;
      int max_index = -1;
      for (int j = 0; j < N; j++) {
        if (prob[i][j] > max_prob) {
          max_prob = prob[i][j];
          max_index = j;
        }
      }
      E[i] = max_index;
    }
  }
};

void input_data(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    cin >> L >> N >> S;
    for (int i = 0; i < N; i++) {
      cin >> Y[i] >> X[i];
    }
  }
  else {
    // ファイル入力
    ifs >> L >> N >> S;
    for (int i = 0; i < N; i++) {
      ifs >> Y[i] >> X[i];
    }

    local_case.init(N);
    for (int i = 0; i < N; i++) {
      ifs >> local_case.A[i];
    }
    for (int i = 0; i < MAX_Q; i++) {
      ifs >> local_case.f[i];
    }
  }
}

ll calculate_local_score(const Layout& layout, const Estimation& estimation)
{
  if (exec_mode == 0) {
    return 0;
  }

  double res = 1e14;
  int seikai_count = 0;
  for (int i = 0; i < N; i++) {
    if (estimation.E[i] == local_case.A[i]) {
      seikai_count++;
    }
    else {
      res *= 0.8;
    }
  }
  double cost = layout.cost + local_case.cost + 1e5;
  cerr
    << "cost = " << cost << ", "
    << "local_case.cost = " << local_case.cost << ", "
    << "layout.cost = " << layout.cost << ", "
    << "seikai_count = " << seikai_count << " / " << N << endl;
  res /= cost;
  return ceil(res);
}

void open_ofs(int case_num, ofstream& ofs)
{
  if (exec_mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofs.open(oss.str());
  }
}

void output_layout(ofstream& ofs, const Layout& layout)
{
  if (exec_mode == 0) {
    // 標準出力
    for (int i = 0; i < L; i++) {
      for (int j = 0; j < L; j++) {
        cout << layout.p[i][j] << " ";
      }
      cout << endl;
    }
    cout.flush();
  }
  else {
    // ファイル出力
    for (int i = 0; i < L; i++) {
      for (int j = 0; j < L; j++) {
        ofs << layout.p[i][j] << " ";
      }
      ofs << endl;
    }
  }
}

void output_estimation(ofstream& ofs, const Estimation& estimation)
{
  if (exec_mode == 0) {
    // 標準出力
    cout << "-1 -1 -1" << endl;
    for (int i = 0; i < N; i++) {
      cout << estimation.E[i] << endl;
    }
    cout.flush();
  }
  else {
    // ファイル出力
    ofs << "-1 -1 -1" << endl;
    for (int i = 0; i < N; i++) {
      ofs << estimation.E[i] << endl;
    }
  }
}

int query(ofstream& ofs, int i, int y, int x, const Layout& layout)
{
  int m = 0;
  if (exec_mode == 0) {
    // 標準出力
    cout << i << " " << y << " " << x << endl;
    cout.flush();
    cin >> m;
  }
  else {
    // ファイル出力
    ofs << i << " " << y << " " << x << endl;
    m = local_case.query(i, y, x, layout);
  }
  return m;
}

void estimate_1(ofstream& ofs, const Layout& layout, Estimation& estimation) {
  int iter = 0;
  while (iter < MAX_Q) {
    iter++;
    int ra_in = rand_xorshift() % N;
    int ra_out = rand_xorshift() % N;

    int y = (layout.mid_y - Y[ra_out] + L) % L;
    int x = (layout.mid_x - X[ra_out] + L) % L;
    int m = query(ofs, ra_in, y, x, layout);

    estimation.update_prob(ra_in, y, x, m, layout);

    if (estimation.get_minimum_prob() > 0.90) {
      break;
    }
  }

  estimation.calculate_estimation();
  cerr << "iter = " << iter << endl;
  //for(int i = 0; i < N; i++) {
  //  for(int j = 0; j < N; j++) {
  //    cerr << estimation.prob[i][j] << " ";
  //  }
  //  cerr << endl;
  //}
}

ll solve_case(int case_num)
{
  start_timer();

  input_data(case_num);

  Layout layout(L);
  layout.layout_1();

  ofstream ofs;
  open_ofs(case_num, ofs);
  output_layout(ofs, layout);

  Estimation estimation(N);

  estimate_1(ofs, layout, estimation);

  output_estimation(ofs, estimation);

  ll score = 0;
  if (exec_mode != 0) {
    score = calculate_local_score(layout, estimation);
  }
  return score;
}

int main()
{
  exec_mode = 2;

  if (exec_mode == 0) {
    solve_case(0);
  }
  else if (exec_mode <= 2) {
    ll sum_score = 0;
    for (int i = 0; i < 15; i++)
    {
      ll score = solve_case(i);
      sum_score += score;
      if (exec_mode == 1) {
        cerr << score << endl;
      }
      else {
        cerr << "case = " << setw(2) << i << ", "
          << "score = " << setw(4) << score << ", "
          << "sum = " << setw(5) << sum_score << ", "
          << "time = " << setw(5) << get_elapsed_time() << ", "
          << endl;
      }
    }
  }

  return 0;
}
