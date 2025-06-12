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

  static double rand_range_double(double l, double r)
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

  std::mt19937 rng_engine;
  std::normal_distribution<> normal_dist;
}

inline void normalize(vector<double>& vec)
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

struct HyperParameters
{
public:
  int value_1, value_2, value_3;
  double value_4;
  int value_5;

  HyperParameters() {
    value_1 = 0;
    value_2 = 0;
    value_3 = 0;
    value_4 = 0;
    value_5 = 0;
  }

  HyperParameters(int v1, int v2, int v3, double v4, int v5) {
    value_1 = v1;
    value_2 = v2;
    value_3 = v3;
    value_4 = v4;
    value_5 = v5;
  }

  bool operator==(const HyperParameters& other) const {
    return value_1 == other.value_1 && value_2 == other.value_2 && value_3 == other.value_3 &&
      value_4 == other.value_4 && value_5 == other.value_5;
  }
};

struct HyperKeys
{
  int s_key;
  int l_key;

  HyperKeys(int s, int l) {
    l_key = l;
    s_key = s;
  }

  bool operator<(const HyperKeys& other) const {
    if (s_key != other.s_key) {
      return s_key < other.s_key;
    }
    return l_key < other.l_key;
  }
};

map<HyperKeys, HyperParameters> hyper_parameters = {
  {HyperKeys(1, 1), HyperParameters(475, 471, 127, 0.518192, 10)},
  {HyperKeys(1, 2), HyperParameters(475, 471, 127, 0.518192, 10)},
  {HyperKeys(1, 3), HyperParameters(475, 471, 134, 0.518192, 10)},
  {HyperKeys(1, 4), HyperParameters(475, 471, 134, 0.518192, 10)},
  {HyperKeys(1, 5), HyperParameters(475, 471, 134, 0.518192, 10)},
  {HyperKeys(4, 1), HyperParameters(421, 402, 149, 0.583516, 10)},
  {HyperKeys(4, 2), HyperParameters(421, 402, 149, 0.583516, 10)},
  {HyperKeys(4, 3), HyperParameters(421, 402, 149, 0.583516, 10)},
  {HyperKeys(4, 4), HyperParameters(419, 398, 173, 0.608158, 10)},
  {HyperKeys(4, 5), HyperParameters(419, 398, 173, 0.608158, 10)},
  {HyperKeys(9, 1), HyperParameters(446, 398, 173, 0.608158, 10)},
  {HyperKeys(9, 2), HyperParameters(446, 398, 173, 0.608158, 10)},
  {HyperKeys(9, 3), HyperParameters(446, 398, 173, 0.608158, 10)},
  {HyperKeys(9, 4), HyperParameters(419, 374, 173, 0.608158, 10)},
  {HyperKeys(9, 5), HyperParameters(419, 374, 173, 0.608158, 10)},
  {HyperKeys(16, 1), HyperParameters(451, 379, 315, 0.598233, 10)},
  {HyperKeys(16, 2), HyperParameters(451, 379, 315, 0.598233, 10)},
  {HyperKeys(16, 3), HyperParameters(451, 379, 315, 0.598233, 10)},
  {HyperKeys(16, 4), HyperParameters(451, 379, 315, 0.598233, 10)},
  {HyperKeys(16, 5), HyperParameters(451, 379, 315, 0.598233, 10)},
  {HyperKeys(25, 1), HyperParameters(458, 375, 315, 0.650053, 10)},
  {HyperKeys(25, 2), HyperParameters(458, 375, 315, 0.650053, 10)},
  {HyperKeys(25, 3), HyperParameters(458, 375, 315, 0.650053, 10)},
  {HyperKeys(25, 4), HyperParameters(458, 375, 315, 0.650053, 20)},
  {HyperKeys(25, 5), HyperParameters(458, 375, 315, 0.650053, 20)},
  {HyperKeys(36, 1), HyperParameters(644, 518, 339, 0.704136, 10)},
  {HyperKeys(36, 2), HyperParameters(464, 334, 302, 0.630076, 20)},
  {HyperKeys(36, 3), HyperParameters(464, 334, 302, 0.630076, 20)},
  {HyperKeys(36, 4), HyperParameters(464, 334, 302, 0.630076, 20)},
  {HyperKeys(36, 5), HyperParameters(476, 334, 320, 0.630076, 20)},
  {HyperKeys(49, 1), HyperParameters(654, 518, 347, 0.679975, 10)},
  {HyperKeys(49, 2), HyperParameters(464, 287, 339, 0.623786, 20)},
  {HyperKeys(49, 3), HyperParameters(464, 287, 339, 0.623786, 20)},
  {HyperKeys(49, 4), HyperParameters(464, 287, 339, 0.623786, 20)},
  {HyperKeys(49, 5), HyperParameters(476, 334, 320, 0.630076, 20)},
  {HyperKeys(64, 1), HyperParameters(641, 431, 339, 0.770498, 10)},
  {HyperKeys(64, 2), HyperParameters(463, 276, 361, 0.700321, 20)},
  {HyperKeys(64, 3), HyperParameters(463, 272, 375, 0.709452, 20)},
  {HyperKeys(64, 4), HyperParameters(463, 284, 361, 0.746696, 20)},
  {HyperKeys(64, 5), HyperParameters(759, 572, 369, 0.757343, 10)},
  {HyperKeys(81, 1), HyperParameters(634, 375, 339, 0.75, 10)},
  {HyperKeys(81, 2), HyperParameters(634, 375, 353, 0.75, 20)},
  {HyperKeys(81, 3), HyperParameters(634, 375, 353, 0.75, 20)},
  {HyperKeys(81, 4), HyperParameters(634, 375, 333, 0.75, 20)},
  {HyperKeys(81, 5), HyperParameters(634, 375, 339, 0.75, 20)},
  {HyperKeys(100, 1), HyperParameters(634, 375, 353, 0.75, 20)},
  {HyperKeys(100, 2), HyperParameters(634, 375, 353, 0.75, 20)},
  {HyperKeys(100, 3), HyperParameters(969, 700, 538, 0.784865, 10)},
  {HyperKeys(100, 4), HyperParameters(977, 709, 534, 0.784865, 10)},
  {HyperKeys(100, 5), HyperParameters(980, 709, 616, 0.730521, 10)},
  {HyperKeys(121, 1), HyperParameters(980, 709, 546, 0.784865, 10)},
  {HyperKeys(121, 2), HyperParameters(980, 709, 546, 0.784865, 10)},
  {HyperKeys(121, 3), HyperParameters(1000, 709, 524, 0.730521, 10)},
  {HyperKeys(121, 4), HyperParameters(980, 709, 613, 0.730521, 10)},
  {HyperKeys(121, 5), HyperParameters(980, 709, 613, 0.730521, 10)},
  {HyperKeys(144, 1), HyperParameters(987, 667, 602, 0.804736, 10)},
  {HyperKeys(144, 2), HyperParameters(987, 667, 602, 0.804736, 10)},
  {HyperKeys(144, 3), HyperParameters(993, 666, 613, 0.730521, 10)},
  {HyperKeys(144, 4), HyperParameters(993, 666, 613, 0.730521, 10)},
  {HyperKeys(144, 5), HyperParameters(982, 667, 602, 0.79957, 10)},
  {HyperKeys(169, 1), HyperParameters(993, 666, 613, 0.730521, 10)},
  {HyperKeys(169, 2), HyperParameters(987, 667, 602, 0.804736, 10)},
  {HyperKeys(169, 3), HyperParameters(987, 619, 599, 0.766846, 10)},
  {HyperKeys(169, 4), HyperParameters(993, 589, 617, 0.730521, 10)},
  {HyperKeys(169, 5), HyperParameters(993, 589, 617, 0.730521, 10)},
  {HyperKeys(196, 1), HyperParameters(1000, 582, 636, 0.798162, 10)},
  {HyperKeys(196, 2), HyperParameters(1000, 582, 636, 0.798162, 10)},
  {HyperKeys(196, 3), HyperParameters(1000, 582, 636, 0.798162, 10)},
  {HyperKeys(196, 4), HyperParameters(1000, 582, 636, 0.798162, 10)},
  {HyperKeys(196, 5), HyperParameters(1000, 582, 636, 0.798162, 10)},
  {HyperKeys(225, 1), HyperParameters(985, 552, 646, 0.752293, 10)},
  {HyperKeys(225, 2), HyperParameters(985, 552, 646, 0.752293, 10)},
  {HyperKeys(225, 3), HyperParameters(1000, 500, 645, 0.75, 10)},
  {HyperKeys(225, 4), HyperParameters(1000, 500, 645, 0.75, 10)},
  {HyperKeys(225, 5), HyperParameters(1000, 500, 645, 0.75, 10)},
  {HyperKeys(256, 1), HyperParameters(985, 552, 646, 0.752293, 10)},
  {HyperKeys(256, 2), HyperParameters(985, 552, 645, 0.75, 10)},
  {HyperKeys(256, 3), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(256, 4), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(256, 5), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(289, 1), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(289, 2), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(289, 3), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(289, 4), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(289, 5), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(324, 1), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(324, 2), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(324, 3), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(324, 4), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(324, 5), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(361, 1), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(361, 2), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(361, 3), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(361, 4), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(361, 5), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(400, 1), HyperParameters(1000, 406, 1000, 0.810731, 10)},
  {HyperKeys(400, 2), HyperParameters(1000, 406, 994, 0.824366, 10)},
  {HyperKeys(400, 3), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(400, 4), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(400, 5), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(441, 1), HyperParameters(997, 406, 1004, 0.806231, 10)},
  {HyperKeys(441, 2), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(441, 3), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(441, 4), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(441, 5), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(484, 1), HyperParameters(997, 406, 1004, 0.769124, 10)},
  {HyperKeys(484, 2), HyperParameters(997, 406, 1004, 0.769124, 10)},
  {HyperKeys(484, 3), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(484, 4), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(484, 5), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(529, 1), HyperParameters(1000, 383, 993, 0.800669, 10)},
  {HyperKeys(529, 2), HyperParameters(1000, 457, 1000, 0.75, 10)},
  {HyperKeys(529, 3), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(529, 4), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(529, 5), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(576, 1), HyperParameters(1000, 408, 1000, 0.779401, 10)},
  {HyperKeys(576, 2), HyperParameters(1000, 408, 1000, 0.779401, 10)},
  {HyperKeys(576, 3), HyperParameters(1000, 408, 1000, 0.779401, 10)},
  {HyperKeys(576, 4), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(576, 5), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(625, 1), HyperParameters(1000, 408, 1000, 0.779401, 10)},
  {HyperKeys(625, 2), HyperParameters(1000, 408, 1000, 0.779401, 10)},
  {HyperKeys(625, 3), HyperParameters(1000, 408, 1000, 0.779401, 10)},
  {HyperKeys(625, 4), HyperParameters(1000, 497, 1003, 0.75, 10)},
  {HyperKeys(625, 5), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(676, 1), HyperParameters(1000, 408, 1000, 0.779401, 10)},
  {HyperKeys(676, 2), HyperParameters(1000, 408, 1000, 0.779401, 10)},
  {HyperKeys(676, 3), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(676, 4), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(676, 5), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(729, 1), HyperParameters(999, 335, 1012, 0.77592, 10)},
  {HyperKeys(729, 2), HyperParameters(1000, 408, 1000, 0.779401, 10)},
  {HyperKeys(729, 3), HyperParameters(1000, 433, 1000, 0.786346, 10)},
  {HyperKeys(729, 4), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(729, 5), HyperParameters(1000, 433, 1000, 0.786346, 10)},
  {HyperKeys(784, 1), HyperParameters(999, 335, 1012, 0.77592, 10)},
  {HyperKeys(784, 2), HyperParameters(999, 335, 1012, 0.77592, 10)},
  {HyperKeys(784, 3), HyperParameters(1000, 395, 1000, 0.786346, 10)},
  {HyperKeys(784, 4), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(784, 5), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(841, 1), HyperParameters(999, 335, 1012, 0.77592, 10)},
  {HyperKeys(841, 2), HyperParameters(999, 335, 1005, 0.77592, 10)},
  {HyperKeys(841, 3), HyperParameters(1000, 395, 1000, 0.786346, 10)},
  {HyperKeys(841, 4), HyperParameters(991, 420, 1003, 0.732333, 10)},
  {HyperKeys(841, 5), HyperParameters(1000, 500, 1000, 0.75, 10)},
  {HyperKeys(900, 1), HyperParameters(999, 335, 1012, 0.77592, 10)},
  {HyperKeys(900, 2), HyperParameters(999, 335, 1005, 0.77592, 10)},
  {HyperKeys(900, 3), HyperParameters(1000, 395, 1000, 0.786346, 10)},
  {HyperKeys(900, 4), HyperParameters(1000, 419, 1000, 0.75, 10)},
  {HyperKeys(900, 5), HyperParameters(1000, 500, 1000, 0.75, 10)},
};

constexpr int INT_INF = 1e9;
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
  int top_y;
  int top_x;

  Layout(int l) {
    L = l;
    p.resize(l, vector<int>(l, 0));
    cost = 0;
    top_y = L / 2;
    top_x = L / 2;
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

  void layout_2() {
    int minimum_sum = INT_INF;
    int minimum_x = -1;
    int minimum_y = -1;
    for (int i = 0; i < L; i++) {
      for (int j = 0; j < L; j++) {
        int sum = 0;
        for (int k = 0; k < N; k++) {
          int dy = min(abs(i - Y[k]), L - abs(i - Y[k]));
          int dx = min(abs(j - X[k]), L - abs(j - X[k]));
          sum += dy + dx;
        }
        if (sum < minimum_sum) {
          minimum_sum = sum;
          minimum_x = j;
          minimum_y = i;
        }
      }
    }

    top_x = minimum_x;
    top_y = minimum_y;

    for (int i = 0; i < L; i++) {
      for (int j = 0; j < L; j++) {
        if (i == top_y && j == top_x) {
          p[i][j] = hyper_parameters[HyperKeys(S, L / 10)].value_1;
        }
        else {
          int dy = min(abs(i - top_y), L - abs(i - top_y));
          int dx = min(abs(j - top_x), L - abs(j - top_x));

          int val_2 = hyper_parameters[HyperKeys(S, L / 10)].value_2;
          int val_3 = hyper_parameters[HyperKeys(S, L / 10)].value_3;
          p[i][j] = max(0, val_2 - (dy + dx - 1) * val_3 / L);
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

  void reset() {
    query_count = 0;
    cost = 0;
  }

  int query(int i, int y, int x, const Layout& layout) {
    if (query_count >= MAX_Q) {
      return -1;
    }

    int qy = (Y[A[i]] + y + L) % L;
    int qx = (X[A[i]] + x + L) % L;
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

  void update_prob(int index_in, int y, int x, int m, const Layout& layout, int dir) {
    y = (y + L) % L;
    x = (x + L) % L;
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

    for(int times = 0; times < 1; times++) {
      // 行方向に正規化
      for (int i = 0; i < N; i++) {
        normalize_matrix(prob, 0, i);
      }
      // 列方向に正規化
      for (int j = 0; j < N; j++) {
        normalize_matrix(prob, 1, j);
      }
    }
    if (dir == 0) {
      // 行方向に正規化
      for (int i = 0; i < N; i++) {
        normalize_matrix(prob, 0, i);
      }
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

void generate_test_case(int l, int n, int s)
{
  L = l;
  N = n;
  S = s;

  local_case.init(N);

  // Y,Xの初期化
  set<pair<int, int>> YX;
  for (int i = 0; i < N; i++) {
    int y = rand_range(0, L - 1);
    int x = rand_range(0, L - 1);
    while (YX.count(make_pair(y, x))) {
      y = rand_range(0, L - 1);
      x = rand_range(0, L - 1);
    }
    YX.insert(make_pair(y, x));
  }
  for (int i = 0; i < N; i++) {
    auto it = YX.begin();
    advance(it, i);
    Y[i] = it->first;
    X[i] = it->second;
  }

  // Aの初期化
  for (int i = 0; i < N; i++) {
    local_case.A[i] = i;
  }
  shuffle_array(local_case.A.data(), N);

  // fの初期化
  std::normal_distribution<>::param_type param(0.0, S);
  normal_dist.param(param);
  for (int i = 0; i < MAX_Q; i++) {
    // 平均0, 標準偏差S の正規分布から値をサンプリングする
    local_case.f[i] = round(normal_dist(rng_engine));
  }
}

void input_data(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (exec_mode == 3) {
    // テストケースジェネレート済のため何もしない
    local_case.reset();
  }
  else if (!ifs.is_open()) {
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

  if (exec_mode != 3 && exec_mode != 1) {
    cerr
      << "cost = " << cost << ", "
      << "local_case.cost = " << local_case.cost << ", "
      << "layout.cost = " << layout.cost << ", "
      << "seikai_count = " << seikai_count << " / " << N << endl;
  }

  res /= cost;
  return ceil(res);
}

void open_ofs(int case_num, ofstream& ofs)
{
  if (exec_mode != 0 && exec_mode != 3) {
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
  else if (exec_mode != 3) {
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
  else if (exec_mode != 3) {
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
    if (exec_mode != 3) {
      ofs << i << " " << y << " " << x << endl;
    }
    m = local_case.query(i, y, x, layout);
  }
  return m;
}

void estimate_1(ofstream& ofs, const Layout& layout, Estimation& estimation) {
  double threshold = hyper_parameters[HyperKeys(S, L / 10)].value_4;
  int max_tate_length = hyper_parameters[HyperKeys(S, L / 10)].value_5;
  int iter = 0;

  // 列で見る
  vector<pair<int, int>> vp;
  for (int i = 0; i < N; i++) {
    int y = (layout.top_y - Y[i] + L) % L;
    if (abs(y - L) < y) {
      y = y - L;
    }
    int x = (layout.top_x - X[i] + L) % L;
    if (abs(x - L) < x) {
      x = x - L;
    }
    vp.push_back(make_pair(abs(x) + abs(y), i));
  }
  sort(vp.begin(), vp.end());
  for (auto it = vp.begin(); it != vp.end(); ++it) {
    int idx_out = it->second;
    int dist = it->first;
    if (dist > max_tate_length) {
      break;
    }

    while (iter < MAX_Q) {
      double max_prob = 0.0;
      int idx_in = 0;
      for (int i = 0; i < N; i++) {
        if (max_prob < estimation.prob[i][idx_out]) {
          max_prob = estimation.prob[i][idx_out];
          idx_in = i;
        }
      }
      if (max_prob > threshold) {
        break;
      }

      int y = (layout.top_y - Y[idx_out] + L) % L;
      if (abs(y - L) < y) {
        y = y - L;
      }
      int x = (layout.top_x - X[idx_out] + L) % L;
      if (abs(x - L) < x) {
        x = x - L;
      }
      int m = query(ofs, idx_in, y, x, layout);

      iter++;
      estimation.update_prob(idx_in, y, x, m, layout, 1);
    }
  }

  // 行で見る
  int idx_in = 0;
  int idx_out = 0;
  while (iter < MAX_Q) {
    if (estimation.get_minimum_prob() > threshold) {
      break;
    }
    for (int i = 0; i < N; i++) {
      idx_in = (idx_in + 1) % N;
      double max_prob = 0.0;
      for (int j = 0; j < N; j++) {
        if (max_prob < estimation.prob[idx_in][j]) {
          max_prob = estimation.prob[idx_in][j];
          idx_out = j;
        }
      }
      if (max_prob <= threshold) {
        break;
      }
    }

    int y = (layout.top_y - Y[idx_out] + L) % L;
    if (abs(y - L) < y) {
      y = y - L;
    }
    int x = (layout.top_x - X[idx_out] + L) % L;
    if (abs(x - L) < x) {
      x = x - L;
    }
    int m = query(ofs, idx_in, y, x, layout);

    iter++;
    estimation.update_prob(idx_in, y, x, m, layout, 0);
  }

  estimation.calculate_estimation();

  if (exec_mode != 3 && exec_mode != 1) {
    cerr << "iter = " << iter << endl;
  }
}

ll solve_case(int case_num)
{
  start_timer();

  input_data(case_num);

  Layout layout(L);
  layout.layout_2();

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

void output_hyper_parameters()
{
  std::ostringstream oss;
  oss << "./hyper_parameters.txt";
  ofstream ofs(oss.str());
  ofs << "map<HyperKeys, HyperParameters> hyper_parameters = {" << endl;
  for (auto it = hyper_parameters.begin(); it != hyper_parameters.end(); ++it) {
    ofs << "  {HyperKeys(" << it->first.s_key << ", " << it->first.l_key << "), "
      << "HyperParameters(" << it->second.value_1 << ", " << it->second.value_2 << ", "
      << it->second.value_3 << ", " << it->second.value_4 << ", " << it->second.value_5 << ")}," << endl;
  }
  ofs << "};" << endl;
}

int main()
{
  for (int _ = 0; _ < 0; _++) {
    rand_xorshift();
  }

  exec_mode = 1;

  if (exec_mode == 0) {
    solve_case(0);
  }
  else if (exec_mode <= 2) {
    ll sum_score = 0;
    for (int i = 0; i < 150; i++)
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
  else if (exec_mode == 3) {
    int iter = 0;
    queue<pair<HyperKeys, HyperParameters>> que;
    while (true)
    {
      iter++;
      // 100回戦わせて70勝以上したらハイパーパラメータを更新する
      int s = rand_range(1, 30);
      int ss = s * s;
      int l_key = rand_xorshift() % 5 + 1;

      HyperKeys h_key(ss, l_key);

      HyperParameters old_params = hyper_parameters[h_key];
      HyperParameters new_params = old_params;

      if (!que.empty()) {
        h_key = que.front().first;
        ss = h_key.s_key;
        for (int i = 1; i <= 30; i++) {
          if (i * i == ss) {
            s = i;
            break;
          }
        }
        l_key = h_key.l_key;
        old_params = hyper_parameters[h_key];
        new_params = que.front().second;
        que.pop();
      }
      else {
        int ra = rand_xorshift() % 20000;
        if (ra < 50) {
          new_params.value_1 = rand_range(1, 1000);
          new_params.value_2 = rand_range(1, 1000);
          new_params.value_3 = rand_range(1, 2000);
          new_params.value_4 = sqrt(rand_range_double(0.6, 0.98));
        }
        else if (ra < 100) {
          new_params.value_1 = rand_range(1, 1000);
        }
        else if (ra < 150) {
          new_params.value_2 = rand_range(1, 1000);
        }
        else if (ra < 200) {
          new_params.value_3 = rand_range(1, 2000);
        }
        else if (ra < 250) {
          new_params.value_4 = sqrt(rand_range_double(0.6, 0.98));
        }
        else if (ra < 400) {
          if (s == 1) {
            new_params = hyper_parameters[HyperKeys(2 * 2, l_key)];
          }
          else if (s == 30) {
            new_params = hyper_parameters[HyperKeys(29 * 29, l_key)];
          }
          else {
            if (rand_xorshift() % 2 == 0) {
              new_params = hyper_parameters[HyperKeys((s - 1) * (s - 1), l_key)];
            }
            else {
              new_params = hyper_parameters[HyperKeys((s + 1) * (s + 1), l_key)];
            }
          }
        }
        else if (ra < 900) {
          new_params.value_1 = min(1000, max(1, old_params.value_1 + (int)rand_range(-25, 25)));
        }
        else if (ra < 1400) {
          new_params.value_2 = min(1000, max(1, old_params.value_2 + (int)rand_range(-25, 25)));
        }
        else if (ra < 1900) {
          new_params.value_3 = min(2000, max(1, old_params.value_3 + (int)rand_range(-25, 25)));
        }
        else if (ra < 2400) {
          new_params.value_4 = min(0.99, old_params.value_4 + rand_range_double(-0.05, 0.05));
        }
        else if (ra < 2600) {
          new_params.value_5 = min(100, max(0, old_params.value_5 + (int)rand_range(-10, 10)));
        }
        else if (ra < 5000) {
          if (rand_xorshift() % 2 == 0) new_params.value_1 = min(1000, max(1, old_params.value_1 + (int)rand_range(-15, 15)));
          if (rand_xorshift() % 2 == 0)new_params.value_2 = min(1000, max(1, old_params.value_2 + (int)rand_range(-200, 200)));
          if (rand_xorshift() % 2 == 0)new_params.value_3 = min(2000, max(1, old_params.value_3 + (int)rand_range(-15, 15)));
          if (rand_xorshift() % 2 == 0)new_params.value_4 = min(0.99, old_params.value_4 + rand_range_double(-0.05, 0.05));
        }
        else if (ra < 10000) {
          if (rand_xorshift() % 2 == 0) new_params.value_1 = min(1000, max(1, old_params.value_1 + (int)rand_range(-150, 150)));
          if (rand_xorshift() % 2 == 0)new_params.value_2 = min(1500, max(1, old_params.value_2 + (int)rand_range(-300, 300)));
          if (rand_xorshift() % 2 == 0)new_params.value_3 = min(2000, max(1, old_params.value_3 + (int)rand_range(-150, 150)));
          if (rand_xorshift() % 2 == 0)new_params.value_4 = min(0.99, old_params.value_4 + rand_range_double(-0.15, 0.15));
        }
        else if (ra < 20000) {
          if (s == 1) {
            new_params = hyper_parameters[HyperKeys(2 * 2, l_key)];
          }
          else if (s == 30) {
            new_params = hyper_parameters[HyperKeys(29 * 29, l_key)];
          }
          else {
            if (rand_xorshift() % 2 == 0) {
              new_params = hyper_parameters[HyperKeys((s - 1) * (s - 1), l_key)];
            }
            else {
              new_params = hyper_parameters[HyperKeys((s + 1) * (s + 1), l_key)];
            }
          }
          if (rand_xorshift() % 2 == 0) new_params.value_1 = min(1000, max(1, new_params.value_1 + (int)rand_range(-15, 15)));
          if (rand_xorshift() % 2 == 0)new_params.value_2 = min(1500, max(1, new_params.value_2 + (int)rand_range(-200, 200)));
          if (rand_xorshift() % 2 == 0)new_params.value_3 = min(2000, max(1, new_params.value_3 + (int)rand_range(-15, 15)));
          if (rand_xorshift() % 2 == 0)new_params.value_4 = min(0.99, new_params.value_4 + rand_range_double(-0.05, 0.05));
        }
      }

      if (new_params == old_params) {
        continue;
      }

      int win_count = 0;
      int draw_count = 0;
      int lose_count = 0;
      for (int i = 0; i < 40; i++) {
        int l = min(50, l_key * 10 + (int)rand_range(0, 9));
        int n = rand_range(60, 100);
        generate_test_case(l, n, ss);

        ll old_score = solve_case(0);
        hyper_parameters[h_key] = new_params;
        ll new_score = solve_case(0);
        hyper_parameters[h_key] = old_params;

        if (new_score > old_score) {
          win_count++;
        }
        else if (new_score == old_score) {
          draw_count++;
        }
        else {
          lose_count++;
        }

        if (win_count <= 4 && lose_count >= 1) {
          break;
        }
        if (win_count >= 20) {
          break;
        }
        if (draw_count + lose_count >= 5) {
          break;
        }
        if (lose_count > win_count) {
          break;
        }
      }

      // 桁をそろえて出力
      cerr << "iter = " << iter << ", "
        << "s = " << std::setw(3) << ss << ", "
        << "l_key = " << std::setw(2) << l_key << ", "
        << "old_params = (" << std::setw(4) << old_params.value_1 << ", " << std::setw(3) << old_params.value_2 << ", " << std::setw(4) << old_params.value_3 << ", " << std::setw(8) << old_params.value_4 << ", " << std::setw(3) << old_params.value_5 << "), "
        << "new_params = (" << std::setw(4) << new_params.value_1 << ", " << std::setw(3) << new_params.value_2 << ", " << std::setw(4) << new_params.value_3 << ", " << std::setw(8) << new_params.value_4 << ", " << std::setw(3) << new_params.value_5 << "), "
        << "win_count = " << std::setw(2) << win_count << ", "
        << "draw_count = " << std::setw(2) << draw_count << ", "
        << "lose_count = " << std::setw(2) << lose_count << ", "
        << "time = " << get_elapsed_time() << endl;

      if (win_count >= 20) {
        hyper_parameters[h_key] = new_params;
        output_hyper_parameters();
        if (1 < s) {
          for (int i = 0; i < 20; i++) {
            que.push(make_pair(HyperKeys((s - 1) * (s - 1), l_key), new_params));
          }
        }
        if (s < 30) {
          for (int i = 0; i < 20; i++) {
            que.push(make_pair(HyperKeys((s + 1) * (s + 1), l_key), new_params));
          }
        }
        if (1 < l_key) {
          for (int i = 0; i < 20; i++) {
            que.push(make_pair(HyperKeys(ss, l_key - 1), new_params));
          }
        }
        if (l_key < 5) {
          for (int i = 0; i < 20; i++) {
            que.push(make_pair(HyperKeys(ss, l_key + 1), new_params));
          }
        }
      }

      if (win_count >= 15) {
        for (int i = 0; i < 10; i++) {
          HyperParameters new_new_params = new_params;
          if (rand_xorshift() % 2 == 0) new_new_params.value_1 = min(1000, max(1, new_params.value_1 + (int)rand_range(-15, 15)));
          if (rand_xorshift() % 2 == 0) new_new_params.value_2 = min(1000, max(1, new_params.value_2 + (int)rand_range(-50, 50)));
          if (rand_xorshift() % 2 == 0) new_new_params.value_3 = min(2000, max(1, new_params.value_3 + (int)rand_range(-15, 15)));
          if (rand_xorshift() % 2 == 0) new_new_params.value_4 = min(0.99, new_params.value_4 + rand_range_double(-0.05, 0.05));
          que.push(make_pair(h_key, new_new_params));
        }
        if (1 < s) {
          for (int i = 0; i < 10; i++) {
            HyperParameters new_new_params = new_params;
            if (rand_xorshift() % 2 == 0) new_new_params.value_1 = min(1000, max(1, new_params.value_1 + (int)rand_range(-15, 15)));
            if (rand_xorshift() % 2 == 0) new_new_params.value_2 = min(1000, max(1, new_params.value_2 + (int)rand_range(-50, 50)));
            if (rand_xorshift() % 2 == 0) new_new_params.value_3 = min(2000, max(1, new_params.value_3 + (int)rand_range(-15, 15)));
            if (rand_xorshift() % 2 == 0) new_new_params.value_4 = min(0.99, new_params.value_4 + rand_range_double(-0.05, 0.05));
            que.push(make_pair(HyperKeys((s - 1) * (s - 1), l_key), new_new_params));
          }
        }
        if (s < 30) {
          for (int i = 0; i < 10; i++) {
            HyperParameters new_new_params = new_params;
            if (rand_xorshift() % 2 == 0) new_new_params.value_1 = min(1000, max(1, new_params.value_1 + (int)rand_range(-15, 15)));
            if (rand_xorshift() % 2 == 0) new_new_params.value_2 = min(1000, max(1, new_params.value_2 + (int)rand_range(-50, 50)));
            if (rand_xorshift() % 2 == 0) new_new_params.value_3 = min(2000, max(1, new_params.value_3 + (int)rand_range(-15, 15)));
            if (rand_xorshift() % 2 == 0) new_new_params.value_4 = min(0.99, new_params.value_4 + rand_range_double(-0.05, 0.05));
            que.push(make_pair(HyperKeys((s + 1) * (s + 1), l_key), new_new_params));
          }
        }
        if (1 < l_key) {
          for (int i = 0; i < 10; i++) {
            HyperParameters new_new_params = new_params;
            if (rand_xorshift() % 2 == 0) new_new_params.value_1 = min(1000, max(1, new_params.value_1 + (int)rand_range(-15, 15)));
            if (rand_xorshift() % 2 == 0) new_new_params.value_2 = min(1000, max(1, new_params.value_2 + (int)rand_range(-50, 50)));
            if (rand_xorshift() % 2 == 0) new_new_params.value_3 = min(2000, max(1, new_params.value_3 + (int)rand_range(-15, 15)));
            if (rand_xorshift() % 2 == 0) new_new_params.value_4 = min(0.99, new_params.value_4 + rand_range_double(-0.05, 0.05));
            que.push(make_pair(HyperKeys(ss, l_key - 1), new_params));
          }
        }
        if (l_key < 5) {
          for (int i = 0; i < 10; i++) {
            HyperParameters new_new_params = new_params;
            if (rand_xorshift() % 2 == 0) new_new_params.value_1 = min(1000, max(1, new_params.value_1 + (int)rand_range(-15, 15)));
            if (rand_xorshift() % 2 == 0) new_new_params.value_2 = min(1000, max(1, new_params.value_2 + (int)rand_range(-50, 50)));
            if (rand_xorshift() % 2 == 0) new_new_params.value_3 = min(2000, max(1, new_params.value_3 + (int)rand_range(-15, 15)));
            if (rand_xorshift() % 2 == 0) new_new_params.value_4 = min(0.99, new_params.value_4 + rand_range_double(-0.05, 0.05));
            que.push(make_pair(HyperKeys(ss, l_key + 1), new_params));
          }
        }
      }
    }
  }

  return 0;
}
