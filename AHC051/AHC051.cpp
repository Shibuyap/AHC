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
#include <memory>
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

typedef pair<int, int> P;
typedef long long int ll;

// タイマー
class Timer
{
private:
  std::chrono::steady_clock::time_point start_time_clock;

public:
  void start()
  {
    start_time_clock = std::chrono::steady_clock::now();
  }

  double get_elapsed_time()
  {
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time_clock;
    return elapsed.count();
  }
};

Timer timer;

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

// ユークリッド距離を計算する関数
double euclidean_distance(int x1, int y1, int x2, int y2)
{
  return sqrt(static_cast<double>((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)));
}

// 符号を返す関数
int sign(int x)
{
  return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}

// 3点の向きを計算する関数（座標版）
// a, b, cの3点に対して、aからbへの線分とaからcへの線分の外積の符号を返す
int orientation(int ax, int ay, int bx, int by, int cx, int cy)
{
  int cross = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
  return sign(cross);
}

// 2つの線分が交差するかどうかを判定する関数（座標版）
bool segments_intersect(int p1x, int p1y, int p2x, int p2y, int q1x, int q1y, int q2x, int q2y)
{
  // バウンディングボックスのチェック
  if (max(p1x, p2x) < min(q1x, q2x) ||
    max(q1x, q2x) < min(p1x, p2x) ||
    max(p1y, p2y) < min(q1y, q2y) ||
    max(q1y, q2y) < min(p1y, p2y)) {
    return false;
  }

  int o1 = orientation(p1x, p1y, p2x, p2y, q1x, q1y);
  int o2 = orientation(p1x, p1y, p2x, p2y, q2x, q2y);
  int o3 = orientation(q1x, q1y, q2x, q2y, p1x, p1y);
  int o4 = orientation(q1x, q1y, q2x, q2y, p2x, p2y);

  return (o1 * o2 <= 0) && (o3 * o4 <= 0);
}

const double TIME_LIMIT = 1.8;
int exec_mode;

const int MAX_V = 1100;

class Board
{
public:
  int n, m, k;
  vector<int> x, y; // 各頂点の座標
  vector<vector<int>> nearest; // 各頂点から別の頂点への近い順
  vector<vector<double>> p; // k * nの行列
};

class Place
{
public:
  int k;
  int v1, v2;
  int keep_v1, keep_v2;
  Place()
    : k(-1), v1(-1), v2(-1), keep_v1(-1), keep_v2(-1)
  {
  }
  Place(int k_, int v1_, int v2_, int keep_v1_, int keep_v2_)
    : k(k_), v1(v1_), v2(v2_), keep_v1(keep_v1_), keep_v2(keep_v2_)
  {
  }

  // 等価比較演算子
  bool operator==(const Place& other) const
  {
    return k == other.k && v1 == other.v1 && v2 == other.v2 &&
      keep_v1 == other.keep_v1 && keep_v2 == other.keep_v2;
  }
};

class State
{
public:
  vector<Place> places; // n個の処理装置、m個の分別器、1つの搬入口を含む

  vector<int> d;
  int s;

  int score;
};

struct GameData
{
  State state;
  Board board;
};

// 分別器評価関数の基底クラス
class SeparatorEvaluator
{
public:
  virtual ~SeparatorEvaluator() = default;

  // 分別器kを頂点vに配置した場合の評価値を計算
  // prob_dist: 頂点vに到達する各ごみ種類の確率分布
  // 返り値: 評価値（大きいほど良い）
  virtual double evaluate(const Board& b, int k, const vector<double>& prob_dist) const = 0;

  virtual string name() const = 0;
};

// エントロピー減少量による評価
class EntropyReductionEvaluator : public SeparatorEvaluator
{
public:
  double evaluate(const Board& b, int k, const vector<double>& prob_dist) const override
  {
    // 元のエントロピーを計算
    double original_entropy = 0.0;
    double total_prob = 0.0;
    for (int i = 0; i < b.n; i++) {
      if (prob_dist[i] > 1e-9) {
        total_prob += prob_dist[i];
        original_entropy -= prob_dist[i] * log2(prob_dist[i]);
      }
    }
    if (total_prob < 1e-9) return 0.0;
    original_entropy /= total_prob;
    original_entropy += log2(total_prob);

    // 分別後の各出口でのエントロピーを計算
    double entropy_after = 0.0;
    vector<double> prob_v1(b.n, 0.0), prob_v2(b.n, 0.0);
    double total_v1 = 0.0, total_v2 = 0.0;

    for (int i = 0; i < b.n; i++) {
      prob_v1[i] = prob_dist[i] * b.p[k][i];
      prob_v2[i] = prob_dist[i] * (1.0 - b.p[k][i]);
      total_v1 += prob_v1[i];
      total_v2 += prob_v2[i];
    }

    // v1のエントロピー
    if (total_v1 > 1e-9) {
      double entropy_v1 = 0.0;
      for (int i = 0; i < b.n; i++) {
        if (prob_v1[i] > 1e-9) {
          double p = prob_v1[i] / total_v1;
          entropy_v1 -= p * log2(p);
        }
      }
      entropy_after += (total_v1 / total_prob) * entropy_v1;
    }

    // v2のエントロピー
    if (total_v2 > 1e-9) {
      double entropy_v2 = 0.0;
      for (int i = 0; i < b.n; i++) {
        if (prob_v2[i] > 1e-9) {
          double p = prob_v2[i] / total_v2;
          entropy_v2 -= p * log2(p);
        }
      }
      entropy_after += (total_v2 / total_prob) * entropy_v2;
    }

    return original_entropy - entropy_after;
  }

  string name() const override { return "EntropyReduction"; }
};

// 分散による評価
class VarianceEvaluator : public SeparatorEvaluator
{
public:
  double evaluate(const Board& b, int k, const vector<double>& prob_dist) const override
  {
    // 各出口での確率分布を計算
    vector<double> prob_v1(b.n, 0.0), prob_v2(b.n, 0.0);
    double total_v1 = 0.0, total_v2 = 0.0;

    for (int i = 0; i < b.n; i++) {
      prob_v1[i] = prob_dist[i] * b.p[k][i];
      prob_v2[i] = prob_dist[i] * (1.0 - b.p[k][i]);
      total_v1 += prob_v1[i];
      total_v2 += prob_v2[i];
    }

    // 各出口での分布の分散を計算（正規化後）
    double variance_v1 = 0.0, variance_v2 = 0.0;

    if (total_v1 > 1e-9) {
      double mean_v1 = 1.0 / b.n; // 理想的な平均値
      for (int i = 0; i < b.n; i++) {
        double p = prob_v1[i] / total_v1;
        variance_v1 += (p - mean_v1) * (p - mean_v1);
      }
    }

    if (total_v2 > 1e-9) {
      double mean_v2 = 1.0 / b.n; // 理想的な平均値
      for (int i = 0; i < b.n; i++) {
        double p = prob_v2[i] / total_v2;
        variance_v2 += (p - mean_v2) * (p - mean_v2);
      }
    }

    // 分散が大きいほど分離性能が高い（特定のごみ種類に偏っている）
    return variance_v1 + variance_v2;
  }

  string name() const override { return "Variance"; }
};

class DirectedAcyclicGraph
{
public:
  // --- 1. 参照を受け取るコンストラクタ -----------------------------
  DirectedAcyclicGraph(State& state, Board& board)
    : state_(state), board_(board)
  {
  }

  // --- 2. コピー禁止（参照は付け替えられないので） -----------------
  DirectedAcyclicGraph(const DirectedAcyclicGraph&) = delete;
  DirectedAcyclicGraph& operator=(const DirectedAcyclicGraph&) = delete;

  //    ムーブは必要ならデフォルト許可
  DirectedAcyclicGraph(DirectedAcyclicGraph&&) = default;
  DirectedAcyclicGraph& operator=(DirectedAcyclicGraph&&) = default;

  // --- 3. 必要に応じて getter ---------------------------------------
  State& state()       noexcept { return state_; }
  const State& state() const noexcept { return state_; }

  Board& board()       noexcept { return board_; }
  const Board& board() const noexcept { return board_; }

  void RecreateTopologicalSort()
  {
    // g_ がまだ生成されていない場合に備えて
    if (g_.empty()) {
      RecreateG();
    }

    const int N = static_cast<int>(g_.size());

    // ---- 1.  indegree 計算 ----
    std::vector<int> indeg(N, 0);
    for (int u = 0; u < N; ++u) {
      for (int v : g_[u]) {
        ++indeg[v];
      }
    }

    // ---- 2.  キュー初期化（根を先頭に）----
    std::queue<int> q;
    q.push(board_.n + board_.m);          // もともと root として使っていた頂点

    // ---- 3.  トポロジカルソート ----
    topo_.clear();
    topo_.reserve(N);

    while (!q.empty()) {
      int v = q.front();
      q.pop();
      topo_.push_back(v);

      for (int nxt : g_[v]) {
        if (--indeg[nxt] == 0) {
          q.push(nxt);
        }
      }
    }

    // ---- 4.  topoPos_ を再構築 ----
    topoPos_.assign(N, -1);
    for (int i = 0; i < static_cast<int>(topo_.size()); ++i) {
      topoPos_[topo_[i]] = i;
    }
  }

  void RecreateG()
  {
    const int N = static_cast<int>(state_.places.size());
    g_.assign(N, {});                    // 全頂点を空で初期化

    const int root = board_.n + board_.m;
    if (root < 0 || root >= N) return;   // 念のため範囲チェック

    std::vector<char> vis(N, 0);
    std::queue<int> q;
    q.push(root);
    vis[root] = 1;

    while (!q.empty()) {
      const int u = q.front();
      q.pop();

      const int v1 = state_.places[u].v1;
      if (v1 != -1) {
        g_[u].push_back(v1);
        if (!vis[v1]) {
          vis[v1] = 1;
          q.push(v1);
        }
      }

      const int v2 = state_.places[u].v2;
      if (v2 != -1) {
        g_[u].push_back(v2);
        if (!vis[v2]) {
          vis[v2] = 1;
          q.push(v2);
        }
      }
    }
  }

  bool CanMaintainDAG(int vertex, int new_v1, int new_v2)
  {
    const int N = static_cast<int>(state_.places.size());
    const int root = board_.n + board_.m;

    int old_v1 = state_.places[vertex].v1;
    int old_v2 = state_.places[vertex].v2;

    state_.places[vertex].v1 = new_v1;
    state_.places[vertex].v2 = new_v2;

    vector<vector<int>> tmp_g(N);
    vector<int> indeg(N, 0);
    vector<char> reachable(N, 0);

    queue<int> q;
    q.push(root);
    reachable[root] = 1;

    while (!q.empty()) {
      int u = q.front();
      q.pop();

      int v1 = state_.places[u].v1;
      if (v1 != -1) {
        tmp_g[u].push_back(v1);
        indeg[v1]++;
        if (!reachable[v1]) {
          reachable[v1] = 1;
          q.push(v1);
        }
      }

      int v2 = state_.places[u].v2;
      if (v2 != -1) {
        tmp_g[u].push_back(v2);
        indeg[v2]++;
        if (!reachable[v2]) {
          reachable[v2] = 1;
          q.push(v2);
        }
      }
    }

    queue<int> topo_q;
    for (int i = 0; i < N; i++) {
      if (reachable[i] && indeg[i] == 0) {
        topo_q.push(i);
      }
    }

    int cnt = 0;
    while (!topo_q.empty()) {
      int u = topo_q.front();
      topo_q.pop();
      cnt++;

      for (int v : tmp_g[u]) {
        if (--indeg[v] == 0) {
          topo_q.push(v);
        }
      }
    }

    state_.places[vertex].v1 = old_v1;
    state_.places[vertex].v2 = old_v2;

    int reachable_cnt = 0;
    for (int i = 0; i < N; i++) {
      if (reachable[i]) reachable_cnt++;
    }

    return cnt == reachable_cnt;
  }

public:
  vector<vector<int>> g_;
  vector<int> topo_;
  vector<int> topoPos_;
  State& state_;
  Board& board_;
};

static GameData read_state(std::istream& is)
{
  State s;
  Board b;
  is >> b.n >> b.m >> b.k;

  b.x.resize(b.n + b.m + 1);
  b.y.resize(b.n + b.m + 1);
  s.places.resize(b.n + b.m + 1);

  auto read_places = [&](int cnt, int offset)
    {
      for (int i = 0; i < cnt; ++i) {
        int x, y;
        is >> x >> y;
        b.x[offset + i] = x;
        b.y[offset + i] = y;
        s.places[offset + i] = Place(0, -1, -1, -1, -1); // kは-1で初期化
      }
    };

  read_places(b.n, 0);
  read_places(b.m, b.n);
  s.places.back() = Place(0, -1, -1, -1, -1);

  b.p.assign(b.k, std::vector<double>(b.n));
  for (auto& row : b.p) {
    for (auto& v : row) {
      is >> v;
    }
  }

  GameData gd;
  gd.state = s;
  gd.board = b;
  return gd;
}

GameData input_data(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());
  if (ifs) {
    // ファイルが存在
    return read_state(ifs);
  }
  // 標準入力へフォールバック
  return read_state(std::cin);
}

static void write_state(std::ostream& os, const Board& b, const State& s)
{
  for (int i = 0; i < b.n; ++i) {
    os << s.d[i] << ' ';
  }
  os << '\n' << s.s << '\n';

  for (int i = b.n; i < b.n + b.m; ++i) {
    const auto& pl = s.places[i];
    if (pl.v1 == -1) {
      os << -1 << '\n';
    }
    else {
      os << pl.k << ' ' << pl.v1 << ' ' << pl.v2 << '\n';
    }
  }
}

void output_data(int case_num, const Board& b, const State& state)
{
  // 標準出力
  if (exec_mode == 0) {
    write_state(std::cout, b, state);
    return;
  }

  // ファイル出力
  std::ostringstream oss;
  oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ofstream ofs(oss.str());
  if (ofs) {
    write_state(ofs, b, state);
  }
}

vector<int> GetNearOrder(int i, const State& s, const Board& b)
{
  vector<pair<double, int>> distances;
  for (int j = 0; j < s.places.size(); j++) {
    if (j == i) {
      continue;
    }
    double dist = euclidean_distance(b.x[i], b.y[i], b.x[j], b.y[j]);
    distances.push_back({ dist, j });
  }
  sort(distances.begin(), distances.end());
  vector<int> near_order;
  for (const auto& p : distances) {
    near_order.push_back(p.second);
  }
  return near_order;
}

/*
 * Hungarian (Kuhn–Munkres) algorithm for ASSIGNMENT PROBLEM
 * --------------------------------------------------------
 * 与えられた n×n の「利益」行列 profit[i][j] (double) について、
 * 各行・各列から 1 要素ずつ選択し、総利益を最大化する。
 *
 * 返り値:
 *   pair<double, vector<int>> (最大利益, 行 i に割り当てた列 idx[i])
 *   割当が存在しない場合 idx[i] = -1
 *
 * 計算量:  O(n^3)   (n ≤ 20 なら μs〜ms オーダ)
 *
 * 実装方針:
 *   最大化 → 最小化へ変換:  cost[i][j] = maxP - profit[i][j]
 *   ※ maxP は各行の最大値を取って +α した十分大きな定数で OK
 *   double の比較誤差を避けるため EPS を使う。
 */
pair<double, vector<int>> hungarian_max(const vector<vector<double>>& profit)
{
  const double INF = 1e100;
  const double EPS = 1e-12;
  int n = (int)profit.size();
  vector<double> u(n + 1), v(n + 1); // ポテンシャル
  vector<int> p(n + 1), way(n + 1);  // p: col→row, way: col→prev col
  vector<int> assignment(n, -1);     // row i に割り当てた col

  /* --- 最大化 → 最小化 (cost) --- */
  double maxP = -INF;
  for (auto& row : profit) {
    for (double x : row) {
      maxP = max(maxP, x);
    }
  }
  vector<vector<double>> cost(n, vector<double>(n));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      cost[i][j] = maxP - profit[i][j];
    }
  }

  /* --- 本体 --- */
  for (int i = 1; i <= n; ++i) { // 行を 1-origin で管理
    p[0] = i;
    int j0 = 0; // 現在の空列
    vector<double> minv(n + 1, INF);
    vector<char> used(n + 1, false);
    do {
      used[j0] = true;
      int i0 = p[j0], j1 = 0;
      double delta = INF;
      for (int j = 1; j <= n; ++j) {
        if (!used[j]) {
          double cur = cost[i0 - 1][j - 1] - u[i0] - v[j];
          if (cur < minv[j]) {
            minv[j] = cur;
            way[j] = j0;
          }
          if (minv[j] < delta) {
            delta = minv[j];
            j1 = j;
          }
        }
      }
      for (int j = 0; j <= n; ++j) {
        if (used[j]) {
          u[p[j]] += delta;
          v[j] -= delta;
        }
        else {
          minv[j] -= delta;
        }
      }
      j0 = j1;
    } while (p[j0] != 0);

    // 増加パスの逆辿り
    do {
      int j1 = way[j0];
      p[j0] = p[j1];
      j0 = j1;
    } while (j0);
  }

  // 結果の取り出し
  double maxSum = 0.0;
  for (int j = 1; j <= n; ++j) {
    int i = p[j];
    if (i) {
      assignment[i - 1] = j - 1;
      maxSum += profit[i - 1][j - 1];
    }
  }
  return { maxSum, assignment };
}

vector<double> dp;
vector<vector<double>> dp2;
int calculate_score(const Board& b, State& s, const DirectedAcyclicGraph& dag)
{
  if (dp.size() != s.places.size()) {
    dp.assign(s.places.size(), 0.0);
  }
  if (dp2.size() != b.n) {
    dp2.assign(b.n, vector<double>(b.n, 0.0));
  }
  for (int i = 0; i < b.n; i++) {
    for (int j = 0; j < b.n; j++) {
      dp2[i][j] = 0.0;
    }
  }

  const auto& topo = dag.topo_;

  for (int i = 0; i < b.n; i++) {
    // 初期化
    for (auto v : topo) {
      dp[v] = 0.0;
    }
    for (int j = 0; j < b.n; j++) {
      dp[j] = 0.0;
    }
    dp[b.n + b.m] = 1.0;

    // 搬入口からスタート
    for (int v : topo) {
      if (s.places[v].k == -1) {
        continue;
      }
      if (s.places[v].v1 == -1 || s.places[v].v2 == -1) {
        continue;
      }
      int kk = s.places[v].k;
      dp[s.places[v].v1] += dp[v] * b.p[kk][i];
      dp[s.places[v].v2] += dp[v] * (1.0 - b.p[kk][i]);
    }

    for (int j = 0; j < b.n; j++) {
      dp2[i][j] = dp[j];
    }
  }

  auto hangarian_result = hungarian_max(dp2);

  for (int i = 0; i < b.n; i++) {
    for (int j = 0; j < b.n; j++) {
      if (hangarian_result.second[j] == i) {
        s.d[i] = j;
        break;
      }
    }
  }

  double missing = b.n - hangarian_result.first;
  int res = round(1e9 / b.n * missing);
  if (res == 0) {
    res = 1; // 最低スコアは1
  }
  return res;
}

void sample_method(const Board& b, State& s, DirectedAcyclicGraph& dag)
{
  s.d.resize(b.n);
  for (int i = 0; i < b.n; i++) {
    s.d[i] = i;
  }

  // s決定
  {
    vector<int> near_order = GetNearOrder(b.n + b.m, s, b);
    s.s = near_order[0]; // 最も近い場所をsに設定
  }

  s.places[s.s].k = 0;
  s.places[s.s].v1 = 0;
  s.places[s.s].v2 = 1;

  dag.RecreateG();
  dag.RecreateTopologicalSort();

  s.score = calculate_score(b, s, dag);
}

void initialize_board(Board& b, State& s)
{
  b.nearest.clear();
  for (int i = 0; i < b.n + b.m + 1; i++) {
    b.nearest.push_back(GetNearOrder(i, s, b));
  }
}

void initialize_State_d(Board& b, State& s)
{
  s.d.resize(b.n);
  for (int i = 0; i < b.n; i++) {
    s.d[i] = i;
  }
}

void initialize_State_s(Board& b, State& s)
{
  for (auto v : b.nearest[b.n + b.m]) {
    if (v >= b.n) {
      s.s = v; // 最も近い場所をsに設定
      s.places[b.n + b.m].k = 0;
      s.places[b.n + b.m].v1 = v;
      s.places[b.n + b.m].v2 = v;
      break;
    }
  }
}

void initialize_State_k(Board& b, State& s)
{
  for (int i = b.n; i < b.n + b.m; i++) {
    s.places[i].k = rand_xorshift() % b.k;
    s.places[i].v1 = -1;
    s.places[i].v2 = -1;
  }
}

void initialize_State_k_and_DAG_pylamid(Board& b, State& s, DirectedAcyclicGraph& dag)
{
  for (int i = b.n; i < b.n + b.m; i++) {
    s.places[i].k = rand_xorshift() % b.k;
    s.places[i].v1 = -1;
    s.places[i].v2 = -1;
  }

  // ピラミッド型に配置
  int cnt = b.n;
  s.s = cnt;
  s.places[b.n + b.m].k = 0;
  s.places[b.n + b.m].v1 = s.s;
  s.places[b.n + b.m].v2 = s.s;
  cnt++;

  vector<int> row;
  row.push_back(s.s);
  while (cnt + row.size() * 4 + 1 <= b.n + b.m) {
    vector<int> next_row;
    for (int i = 0; i < row.size(); i++) {
      int current = row[i];
      if (i == 0) {
        s.places[current].k = 0;
        s.places[current].v1 = cnt;
        next_row.push_back(cnt);
        cnt++;
        s.places[current].v2 = cnt;
        next_row.push_back(cnt);
        cnt++;
      }
      else {
        s.places[current].k = 0;
        s.places[current].v1 = next_row.back();
        s.places[current].v2 = cnt;
        next_row.push_back(cnt);
        cnt++;
      }
    }
    row = next_row;
  }

  // 最後の行の処理
  int cnt2 = 0;
  for (int i = 0; i < row.size(); i++) {
    int current = row[i];
    s.places[current].k = 0;
    s.places[current].v1 = cnt2 % b.n;
    cnt2++;
    s.places[current].v2 = cnt2 % b.n;
    cnt2++;
  }

  dag.RecreateG();
  dag.RecreateTopologicalSort();
}

bool can_use(int current, int next, const Board& b, const State& s, const vector<P>& edges, const vector<bool>& used)
{
  if (next == b.n + b.m) {
    return false; // INLETはスキップ
  }

  if (b.n <= next && next < b.n + b.m && used[next]) {
    return false;
  }

  int ok = 1;
  for (P edge : edges) {
    if (edge.first == current || edge.second == current) {
      continue;
    }
    if (edge.first == next || edge.second == next) {
      continue;
    }
    //continue;
    if (segments_intersect(b.x[edge.first], b.y[edge.first], b.x[edge.second], b.y[edge.second], b.x[current], b.y[current], b.x[next], b.y[next])) {
      return false;
    }
  }
  return true;
}

void initialize_DAG(Board& b, State& s, DirectedAcyclicGraph& dag)
{
  vector<int> visited(b.n + b.m + 1, 0);
  vector<bool> used(b.n + b.m + 1, false);
  vector<P> edges;
  edges.push_back({ b.n + b.m, s.s }); // sとINLETを結ぶエッジを追加

  vector<vector<int>> parents(b.n + b.m + 1);
  queue<P> q;
  q.push({ s.s, 1 });
  while (!q.empty()) {
    int current = q.front().first;
    int depth = q.front().second;
    q.pop();

    if (used[current]) {
      continue;
    }
    used[current] = true;

    vector<int> nexts;
    //for (int i = 0; i < 10; i++) {
    //  if (nexts.size() >= 2) {
    //    break; // 2つの近い頂点を見つけたら終了
    //  }
    //  int next = b.nearest[current][i];
    //  bool skip = false;
    //  for (auto found : nexts) {
    //    if (found == next) {
    //      skip = true; // 既に見つけた頂点はスキップ
    //    }
    //  }
    //  if (skip) {
    //    continue;
    //  }
    //  if ((b.n <= next && next < b.n + b.m && visited[next] >= depth + 1)) {
    //    if (can_use(current, next, b, s, edges, used)) {
    //      nexts.push_back(next);
    //    }
    //  }
    //}
    //for (int i = 0; i < 10; i++) {
    //  if (nexts.size() >= 2) {
    //    break; // 2つの近い頂点を見つけたら終了
    //  }
    //  int next = b.nearest[current][i];
    //  bool skip = false;
    //  for (auto found : nexts) {
    //    if (found == next) {
    //      skip = true; // 既に見つけた頂点はスキップ
    //    }
    //  }
    //  if (skip) {
    //    continue;
    //  }
    //  if (next < b.n) {
    //    if (can_use(current, next, b, s, edges, used)) {
    //      nexts.push_back(next);
    //    }
    //  }
    //}
    for (int i = 0; i < b.nearest[current].size(); i++) {
      if (nexts.size() >= 2) {
        break; // 2つの近い頂点を見つけたら終了
      }
      int next = b.nearest[current][i];
      bool skip = false;
      for (auto found : nexts) {
        if (found == next) {
          skip = true; // 既に見つけた頂点はスキップ
        }
      }
      if (skip) {
        continue;
      }
      if (next < b.n && i < 3 && depth < 5) {
        continue;
      }
      if (can_use(current, next, b, s, edges, used)) {
        nexts.push_back(next);
      }
    }

    for (auto next : nexts) {
      //for (int i = 0; i < 10; i++) {
      //  int next = b.nearest[current][i];
      //if (!can_use(current, next, b, s, edges, used)) {
      //  continue;
      //}
      if (s.places[current].v1 == -1) {
        s.places[current].v1 = next;
      }
      else {
        s.places[current].v2 = next;
      }
      edges.push_back({ current, next });
      parents[next].push_back(current);

      if (b.n <= next && next < b.n + b.m && visited[next] == 0) {
        q.push({ next,depth + 1 });
        visited[next] = depth + 1;
      }

      if (s.places[current].v2 != -1) {
        break;
      }
    }

    if (s.places[current].v2 == -1) {
      s.places[current].v2 = s.places[current].v1;
    }
  }

  queue<int> q2;
  for (int i = b.n; i < b.n + b.m; i++) {
    if (s.places[i].v1 == -1) {
      q2.push(i);
    }
  }
  while (!q2.empty()) {
    int current = q2.front();
    q2.pop();
    for (int parent : parents[current]) {
      if (s.places[parent].v1 == current && s.places[parent].v2 == current) {
        s.places[parent].v1 = -1;
        s.places[parent].v2 = -1;
        q2.push(parent);
      }
      else if (s.places[parent].v1 == current) {
        s.places[parent].v1 = s.places[parent].v2;
      }
      else {
        s.places[parent].v2 = s.places[parent].v1;
      }
    }
  }

  dag.RecreateG();
  dag.RecreateTopologicalSort();
}

void initialize(Board& b, State& s, DirectedAcyclicGraph& dag)
{
  initialize_board(b, s);
  initialize_State_d(b, s);
  initialize_State_s(b, s);
  if (false) {
    initialize_State_k(b, s);
    initialize_DAG(b, s, dag);
  }
  else {
    initialize_State_k_and_DAG_pylamid(b, s, dag);
  }

  s.score = calculate_score(b, s, dag);
}

// 特定の頂点に到達する各ごみ種類の確率分布を計算
vector<double> calculate_prob_dist_at_vertex(const Board& b, const State& s, int vertex, const vector<int>& topo)
{
  vector<double> result(b.n, 0.0);

  // 各ごみ種類について確率を計算
  for (int i = 0; i < b.n; i++) {
    vector<double> dp(s.places.size(), 0.0);
    dp[s.s] = 1.0;

    for (int v : topo) {
      if (v == vertex) break; // 目的の頂点に到達
      if (s.places[v].k == -1) continue;
      if (s.places[v].v1 == -1 || s.places[v].v2 == -1) continue;

      int kk = s.places[v].k;
      dp[s.places[v].v1] += dp[v] * b.p[kk][i];
      dp[s.places[v].v2] += dp[v] * (1.0 - b.p[kk][i]);
    }

    result[i] = dp[vertex];
  }

  return result;
}

// 評価関数を使って最適な分別器を選択
int select_best_separator(const Board& b, const State& s, int vertex,
  const vector<double>& prob_dist,
  const SeparatorEvaluator* evaluator)
{
  int best_k = 0;
  double best_score = -1e9;

  for (int k = 0; k < b.k; k++) {
    double score = evaluator->evaluate(b, k, prob_dist);
    if (score > best_score) {
      best_score = score;
      best_k = k;
    }
  }

  return best_k;
}

bool try_update_state(State& s, State& best_state, const Board& b, DirectedAcyclicGraph& dag, double temp)
{
  // スコアを計算
  int new_score = calculate_score(b, s, dag);

  double diff_score = (s.score - new_score) * 123.6;
  double prob = exp(diff_score / temp);
  if (prob > rand_01()) {
    if (diff_score > 0) {
      s.score = new_score;
      if (s.score < best_state.score) {
        //cerr << "New score: " << s.score << endl;
        dag.RecreateG();
        dag.RecreateTopologicalSort();
        //cerr << s.score << ' ' << calculate_score(b, s, dag) << endl;
        best_state = s;
      }
      return true;
    }
  }
  return false;
}


void method1(const Board& b, State& s, DirectedAcyclicGraph& dag)
{
  for (int i = 0; i < b.n + b.m + 1; i++) {
    s.places[i].keep_v1 = s.places[i].v1;
    s.places[i].keep_v2 = s.places[i].v2;
  }

  dag.RecreateG();
  dag.RecreateTopologicalSort();

  State best_state = s;

  vector<int> use_places;
  {
    auto topo = dag.topo_;
    for (auto v : topo) {
      if (b.n <= v && v < b.n + b.m) {
        use_places.push_back(v);
      }
    }
  }

  vector<int> goal_places;
  {
    auto topo = dag.topo_;
    for (auto v : topo) {
      if (b.n <= v && v < b.n + b.m && s.places[v].v1 < b.n) {
        goal_places.push_back(v);
      }
    }
  }

  // 初期解を評価関数を使って生成
  std::vector<std::unique_ptr<SeparatorEvaluator>> evaluators;
  evaluators.emplace_back(std::make_unique<EntropyReductionEvaluator>());
  evaluators.emplace_back(std::make_unique<VarianceEvaluator>());

  for (auto& evaluator : evaluators) {
    auto topo = dag.topo_;
    for (auto v : topo) {
      if (b.n <= v && v < b.n + b.m) {
        auto prob_dist = calculate_prob_dist_at_vertex(b, s, v, topo);
        s.places[v].k = select_best_separator(b, s, v, prob_dist, evaluator.get());
      }
    }
    s.score = calculate_score(b, s, dag);
    if (s.score < best_state.score) {
      //cerr << "New score with " << evaluator->name() << ": " << s.score << endl;
      best_state = s;
    }
  }

  // その後、ランダムな初期解も試す
  for (int _ = 0; _ < 10; _++) {
    for (auto v : use_places) {
      s.places[v].k = rand_xorshift() % b.k;
    }
    s.score = calculate_score(b, s, dag);
    if (s.score < best_state.score) {
      //cerr << "New random score: " << s.score << endl;
      best_state = s;
    }
  }

  s = best_state; // ベスト状態を適用
  //std::cerr << "Best initial score: " << s.score << endl;
  //std::cerr << "Starting hill climbing method..." << timer.get_elapsed_time() << " seconds" << endl;

  double now_time = timer.get_elapsed_time();
  const double START_TEMP = 1e1;
  const double END_TEMP = 1e-5;

  int loop = 0;
  while (true) {
    loop++;

    if (loop % 100 == 0) {
      now_time = timer.get_elapsed_time();
      if (now_time > TIME_LIMIT) {
        break;
      }
    }

    if (rand_xorshift() % 12345 == 0) {
      s = best_state;
      dag.RecreateG();
      dag.RecreateTopologicalSort();
    }

    double progress_ratio = now_time / TIME_LIMIT;
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progress_ratio;

    int ra = rand_xorshift() % 250;
    if (ra < 95) {
      int v = rand_xorshift() % use_places.size();
      while (s.places[use_places[v]].v1 == s.places[use_places[v]].v2) {
        v = rand_xorshift() % use_places.size();
      }
      int place_id = use_places[v];
      int keep_k = s.places[place_id].k;
      s.places[place_id].k = rand_xorshift() % b.k; // ランダムにkを変更
      while (s.places[place_id].k == keep_k) {
        s.places[place_id].k = rand_xorshift() % b.k;
      }

      bool is_update = try_update_state(s, best_state, b, dag, temp);

      if (!is_update) {
        // 元に戻す
        s.places[place_id].k = keep_k; // 元のkに戻す
      }
    }
    else if (ra < 100) {
      int v = rand_xorshift() % use_places.size();
      while (s.places[use_places[v]].v1 == s.places[use_places[v]].v2) {
        v = rand_xorshift() % use_places.size();
      }
      int place_id = use_places[v];
      swap(s.places[place_id].v1, s.places[place_id].v2); // v1とv2を入れ替える

      bool is_update = try_update_state(s, best_state, b, dag, temp);

      if (!is_update) {
        // 元に戻す
        swap(s.places[place_id].v1, s.places[place_id].v2); // v1とv2を元に戻す
      }
    }
    else if (ra < 200) {
      //State before_state = s; // 元の状態を保存
      ////cerr << calculate_score(b, s, dag) << endl;
      //auto keep_dp2 = dp2; // dp2を保存
      //auto before_g = dag.g_;
      //auto before_topo = dag.topo_;
      //auto before_topoPos = dag.topoPos_;
      //auto before_places = s.places;
      //auto before_s_d = s.d;
      int v = rand_xorshift() % use_places.size();
      while (s.places[use_places[v]].keep_v1 == s.places[use_places[v]].keep_v2) {
        v = rand_xorshift() % use_places.size();
      }
      int place_id = use_places[v];
      int keep_v1 = s.places[place_id].v1;
      int keep_v2 = s.places[place_id].v2;

      if (s.places[place_id].v1 == s.places[place_id].v2) {
        //cerr << "a " << s.places[place_id].v1 << ' ' << s.places[place_id].v2 << endl;
        s.places[place_id].v1 = s.places[place_id].keep_v1; // v1を元に戻す
        s.places[place_id].v2 = s.places[place_id].keep_v2; // v2を元に戻す
        //cerr << "b " << s.places[place_id].v1 << ' ' << s.places[place_id].v2 << endl;
        if (rand_xorshift() % 2 == 0) {
          swap(s.places[place_id].v1, s.places[place_id].v2); // v1とv2を入れ替える
        }
      }
      else {
        if (rand_xorshift() % 2 == 0) {
          s.places[place_id].v1 = s.places[place_id].v2; // v1をv2に置き換える
        }
        else {
          s.places[place_id].v2 = s.places[place_id].v1; // v2をv1に置き換える
        }
      }
      dag.RecreateG();
      dag.RecreateTopologicalSort();

      bool is_update = try_update_state(s, best_state, b, dag, temp);

      if (!is_update) {
        // 元に戻す
        //s = before_state; // 元の状態に戻す
        s.places[place_id].v1 = keep_v1; // 元のv1に戻す
        s.places[place_id].v2 = keep_v2; // 元のv2に戻す
        //s = best_state;
        //cerr << "c " << s.places[place_id].v1 << ' ' << s.places[place_id].v2 << endl;
        dag.RecreateG();
        dag.RecreateTopologicalSort();
        //if (dag.g_ != before_g) {
        //  cerr << "Graph changed!" << endl;
        //}
        //if (dag.topo_ != before_topo) {
        //  cerr << "Topological sort changed!" << endl;
        //}
        //if (dag.topoPos_ != before_topoPos) {
        //  cerr << "Topological position changed!" << endl;
        //}
        //if (dag.state_.places != before_places) {
        //  cerr << "Places changed!" << endl;
        //}
        //if (dag.state_.d != before_s_d) {
        //  cerr << "s.d changed!" << endl;
        //}
        //cerr << calculate_score(b, s, dag) << endl;
        //calculate_score(b, s, dag);
        //if (dp2 != keep_dp2) {
        //  cerr << "dp2 changed!" << endl;
        //  for (int i = 0; i < b.n; i++) {
        //    for (int j = 0; j < b.n; j++) {
        //      cerr << keep_dp2[i][j] << ' ';
        //    }
        //    cerr << endl;
        //  }
        //  cerr << endl;
        //  for (int i = 0; i < b.n; i++) {
        //    for (int j = 0; j < b.n; j++) {
        //      cerr << dp2[i][j] << ' ';
        //    }
        //    cerr << endl;
        //  }
        //}
        //cerr << s.score << ' ' << calculate_score(b, s, dag);
        //for (int i = 0; i < b.n; i++) {
        //  cerr << ' ' << s.d[i];
        //}
        //cerr << endl;
        //if (s.score != calculate_score(b, s, dag)) {
        //  cerr << "Score mismatch!" << endl;
        //  s = best_state; // 元の状態に戻す
        //  dag.RecreateG();
        //  dag.RecreateTopologicalSort();
        //  cerr << s.score << ' ' << calculate_score(b, s, dag) << endl;
        //}
      }
      else {
        //cerr << "Updated state with v1/v2 swap: " << s.score << ' ' << calculate_score(b, s, dag) << endl;
      }
    }
    else if (ra < 250) {
      int v = rand_xorshift() % goal_places.size();
      int place_id = goal_places[v];
      int keep_v1 = s.places[place_id].v1;
      int keep_v2 = s.places[place_id].v2;
      if (rand_xorshift() % 2 == 0) {
        s.places[place_id].v1 = rand_xorshift() % b.n;
        while (s.places[place_id].v1 == keep_v1) {
          s.places[place_id].v1 = rand_xorshift() % b.n;
        }
      }
      else {
        s.places[place_id].v2 = rand_xorshift() % b.n;
        while (s.places[place_id].v2 == keep_v2) {
          s.places[place_id].v2 = rand_xorshift() % b.n;
        }
      }
      dag.RecreateG();
      dag.RecreateTopologicalSort();

      bool is_update = try_update_state(s, best_state, b, dag, temp);

      if (!is_update) {
        // 元に戻す
        s.places[place_id].v1 = keep_v1; // 元のv1に戻す
        s.places[place_id].v2 = keep_v2; // 元のv2に戻す
        dag.RecreateG();
        dag.RecreateTopologicalSort();
      }
    }
  }

  s = best_state;

  std::cerr << "loop = " << loop << ", score = " << s.score << endl;
}

double graph_distance_sum(const Board& b, const State& s, const vector<int>& indices)
{
  const int N = static_cast<int>(s.places.size());

  vector<bool> visited(N, false);
  queue<int> q;
  q.push(indices[b.n + b.m]);
  visited[indices[b.n + b.m]] = true;
  double total_distance = 0.0;

  while (!q.empty()) {
    int current = q.front();
    q.pop();
    if (s.places[current].v1 == -1 || s.places[current].v2 == -1) {
      continue;
    }

    int v1 = indices[s.places[current].v1];
    int v2 = indices[s.places[current].v2];

    double dist1 = euclidean_distance(b.x[current], b.y[current], b.x[v1], b.y[v1]);
    total_distance += dist1;
    if (!visited[v1]) {
      q.push(v1);
      visited[v1] = true;
    }

    if (v1 != v2) {
      double dist2 = euclidean_distance(b.x[current], b.y[current], b.x[v2], b.y[v2]);
      total_distance += dist2;
      if (!visited[v2]) {
        q.push(v2);
        visited[v2] = true;
      }
    }
  }
  return total_distance;
}

void optimize_vertex(const Board& b, State& s, DirectedAcyclicGraph& dag)
{
  const int N = static_cast<int>(s.places.size());
  vector<int> indices;
  for (int i = 0; i < N; i++) {
    indices.push_back(i);
  }

  double distance_sum = graph_distance_sum(b, s, indices);
  cerr << "Initial distance_sum = " << distance_sum << endl;

  int loop = 0;
  while (true) {
    loop++;

    if (loop % 100 == 0) {
      if (timer.get_elapsed_time() > TIME_LIMIT + 1.0) {
        break;
      }
    }

    int ra = rand_xorshift() % 100;
    int swa1 = -1, swa2 = -1;
    if (ra < 90) {
      swa1 = rand_xorshift() % b.m + b.n;
      swa2 = rand_xorshift() % b.m + b.n;
      while (swa1 == swa2) {
        swa2 = rand_xorshift() % b.m + b.n;
      }
    }
    else if (ra < 100) {
      swa1 = rand_xorshift() % b.n;
      swa2 = rand_xorshift() % b.n;
      while (swa1 == swa2) {
        swa2 = rand_xorshift() % b.n;
      }
    }

    swap(indices[swa1], indices[swa2]); // インデックスを入れ替える

    double new_distance_sum = graph_distance_sum(b, s, indices);

    if (new_distance_sum < distance_sum) {
      distance_sum = new_distance_sum;
    }
    else {
      swap(indices[swa1], indices[swa2]); // 元に戻す
    }
  }

  cerr << "optimize_vertex loop = " << loop << ", distance_sum = " << distance_sum << endl;

  State new_s = s;
  new_s.s = indices[s.s];
  for (int i = 0; i < N; i++) {
    int new_i = indices[i];
    new_s.places[new_i] = s.places[i];
    if (new_s.places[new_i].v1 != -1) {
      new_s.places[new_i].v1 = indices[new_s.places[new_i].v1];
    }
    if (new_s.places[new_i].v2 != -1) {
      new_s.places[new_i].v2 = indices[new_s.places[new_i].v2];
    }
    if (new_s.places[new_i].keep_v1 != -1) {
      new_s.places[new_i].keep_v1 = indices[new_s.places[new_i].keep_v1];
    }
    if (new_s.places[new_i].keep_v2 != -1) {
      new_s.places[new_i].keep_v2 = indices[new_s.places[new_i].keep_v2];
    }
  }

  s = new_s;
}

int solve_case(int case_num)
{
  timer.start();

  GameData gameData = input_data(case_num);
  State& state = gameData.state;
  Board& board = gameData.board;
  DirectedAcyclicGraph dag(state, board);

  //sample_method(state);

  initialize(board, state, dag);
  method1(board, state, dag);

  optimize_vertex(board, state, dag);

  dag.RecreateG();
  dag.RecreateTopologicalSort();
  state.score = calculate_score(board, state, dag);

  output_data(case_num, board, state);

  return state.score;
}

int main()
{
  exec_mode = 2;

  if (exec_mode == 0) {
    solve_case(0);
  }
  else if (exec_mode <= 2) {
    ll sum_score = 0;
    for (int i = 0; i < 10; i++) {
      ll score = solve_case(i);
      sum_score += score;
      if (exec_mode == 1) {
        std::cerr << score << endl;
      }
      else {
        std::cerr << "case = " << setw(2) << i << ", "
          << "score = " << setw(4) << score << ", "
          << "sum = " << setw(5) << sum_score << ", "
          << "time = " << setw(5) << timer.get_elapsed_time() << ", "
          << endl;
      }
    }
  }

  return 0;
}
