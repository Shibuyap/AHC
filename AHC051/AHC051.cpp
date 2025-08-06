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
#include <format>
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

class Board
{
public:
  int n, m, k;
  vector<vector<int>> nearest; // 各頂点から別の頂点への近い順
  vector<vector<double>> p; // k * nの行列
};

class Place
{
public:
  int x, y;
  int k;
  int v1, v2;
  int keep_v1, keep_v2;
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

// Place構造体を使った向き計算のラッパー関数
int orientation(const Place& a, const Place& b, const Place& c)
{
  return orientation(a.x, a.y, b.x, b.y, c.x, c.y);
}

// Place構造体を使った線分交差判定のラッパー関数
bool segments_intersect(const Place& p1, const Place& p2, const Place& q1, const Place& q2)
{
  return segments_intersect(p1.x, p1.y, p2.x, p2.y, q1.x, q1.y, q2.x, q2.y);
}

vector<int> topological_sort(const Board& b, const State& s)
{
  int N = (int)s.places.size();
  vector<int> indeg(N, 0);
  for (int i = b.n; i < b.n + b.m; ++i) {
    if (s.places[i].k == -1) {
      continue;
    }
    ++indeg[s.places[i].v1];
    ++indeg[s.places[i].v2];
  }

  queue<int> q;
  q.push(s.s);

  vector<int> order;
  order.reserve(N);
  while (!q.empty()) {
    int v = q.front();
    q.pop();
    order.push_back(v);
    if (s.places[v].k == -1) {
      continue;
    }
    if (--indeg[s.places[v].v1] == 0) {
      q.push(s.places[v].v1);
    }
    if (--indeg[s.places[v].v2] == 0) {
      q.push(s.places[v].v2);
    }
  }
  return order;
}

static GameData read_state(std::istream& is)
{
  State s;
  Board b;
  is >> b.n >> b.m >> b.k;

  s.places.resize(b.n + b.m + 1);

  auto read_places = [&](int cnt, int offset)
    {
      for (int i = 0; i < cnt; ++i) {
        int x, y;
        is >> x >> y;
        s.places[offset + i] = { x, y, -1, 0, 0 };
      }
    };

  read_places(b.n, 0);           // 拠点
  read_places(b.m, b.n);        // 工場
  s.places.back() = { 0, 5000, -1, 0, 0 };   // ダミー

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
  const std::string path = std::format("./in/{:04}.txt", case_num);
  if (std::ifstream ifs{ path }; ifs) {
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
    if (pl.k == -1) {
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
  std::ofstream ofs(std::format("./out/{:04}.txt", case_num));
  if (ofs) {
    write_state(ofs, b, state);
  }
}

vector<int> GetNearOrder(int i, const State& s)
{
  vector<pair<double, int>> distances;
  for (int j = 0; j < s.places.size() - 1; j++) {
    if (j == i) {
      continue;
    }
    double dist = euclidean_distance(s.places[i].x, s.places[i].y, s.places[j].x, s.places[j].y);
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
  for (auto& row : profit)
    for (double x : row)
      maxP = max(maxP, x);
  vector<vector<double>> cost(n, vector<double>(n));
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      cost[i][j] = maxP - profit[i][j];

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
      for (int j = 1; j <= n; ++j)
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
      for (int j = 0; j <= n; ++j) {
        if (used[j]) {
          u[p[j]] += delta;
          v[j] -= delta;
        }
        else
          minv[j] -= delta;
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

int calculate_score(const Board& b, State& s)
{
  auto topo = topological_sort(b, s);

  vector<vector<double>> dp2(b.n, vector<double>(b.n, 0.0));

  for (int i = 0; i < b.n; i++) {
    vector<double> dp(s.places.size(), 0.0);
    dp[s.s] = 1.0;
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
  return res;
}

void initialize(const Board& b, State& s)
{
  s.d.resize(b.n);
  for (int i = 0; i < b.n; i++) {
    s.d[i] = i;
  }

  vector<bool> visited(b.n + b.m, false);
  vector<P> edges;

  // s決定
  {
    vector<int> near_order = GetNearOrder(b.n + b.m, s);
    s.s = near_order[0]; // 最も近い場所をsに設定
    edges.push_back({ b.n + b.m, s.s }); // sとINLETを結ぶエッジを追加
  }

  vector<vector<int>> parents(b.n + b.m + 1);
  queue<int> q;
  if (s.s >= b.n) {
    q.push(s.s);
  }

  while (!q.empty()) {
    int current = q.front();
    q.pop();
    if (visited[current]) {
      continue;
    }
    visited[current] = true;
    auto near_order = GetNearOrder(current, s);
    s.places[current].k = rand_xorshift() % b.k;
    s.places[current].v1 = -1;
    s.places[current].v2 = -1;
    for (int next : near_order) {
      if (next < b.n) {
        int ok = 1;
        for (P edge : edges) {
          if (edge.first == current || edge.second == current) {
            continue;
          }
          if (edge.first == next || edge.second == next) {
            continue;
          }
          if (segments_intersect(s.places[edge.first], s.places[edge.second], s.places[current], s.places[next])) {
            ok = 0;
            break;
          }
        }
        if (!ok) {
          continue;
        }
        if (s.places[current].v1 == -1) {
          s.places[current].v1 = next;
        }
        else {
          s.places[current].v2 = next;
        }
        edges.push_back({ current, next });
        parents[next].push_back(current);
      }
      else {
        if (s.places[next].k != -1) {
          continue;
        }
        int ok = 1;
        for (P edge : edges) {
          if (edge.first == current || edge.second == current) {
            continue;
          }
          if (edge.first == next || edge.second == next) {
            continue;
          }
          if (segments_intersect(s.places[edge.first], s.places[edge.second], s.places[current], s.places[next])) {
            ok = 0;
            break;
          }
        }
        if (!ok) {
          continue;
        }
        if (s.places[current].v1 == -1) {
          s.places[current].v1 = next;
        }
        else {
          s.places[current].v2 = next;
        }
        q.push(next);
        edges.push_back({ current, next });
        parents[next].push_back(current);
      }

      if (s.places[current].v2 != -1) {
        break;
      }
    }

    if (s.places[current].v1 == -1) {
      s.places[current].k = -1;
    }
    else if (s.places[current].v2 == -1) {
      s.places[current].v2 = s.places[current].v1;
    }
  }

  for (int i = b.n; i < b.n + b.m; i++) {
    if (s.places[i].k == -1) {
      q.push(i);
    }
  }
  while (!q.empty()) {
    int current = q.front();
    q.pop();
    for (int parent : parents[current]) {
      if (s.places[parent].v1 == current && s.places[parent].v2 == current) {
        s.places[parent].k = -1;
        q.push(parent);
      }
      else if (s.places[parent].v1 == current) {
        s.places[parent].v1 = s.places[parent].v2;
      }
      else {
        s.places[parent].v2 = s.places[parent].v1;
      }
    }
  }

  s.score = calculate_score(b, s);
}

void sample_method(const Board& b, State& s)
{
  s.d.resize(b.n);
  for (int i = 0; i < b.n; i++) {
    s.d[i] = i;
  }

  // s決定
  {
    vector<int> near_order = GetNearOrder(b.n + b.m, s);
    s.s = near_order[0]; // 最も近い場所をsに設定
  }

  s.places[s.s].k = 0;
  s.places[s.s].v1 = 0;
  s.places[s.s].v2 = 1;

  s.score = calculate_score(b, s);
}

// 山登り
void method1(const Board& b, State& s)
{
  for (int i = 0; i < b.n + b.m; i++) {
    s.places[i].keep_v1 = s.places[i].v1;
    s.places[i].keep_v2 = s.places[i].v2;
  }
  State best_state = s;

  vector<int> use_places;
  {
    auto topo = topological_sort(b, s);
    for (auto v : topo) {
      if (b.n <= v && v < b.n + b.m && s.places[v].k != -1) {
        use_places.push_back(v);
      }
    }
  }

  for (int _ = 0; _ < 10000; _++) {
    for (auto v : use_places) {
      s.places[v].k = rand_xorshift() % b.k;
    }
    s.score = calculate_score(b, s);
    if (s.score < best_state.score) {
      //cerr << "New score: " << s.score << endl;
      best_state = s; // ベスト状態を更新
    }
  }

  cerr << "Starting hill climbing method..." << timer.get_elapsed_time() << " seconds" << endl;


  double now_time = timer.get_elapsed_time();
  const double START_TEMP = 1e6;
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

    double progress_ratio = now_time / TIME_LIMIT;
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progress_ratio;

    int ra = rand_xorshift() % 200;
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

      // スコアを計算
      int new_score = calculate_score(b, s);

      double diff_score = (s.score - new_score) * 123.6;
      double prob = exp(diff_score / temp);
      if (prob > rand_01()) {
        s.score = new_score;
        if (s.score < best_state.score) {
          //cerr << "New score1: " << s.score << endl;
          best_state = s; // ベスト状態を更新
        }
      }
      else {
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

      // スコアを計算
      int new_score = calculate_score(b, s);

      double diff_score = (s.score - new_score) * 123.6;
      double prob = exp(diff_score / temp);
      if (prob > rand_01()) {
        s.score = new_score;
        if (s.score < best_state.score) {
          //cerr << "New score2: " << s.score << endl;
          best_state = s; // ベスト状態を更新
        }
      }
      else {
        // 元に戻す
        swap(s.places[place_id].v1, s.places[place_id].v2); // v1とv2を入れ替える
      }
    }
    else if (ra < 200) {
      int v = rand_xorshift() % use_places.size();
      while (s.places[use_places[v]].keep_v1 == s.places[use_places[v]].keep_v2) {
        v = rand_xorshift() % use_places.size();
      }
      int place_id = use_places[v];
      // v1とv2を入れ替える前に、元の値を保存
      int keep_v1 = s.places[place_id].v1;
      int keep_v2 = s.places[place_id].v2;
      if (s.places[place_id].v1 == s.places[place_id].v2) {
        s.places[place_id].v1 = s.places[place_id].keep_v1; // v1を元に戻す
        s.places[place_id].v2 = s.places[place_id].keep_v2; // v2を元に戻す
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

      // スコアを計算
      int new_score = calculate_score(b, s);

      double diff_score = (s.score - new_score) * 123.6;
      double prob = exp(diff_score / temp);
      if (prob > rand_01()) {
        s.score = new_score;
        if (s.score < best_state.score) {
          //cerr << "New score2: " << s.score << endl;
          best_state = s; // ベスト状態を更新
        }
      }
      else {
        // 元に戻す
        s.places[place_id].v1 = keep_v1; // 元のv1に戻す
        s.places[place_id].v2 = keep_v2; // 元のv2に戻す
      }
    }
  }

  //cerr << s.score << ' ' << best_state.score << endl;
  s = best_state; // ベスト状態を最終的な状態に設定

  cerr << "loop = " << loop << ", score = " << s.score << endl;
}

int solve_case(int case_num)
{
  timer.start();

  GameData gameData = input_data(case_num);
  State& state = gameData.state;
  Board& board = gameData.board;

  //sample_method(state);

  initialize(board, state);
  method1(board, state);

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
    for (int i = 0; i < 100; i++) {
      ll score = solve_case(i);
      sum_score += score;
      if (exec_mode == 1) {
        cerr << score << endl;
      }
      else {
        cerr << "case = " << setw(2) << i << ", "
          << "score = " << setw(4) << score << ", "
          << "sum = " << setw(5) << sum_score << ", "
          << "time = " << setw(5) << timer.get_elapsed_time() << ", "
          << endl;
      }
    }
  }

  return 0;
}
