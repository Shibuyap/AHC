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

enum class PlaceType
{
  PROCESSOR = 0,
  SORTER = 1,
  INLET = 2,
};

class Place
{
public:
  PlaceType type;
  int id;
  int x, y;
  int k;
  int v1, v2;
};

class State
{
public:
  // 入力データ
  int n, m, k;
  vector<Place> places;
  vector<vector<double>> p; // k * nの行列

  // 出力データ
  vector<int> d;
  int s;

  // その他
  int score;
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

vector<int> topological_sort(const State& s)
{
  int N = (int)s.places.size();
  vector<int> indeg(N, 0);
  for (int i = s.n; i < s.n + s.m; ++i) {
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

// 入力データの読み込み
State input_data(int case_num)
{
  State state;

  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    cin >> state.n >> state.m >> state.k;

    state.places.resize(state.n + state.m + 1);
    for (int i = 0; i < state.n; i++) {
      int x, y;
      cin >> x >> y;
      state.places[i] = { PlaceType::PROCESSOR, i, x, y, -1, 0, 0 };
    }
    for (int i = 0; i < state.m; i++) {
      int x, y;
      cin >> x >> y;
      state.places[state.n + i] = { PlaceType::SORTER, i, x, y, -1, 0, 0 };
    }
    state.places[state.n + state.m] = { PlaceType::INLET, 0, 0, 5000, -1, 0, 0 };

    state.p.assign(state.k, vector<double>(state.n));
    for (int i = 0; i < state.k; i++) {
      for (int j = 0; j < state.n; j++) {
        cin >> state.p[i][j];
      }
    }
  }
  else {
    // ファイル入力
    ifs >> state.n >> state.m >> state.k;

    state.places.resize(state.n + state.m + 1);
    for (int i = 0; i < state.n; i++) {
      int x, y;
      ifs >> x >> y;
      state.places[i] = { PlaceType::PROCESSOR, i, x, y, -1, 0, 0 };
    }
    for (int i = 0; i < state.m; i++) {
      int x, y;
      ifs >> x >> y;
      state.places[state.n + i] = { PlaceType::SORTER, i, x, y, -1, 0, 0 };
    }
    state.places[state.n + state.m] = { PlaceType::INLET, 0, 0, 5000, -1, 0, 0 };

    state.p.assign(state.k, vector<double>(state.n));
    for (int i = 0; i < state.k; i++) {
      for (int j = 0; j < state.n; j++) {
        ifs >> state.p[i][j];
      }
    }

    ifs.close();
  }

  return state;
}

void output_data(int case_num, const State& state)
{
  if (exec_mode == 0) {
    // 標準出力
    for (int i = 0; i < state.n; i++) {
      cout << state.d[i] << " ";
    }
    cout << endl;
    cout << state.s << endl;
    for (int i = state.n; i < state.n + state.m; i++) {
      if (state.places[i].k == -1) {
        cout << -1 << endl;
      }
      else {
        cout << state.places[i].k << " " << state.places[i].v1 << " " << state.places[i].v2 << endl;
      }
    }
  }
  else {
    // ファイル出力
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofstream ofs(oss.str());

    if (ofs.is_open()) {
      for (int i = 0; i < state.n; i++) {
        ofs << state.d[i] << " ";
      }
      ofs << endl;
      ofs << state.s << endl;
      for (int i = state.n; i < state.n + state.m; i++) {
        if (state.places[i].k == -1) {
          ofs << -1 << endl;
        }
        else {
          ofs << state.places[i].k << " " << state.places[i].v1 << " " << state.places[i].v2 << endl;
        }
      }

      ofs.close();
    }
  }
}

vector<int> GetNearOrder(int i, const State& s)
{
  vector<pair<double, int>> distances;
  for (int j = 0; j < s.n + s.m; j++) {
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

int calculate_score(const State& s)
{
  double missing = 0.0;
  auto topo = topological_sort(s);

  for (int i = 0; i < s.n; i++) {
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
      dp[s.places[v].v1] += dp[v] * s.p[kk][i];
      dp[s.places[v].v2] += dp[v] * (1.0 - s.p[kk][i]);
    }

    missing += 1.0;
    for (int j = 0; j < s.n; j++) {
      if (s.d[j] == i) {
        missing -= dp[j];
      }
    }
  }

  int res = round(1e9 / s.n * missing);
  return res;
}

void initialize(State& s)
{
  s.d.resize(s.n);
  for (int i = 0; i < s.n; i++) {
    s.d[i] = i;
  }

  vector<P> edges;

  // s決定
  {
    vector<int> near_order = GetNearOrder(s.n + s.m, s);
    s.s = near_order[0]; // 最も近い場所をsに設定
    edges.push_back({ s.n + s.m, s.s }); // sとINLETを結ぶエッジを追加
  }

  vector<vector<int>> parents(s.n + s.m + 1);
  queue<int> q;
  if (s.s >= s.n) {
    q.push(s.s);
    //parents[s].push_back(n + m);
  }

  while (!q.empty()) {
    int current = q.front();
    q.pop();
    auto near_order = GetNearOrder(current, s);
    s.places[current].k = rand_xorshift() % s.k;
    s.places[current].v1 = -1;
    s.places[current].v2 = -1;
    for (int next : near_order) {
      if (next < s.n) {
        int ok = 1;
        for (P edge : edges) {
          if (edge.first == current || edge.second == current) {
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

  for (int i = s.n; i < s.n + s.m; i++) {
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

  s.score = calculate_score(s);
}

// 山登り
void method1(State& s)
{
  State best_state = s;

  vector<int> use_places;
  {
    auto topo = topological_sort(s);
    for (auto v : topo) {
      if (s.n <= v && v < s.n + s.m && s.places[v].k != -1) {
        use_places.push_back(v);
      }
    }
  }

  double now_time = timer.get_elapsed_time();
  const double START_TEMP = 1e9;
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

    int ra = rand_xorshift() % 100;
    if (ra < 50) {
      int n1 = rand_xorshift() % s.n;
      int n2 = rand_xorshift() % s.n;
      while (n1 == n2) {
        n2 = rand_xorshift() % s.n;
      }
      // dを入れ替える
      swap(s.d[n1], s.d[n2]);
      // スコアを計算
      int new_score = calculate_score(s);

      double diff_score = (s.score - new_score) * 123.6;
      double prob = exp(diff_score / temp);
      if (prob > rand_01()) {
        s.score = new_score;
        //cerr << "New score: " << s.score << endl;
        if (s.score < best_state.score) {
          best_state = s; // ベスト状態を更新
        }
      }
      else {
        // 元に戻す
        swap(s.d[n1], s.d[n2]);
      }
    }
    else {
      int v = rand_xorshift() % use_places.size();
      while (s.places[use_places[v]].v1 == s.places[use_places[v]].v2) {
        v = rand_xorshift() % use_places.size();
      }
      int place_id = use_places[v];
      int keep_k = s.places[place_id].k;
      s.places[place_id].k = rand_xorshift() % s.k; // ランダムにkを変更
      while (s.places[place_id].k == keep_k) {
        s.places[place_id].k = rand_xorshift() % s.k;
      }

      // スコアを計算
      int new_score = calculate_score(s);

      double diff_score = (s.score - new_score) * 12345.6;
      double prob = exp(diff_score / temp);
      if (prob > rand_01()) {
        s.score = new_score;
        //cerr << "New score: " << s.score << endl;
        if (s.score < best_state.score) {
          best_state = s; // ベスト状態を更新
        }
      }
      else {
        // 元に戻す
        s.places[place_id].k = keep_k; // 元のkに戻す
      }
    }
  }

  cerr << "loop = " << loop << ", score = " << s.score << endl;
}

int solve_case(int case_num)
{
  timer.start();

  State state = input_data(case_num);

  initialize(state);

  method1(state);

  output_data(case_num, state);

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
