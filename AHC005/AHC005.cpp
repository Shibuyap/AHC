#include <algorithm>
#include <chrono>
#include <climits>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <math.h>
#include <queue>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace std;

typedef long long int ll;
using Point = pair<int, int>;
using Path = vector<int>;

// 定数
namespace Constants
{
  constexpr double TIME_LIMIT = 2.9;
  constexpr double START_TEMP = 1800.0;
  constexpr double END_TEMP = 0.0;
  constexpr int MAX_ITERATIONS = 2000000;
  constexpr int NEIGHBOR_SWAP = 50;
  constexpr int NEIGHBOR_REVERSE = 100;
  constexpr int NEIGHBOR_INSERT = 300;
  constexpr double SCORE_MULTIPLIER = 12345.6;  // 焼きなましの温度計算用スケーリング係数（経験的に良い値）
  constexpr int WALL = -1;
}

// タイマー
namespace Timer
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
namespace Random
{
  uint32_t rand_xorshift()
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

  double rand_01()
  {
    return (rand_xorshift() + 0.5) * (1.0 / UINT_MAX);
  }

  double rand_range(double l, double r)
  {
    return l + (r - l) * rand_01();
  }

  uint32_t rand_range(uint32_t l, uint32_t r)
  {
    return l + rand_xorshift() % (r - l + 1);
  }

  void shuffle_array(int* arr, int n)
  {
    for (int i = n - 1; i >= 0; i--) {
      int j = rand_xorshift() % (i + 1);
      std::swap(arr[i], arr[j]);
    }
  }
}

int mode;

// 座標変換ヘルパー関数
inline int to_index(int x, int y, int n)
{
  return x * n + y;
}

inline Point to_coord(int index, int n)
{
  return { index / n, index % n };
}

inline int get_x(int index, int n)
{
  return index / n;
}

inline int get_y(int index, int n)
{
  return index % n;
}

// 辺を頂点とみなす
class Vertex
{
public:
  vector<int> indices;
};

class DistToVertex
{
public:
  int distance;
  int point_index;

  DistToVertex(int d, int idx) : distance(d), point_index(idx) {}
};

class Board
{
public:
  int n;
  int si, sj;
  vector<vector<int>> grid;
  vector<vector<int>> graph;
  vector<vector<int>> dist;

  vector<Vertex> vertices;
  vector<vector<DistToVertex>> dist_to_vertices; // 各頂点への距離と対応する点のインデックス

  void init_vertices()
  {
    for (int i = 0; i < n * n; i++) {
      int x = get_x(i, n), y = get_y(i, n);
      if (grid[x][y] == Constants::WALL) { continue; }
      if (x == 0 || grid[x - 1][y] == Constants::WALL) {
        Vertex v;
        int cur_x = x;
        while (cur_x < n && grid[cur_x][y] != Constants::WALL) {
          v.indices.push_back(to_index(cur_x, y, n));
          cur_x++;
        }
        if (v.indices.size() > 1) {
          vertices.push_back(v);
        }
      }
      if (y == 0 || grid[x][y - 1] == Constants::WALL) {
        Vertex v;
        int cur_y = y;
        while (cur_y < n && grid[x][cur_y] != Constants::WALL) {
          v.indices.push_back(to_index(x, cur_y, n));
          cur_y++;
        }
        if (v.indices.size() > 1) {
          vertices.push_back(v);
        }
      }
    }
  }

  void init_dist_to_vertices()
  {
    dist_to_vertices.clear();
    dist_to_vertices.resize(n * n);
    for (int i = 0; i < n * n; i++) {
      if (grid[get_x(i, n)][get_y(i, n)] == Constants::WALL) { continue; }
      for (size_t j = 0; j < vertices.size(); j++) {
        const Vertex& v = vertices[j];
        int min_dist = INT_MAX;
        int point_index = -1;
        for (int idx : v.indices) {
          int d = dist[i][idx];
          if (d < min_dist) {
            min_dist = d;
            point_index = idx;
          }
        }
        if (min_dist < INT_MAX) {
          dist_to_vertices[i].push_back(DistToVertex(min_dist, point_index));
        }
      }
    }
  }

  void init_graph()
  {
    // グラフの初期化
    graph.resize(n * n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (grid[i][j] == Constants::WALL) { continue; }
        int idx = to_index(i, j, n);
        // 上下左右の隣接点を追加
        if (i > 0 && grid[i - 1][j] != Constants::WALL) graph[idx].push_back(to_index(i - 1, j, n));
        if (i < n - 1 && grid[i + 1][j] != Constants::WALL) graph[idx].push_back(to_index(i + 1, j, n));
        if (j > 0 && grid[i][j - 1] != Constants::WALL) graph[idx].push_back(to_index(i, j - 1, n));
        if (j < n - 1 && grid[i][j + 1] != Constants::WALL) graph[idx].push_back(to_index(i, j + 1, n));
      }
    }
  }

  inline int get_cost(int idx) const
  {
    int i = get_x(idx, n), j = get_y(idx, n);
    if (grid[i][j] == Constants::WALL) return INT_MAX;
    return grid[i][j];
  }

  // 全点対最短経路を計算する
  void calc_dist()
  {
    dist.resize(n * n, vector<int>(n * n, INT_MAX));
    for (int start = 0; start < n * n; start++) {
      if (grid[get_x(start, n)][get_y(start, n)] == Constants::WALL) { continue; }
      dist[start] = dijkstra(start);
    }
  }

  void init()
  {
    init_vertices();
    init_graph();
    calc_dist();
    init_dist_to_vertices();
  }

  // ダイクストラ法で最短経路を計算
  vector<int> dijkstra(int start) const
  {
    vector<int> d(n * n, INT_MAX);
    d[start] = 0;
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    pq.push({ 0, start });
    while (!pq.empty()) {
      auto [cost, u] = pq.top();
      pq.pop();
      if (cost > d[u]) { continue; }
      for (int v : graph[u]) {
        if (d[v] > d[u] + get_cost(v)) {
          d[v] = d[u] + get_cost(v);
          pq.push({ d[v], v });
        }
      }
    }
    return d;
  }

  // 経路復元付きダイクストラ
  pair<vector<int>, vector<int>> dijkstra_with_path(int start) const
  {
    vector<int> d(n * n, INT_MAX);
    vector<int> prev(n * n, -1);
    d[start] = 0;
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    pq.push({ 0, start });
    while (!pq.empty()) {
      auto [cost, u] = pq.top();
      pq.pop();
      if (cost > d[u]) { continue; }
      for (int v : graph[u]) {
        if (d[v] > d[u] + get_cost(v)) {
          d[v] = d[u] + get_cost(v);
          prev[v] = u;
          pq.push({ d[v], v });
        }
      }
    }
    return { d, prev };
  }

  // パスを復元
  Path restore_path(const vector<int>& prev, int goal) const
  {
    Path path;
    for (int v = goal; v != -1; v = prev[v]) {
      path.push_back(v);
    }
    reverse(path.begin(), path.end());
    return path;
  }

  Path get_path(int start, int goal) const
  {
    auto [d, prev] = dijkstra_with_path(start);
    return restore_path(prev, goal);
  }
};

class Answer
{
public:
  vector<int> vertices;
  vector<int> points;
  int score;

  Answer() : score(0) {}
  void clear()
  {
    points.clear();
    vertices.clear();
    score = 0;
  }

  int recalc_points(const Board& board, int start_index = 0, int end_index = 1001001)
  {
    if (points.size() > vertices.size()) {
      points.clear();
    }
    if (points.size() < vertices.size()) {
      points.resize(vertices.size());
    }
    points[0] = to_index(board.si, board.sj, board.n);
    int cur_index = to_index(board.si, board.sj, board.n);
    if (start_index > 0) {
      cur_index = points[start_index];
    }
    for (int i = start_index + 1; i < vertices.size() - 1; i++) {
      int next_index = board.dist_to_vertices[cur_index][vertices[i]].point_index;
      if (i >= end_index && points[i] == next_index) {
        break;
      }
      points[i] = next_index;
      cur_index = next_index;
    }
    points[points.size() - 1] = to_index(board.si, board.sj, board.n);
    return 0;
  }

  int calc_score(const Board& board)
  {
    int cost = 0;
    for (size_t i = 0; i < points.size() - 1; i++) {
      int start = points[i];
      int goal = points[i + 1];
      if (start < 0 || start >= board.n * board.n || goal < 0 || goal >= board.n * board.n) {
        return -1;
      }
      cost += board.dist[start][goal];
    }
    score = round(1e4 + 1e7 * board.n / cost);
    return score;
  }
};

Board input_data(int cn)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << cn << ".txt";
  ifstream ifs(oss.str());

  Board board;

  if (!ifs.is_open()) {
    // 標準入力
    cin >> board.n >> board.si >> board.sj;
    board.grid.resize(board.n, vector<int>(board.n));
    for (int i = 0; i < board.n; i++) {
      string line;
      cin >> line;
      for (int j = 0; j < board.n; j++) {
        board.grid[i][j] = (line[j] == '#') ? Constants::WALL : line[j] - '0';
      }
    }
  }
  else {
    // ファイル入力
    ifs >> board.n >> board.si >> board.sj;
    board.grid.resize(board.n, vector<int>(board.n));
    for (int i = 0; i < board.n; i++) {
      string line;
      ifs >> line;
      for (int j = 0; j < board.n; j++) {
        board.grid[i][j] = (line[j] == '#') ? Constants::WALL : line[j] - '0';
      }
    }
  }

  board.init();

  return board;
}

// パスから方向文字列に変換
string path_to_directions(const Path& path, int n)
{
  string directions;
  for (size_t i = 0; i < path.size() - 1; i++) {
    int u = path[i], v = path[i + 1];
    int di = get_x(v, n) - get_x(u, n);
    int dj = get_y(v, n) - get_y(u, n);

    if (di == -1 && dj == 0) {
      directions += 'U';
    }
    else if (di == 1 && dj == 0) {
      directions += 'D';
    }
    else if (di == 0 && dj == -1) {
      directions += 'L';
    }
    else if (di == 0 && dj == 1) {
      directions += 'R';
    }
    else {
      cerr << "Invalid move from " << u << " to " << v << endl;
    }
  }
  return directions;
}

// 完全なパスを構築
Path build_complete_path(const Board& board, const Answer& ans)
{
  Path path;
  path.push_back(ans.points[0]);
  for (size_t i = 0; i < ans.points.size() - 1; i++) {
    auto tmp = board.get_path(ans.points[i], ans.points[i + 1]);
    tmp.erase(tmp.begin());
    path.insert(path.end(), tmp.begin(), tmp.end());
  }
  return path;
}

void output_data(int cn, const Board& board, const Answer& ans)
{
  Path path = build_complete_path(board, ans);
  string directions = path_to_directions(path, board.n);

  if (mode == 0) {
    // 標準出力
    cout << directions << endl;
  }
  else {
    // ファイル出力
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << cn << ".txt";
    ofstream ofs(oss.str());

    ofs << directions << endl;

    if (ofs.is_open()) {
      ofs.close();
    }
  }
}

void build_initial_path(const Board& board, Answer& ans)
{
  ans.clear();
  ans.vertices.clear();
  ans.points.push_back(to_index(board.si, board.sj, board.n));
  ans.vertices.push_back(-1);

  int cur_index = board.si * board.n + board.sj;
  for (int i = 0; i < board.vertices.size(); i++) {
    int next_index = board.dist_to_vertices[cur_index][i].point_index;
    ans.points.push_back(next_index);
    ans.vertices.push_back(i);
    cur_index = next_index;
  }
  ans.points.push_back(board.si * board.n + board.sj); // スタート位置に戻る
  ans.vertices.push_back(-1);
}

// 近傍操作の情報を保持する構造体
struct NeighborInfo
{
  int type;
  int ra1, ra2;
  int start_index;
  int end_index;

  NeighborInfo() : type(0), ra1(0), ra2(0), start_index(0), end_index(0) {}
};

// 2点入れ替え
void apply_swap(Answer& ans, int ra1, int ra2)
{
  swap(ans.vertices[ra1], ans.vertices[ra2]);
}

// 区間反転
void apply_reverse(Answer& ans, int ra1, int ra2)
{
  reverse(ans.vertices.begin() + ra1, ans.vertices.begin() + ra2 + 1);
}

// 1点挿入
void apply_insert(Answer& ans, int ra1, int ra2)
{
  if (ra1 < ra2) {
    for (int i = ra1; i < ra2; i++) {
      swap(ans.vertices[i], ans.vertices[i + 1]);
    }
  }
  else {
    for (int i = ra1; i > ra2; i--) {
      swap(ans.vertices[i], ans.vertices[i - 1]);
    }
  }
}

// 近傍操作を生成
NeighborInfo generate_neighbor(Answer& ans)
{
  NeighborInfo info;
  info.type = Random::rand_xorshift() % 300;

  info.ra1 = Random::rand_xorshift() % (ans.points.size() - 2) + 1;
  info.ra2 = Random::rand_xorshift() % (ans.points.size() - 2) + 1;
  while (info.ra1 == info.ra2) {
    info.ra2 = Random::rand_xorshift() % (ans.points.size() - 2) + 1;
  }

  if (info.type < Constants::NEIGHBOR_SWAP) {
    if (info.ra1 > info.ra2) swap(info.ra1, info.ra2);
    apply_swap(ans, info.ra1, info.ra2);
    info.start_index = info.ra1 - 1;
    info.end_index = info.ra2 + 1;
  }
  else if (info.type < Constants::NEIGHBOR_REVERSE) {
    if (info.ra1 > info.ra2) swap(info.ra1, info.ra2);
    apply_reverse(ans, info.ra1, info.ra2);
    info.start_index = info.ra1 - 1;
    info.end_index = info.ra2 + 1;
  }
  else if (info.type < Constants::NEIGHBOR_INSERT) {
    apply_insert(ans, info.ra1, info.ra2);
    info.start_index = min(info.ra1, info.ra2) - 1;
    info.end_index = max(info.ra1, info.ra2) + 1;
  }

  return info;
}

// 近傍操作を元に戻す
void revert_neighbor(Answer& ans, const NeighborInfo& info)
{
  if (info.type < Constants::NEIGHBOR_SWAP) {
    apply_swap(ans, info.ra1, info.ra2);
  }
  else if (info.type < Constants::NEIGHBOR_REVERSE) {
    apply_reverse(ans, info.ra1, info.ra2);
  }
  else if (info.type < Constants::NEIGHBOR_INSERT) {
    if (info.ra1 < info.ra2) {
      for (int i = info.ra2; i > info.ra1; i--) {
        swap(ans.vertices[i], ans.vertices[i - 1]);
      }
    }
    else {
      for (int i = info.ra2; i < info.ra1; i++) {
        swap(ans.vertices[i], ans.vertices[i + 1]);
      }
    }
  }
}

void run_simulated_annealing(double time_limit, const Board& board, Answer& ans)
{
  ans.calc_score(board);
  Answer best_ans = ans;

  double start_time = Timer::get_elapsed_time();
  const double START_TEMP = Constants::START_TEMP;
  const double END_TEMP = Constants::END_TEMP;

  int loop = 0;
  while (loop < Constants::MAX_ITERATIONS) {
    loop++;

    double progress = (double)loop / Constants::MAX_ITERATIONS;
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progress;

    // 近傍解作成
    NeighborInfo nb_info = generate_neighbor(ans);

    // スコア計算
    double cur_score = ans.score;
    ans.recalc_points(board, nb_info.start_index, nb_info.end_index);
    double new_score = ans.calc_score(board);

    // 焼きなましで採用判定
    double diff = (new_score - cur_score) * Constants::SCORE_MULTIPLIER;
    double prob = exp(diff / temp);
    if (prob > Random::rand_01() || Random::rand_xorshift() % 10000 == 0) {
      // 採用
      cur_score = new_score;

      // ベスト更新
      if (cur_score > best_ans.score) {
        best_ans = ans;
      }
    }
    else {
      // 元に戻す
      revert_neighbor(ans, nb_info);
      ans.recalc_points(board, nb_info.start_index, nb_info.end_index);
      ans.score = cur_score;
    }
  }

  if (mode != 3) {
    cerr << loop << endl;
  }

  ans = best_ans;
}


ll solve_case(int cn)
{
  const double TIME_LIMIT = Constants::TIME_LIMIT;
  Timer::start_timer();

  Board board = input_data(cn);

  Answer ans;
  build_initial_path(board, ans);

  Answer initial_ans = ans;
  Answer best_ans = ans;
  const int SET_COUNT = 1;
  double start_time = Timer::get_elapsed_time();
  for (int i = 0; i < SET_COUNT; i++) {
    double time_limit = start_time + (TIME_LIMIT - start_time) / SET_COUNT * (i + 1);
    ans = initial_ans;
    run_simulated_annealing(time_limit, board, ans);
    cerr << "Set " << i + 1 << ": score = " << ans.score << ", time = " << Timer::get_elapsed_time() << endl;
    if (ans.score > best_ans.score) {
      best_ans = ans;
    }
  }
  ans = best_ans;

  output_data(cn, board, ans);

  ll score = 0;
  if (mode != 0) {
    score = ans.calc_score(board);
  }
  return score;
}

int main()
{
  mode = 2;

  if (mode == 0) {
    solve_case(0);
  }
  else if (mode <= 2) {
    ll sum = 0;
    for (int i = 0; i < 5; i++) {
      ll score = solve_case(i);
      sum += score;
      if (mode == 1) {
        cerr << score << endl;
      }
      else {
        cerr << "case = " << setw(2) << i << ", "
          << "score = " << setw(4) << score << ", "
          << "sum = " << setw(5) << sum << ", "
          << "time = " << setw(5) << Timer::get_elapsed_time() << ", "
          << endl;
      }
    }
  }

  return 0;
}
