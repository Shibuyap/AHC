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

// タイマー
namespace {
  std::chrono::steady_clock::time_point start_time_clock;

  void start_timer() {
    start_time_clock = std::chrono::steady_clock::now();
  }

  double get_elapsed_time() {
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time_clock;
    return elapsed.count();
  }
}

// 乱数
namespace {
  static uint32_t rand_xorshift() {
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

  static double rand_01() {
    return (rand_xorshift() + 0.5) * (1.0 / UINT_MAX);
  }

  static double rand_range(double l, double r) {
    return l + (r - l) * rand_01();
  }

  static uint32_t rand_range(uint32_t l, uint32_t r) {
    return l + rand_xorshift() % (r - l + 1); // [l, r]
  }

  void shuffle_array(int* arr, int n) {
    for (int i = n - 1; i >= 0; i--) {
      int j = rand_xorshift() % (i + 1);
      int swa = arr[i];
      arr[i] = arr[j];
      arr[j] = swa;
    }
  }
}

int exec_mode;

// 辺を頂点とみなす
class Vertex {
public:
  vector<int> indices;
};

class DistToVertex {
public:
  int distance;
  int point_index;
};

class Board {
public:
  int n;
  int si, sj;
  vector<vector<int>> grid;
  vector<vector<int>> graph;
  vector<vector<int>> dist;

  vector<Vertex> vertices;
  vector<vector<DistToVertex>> dist_to_vertices; // 各頂点への距離と対応する点のインデックス

  void init_vertices() {
    for (int i = 0; i < n * n; i++) {
      int x = i / n, y = i % n;
      if (grid[x][y] == -1) continue; // 壁は無視
      if (x == 0 || grid[x - 1][y] == -1) {
        Vertex v;
        int cur_x = x;
        while (cur_x < n && grid[cur_x][y] != -1) {
          v.indices.push_back(cur_x * n + y);
          cur_x++;
        }
        if (v.indices.size() > 1) {
          vertices.push_back(v);
        }
      }
      if (y == 0 || grid[x][y - 1] == -1) {
        Vertex v;
        int cur_y = y;
        while (cur_y < n && grid[x][cur_y] != -1) {
          v.indices.push_back(x * n + cur_y);
          cur_y++;
        }
        if (v.indices.size() > 1) {
          vertices.push_back(v);
        }
      }
    }
  }

  void init_dist_to_vertices() {
    dist_to_vertices.clear();
    dist_to_vertices.resize(n * n);
    for (int i = 0; i < n * n; i++) {
      if (grid[i / n][i % n] == -1) continue; // 壁は無視
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
          dist_to_vertices[i].push_back({ min_dist, point_index });
        }
      }
    }
  }

  void init_graph() {
    // グラフの初期化
    graph.resize(n * n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (grid[i][j] == -1) continue; // 壁は無視
        int idx = i * n + j;
        // 上下左右の隣接点を追加
        if (i > 0 && grid[i - 1][j] != -1) graph[idx].push_back((i - 1) * n + j);
        if (i < n - 1 && grid[i + 1][j] != -1) graph[idx].push_back((i + 1) * n + j);
        if (j > 0 && grid[i][j - 1] != -1) graph[idx].push_back(i * n + (j - 1));
        if (j < n - 1 && grid[i][j + 1] != -1) graph[idx].push_back(i * n + (j + 1));
      }
    }
  }

  inline int get_cost(int idx) const {
    int i = idx / n, j = idx % n;
    if (grid[i][j] == -1) return INT_MAX; // 壁は無視
    return grid[i][j]; // コストはグリッドの値
  }

  // 全点対最短経路を計算する
  void calc_dist() {
    // n*n回ダイクストラ
    dist.resize(n * n, vector<int>(n * n, INT_MAX));
    for (int start = 0; start < n * n; start++) {
      if (grid[start / n][start % n] == -1) continue; // 壁は無視
      vector<int> d(n * n, INT_MAX);
      d[start] = 0;
      priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
      pq.push({ 0, start });
      while (!pq.empty()) {
        auto p = pq.top();
        pq.pop();
        int cost = p.first;
        int u = p.second;
        if (cost > d[u]) continue;
        for (int v : graph[u]) {
          if (d[v] > d[u] + get_cost(v)) {
            d[v] = d[u] + get_cost(v);
            pq.push({ d[v], v });
          }
        }
      }
      dist[start] = d;
    }
  }

  void init() {
    init_vertices();
    init_graph();
    calc_dist();
    init_dist_to_vertices();
  }

  vector<int> get_path(int start, int goal) const {
    //cerr << "get_path: start = " << start << ", goal = " << goal << endl;

    // 復元ダイクストラ
    vector<int> prev(n * n, -1);
    vector<int> d(n * n, INT_MAX);
    d[start] = 0;
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    pq.push({ 0, start });
    while (!pq.empty()) {
      auto p = pq.top();
      pq.pop();
      int cost = p.first;
      int u = p.second;
      if (cost > d[u]) continue;
      for (int v : graph[u]) {
        if (d[v] > d[u] + get_cost(v)) {
          d[v] = d[u] + get_cost(v);
          prev[v] = u;
          pq.push({ d[v], v });
        }
      }
    }

    // パスを復元
    vector<int> path;
    for (int v = goal; v != -1; v = prev[v]) {
      path.push_back(v);
    }
    reverse(path.begin(), path.end());

    return path;
  }
};

class Answer {
public:
  vector<int> vertices;
  vector<int> points;
  int score;
  Answer() : score(0) {}
  void clear() {
    points.clear();
    vertices.clear();
    score = 0;
  }

  int recalc_points(const Board& board, int start_index = 0, int end_index = 1001001) {
    if (points.size() > vertices.size()) {
      points.clear();
    }
    if (points.size() < vertices.size()) {
      points.resize(vertices.size());
    }
    points[0] = board.si * board.n + board.sj; // スタート位置を追加
    int cur_index = board.si * board.n + board.sj;
    if (start_index > 0) {
      cur_index = points[start_index]; // スタート位置を指定
    }
    for (int i = start_index + 1; i < vertices.size() - 1; i++) {
      int next_index = board.dist_to_vertices[cur_index][vertices[i]].point_index;
      if (i >= end_index && points[i] == next_index) {
        break;
      }
      points[i] = next_index; // 次の頂点への最短距離を追加
      cur_index = next_index; // 現在のインデックスを更新
    }
    points[points.size() - 1] = board.si * board.n + board.sj; // スタート位置に戻る
    return 0;
  }

  int calc_score(const Board& board) {
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

Board input_data(int case_num) {
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
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
        if (line[j] == '#') {
          board.grid[i][j] = -1;
        }
        else {
          board.grid[i][j] = line[j] - '0';
        }
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
        if (line[j] == '#') {
          board.grid[i][j] = -1;
        }
        else {
          board.grid[i][j] = line[j] - '0';
        }
      }
    }
  }

  board.init();

  return board;
}

void output_data(int case_num, const Board& board, const Answer& ans) {
  //cerr << "output_data: case_num = " << case_num << endl;
  vector<int> path;
  path.push_back(ans.points[0]);
  for (size_t i = 0; i < ans.points.size() - 1; i++) {
    auto tmp = board.get_path(ans.points[i], ans.points[i + 1]);
    tmp.erase(tmp.begin()); // 最初の点はすでに追加済み
    path.insert(path.end(), tmp.begin(), tmp.end());
  }

  //cerr << "path.size() = " << path.size() << endl;

  const char dc[4] = { 'U', 'D', 'L', 'R' };
  string directions;
  for (size_t i = 0; i < path.size() - 1; i++) {
    int u = path[i], v = path[i + 1];
    int di = v / board.n - u / board.n;
    int dj = v % board.n - u % board.n;
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
      return;
    }
  }

  if (exec_mode == 0) {
    // 標準出力
    cout << directions << endl;
  }
  else {
    // ファイル出力
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofstream ofs(oss.str());

    ofs << directions << endl;

    if (ofs.is_open()) {
      ofs.close();
    }
  }
}

void build_initial_path(const Board& board, Answer& ans) {
  ans.clear();
  ans.vertices.clear();
  ans.points.push_back(board.si * board.n + board.sj); // スタート位置を追加
  ans.vertices.push_back(-1); // スタート位置の頂点は-1とする

  int cur_index = board.si * board.n + board.sj;
  for (int i = 0; i < board.vertices.size(); i++) {
    int next_index = board.dist_to_vertices[cur_index][i].point_index;
    ans.points.push_back(next_index); // 次の頂点への最短距離を追加
    ans.vertices.push_back(i); // 次の頂点のインデックスを追加
    cur_index = next_index; // 現在のインデックスを更新
  }
  ans.points.push_back(board.si * board.n + board.sj); // スタート位置に戻る
  ans.vertices.push_back(-1); // スタート位置の頂点は-1とする
}

void run_simulated_annealing(double time_limit, const Board& board, Answer& ans) {
  ans.calc_score(board);
  Answer best_ans = ans; // ベスト解を初期化

  double start_time = get_elapsed_time();
  double now_time = get_elapsed_time();
  const double START_TEMP = 1800;
  const double END_TEMP = 0.0;

  vector<int> keep_vec(ans.vertices.size());

  int loop = 0;
  while (true) {
    loop++;

    if (loop % 100 == 0) {
      now_time = get_elapsed_time();
      if (now_time > time_limit) break;
    }

    double progress_ratio = (now_time - start_time) / (time_limit - start_time);
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progress_ratio;

    // 近傍解作成
    int ra_exec_mode = rand_xorshift() % 300;
    int ra1, ra2, ra3, ra4, ra5;
    int keep1, keep2, keep3, keep4, keep5;

    int start_index = 0;
    int end_index = 1001001;
    if (ra_exec_mode < 50) {
      // 2点の入れ替え
      ra1 = rand_xorshift() % (ans.points.size() - 2) + 1;
      ra2 = rand_xorshift() % (ans.points.size() - 2) + 1;
      while (ra1 == ra2) {
        ra2 = rand_xorshift() % (ans.points.size() - 2) + 1;
      }
      if (ra1 > ra2) swap(ra1, ra2);
      swap(ans.vertices[ra1], ans.vertices[ra2]); // 2点の入れ替え
      start_index = ra1 - 1;
      end_index = ra2 + 1;
    }
    else if (ra_exec_mode < 100) {
      // 区間reverse
      ra1 = rand_xorshift() % (ans.points.size() - 2) + 1;
      ra2 = rand_xorshift() % (ans.points.size() - 2) + 1;
      while (ra1 == ra2) {
        ra2 = rand_xorshift() % (ans.points.size() - 2) + 1;
      }
      if (ra1 > ra2) swap(ra1, ra2);

      reverse(ans.vertices.begin() + ra1, ans.vertices.begin() + ra2 + 1); // 区間reverse
      start_index = ra1 - 1;
      end_index = ra2 + 1;
    }
    else if (ra_exec_mode < 300) {
      // 1点抜き挿し
      ra1 = rand_xorshift() % (ans.points.size() - 2) + 1;
      ra2 = rand_xorshift() % (ans.points.size() - 2) + 1;
      while (ra1 == ra2) {
        ra2 = rand_xorshift() % (ans.points.size() - 2) + 1;
      }
      // ra1をra2の位置に移動
      if (ra1 < ra2) {
        for (int i = ra1; i < ra2; i++) {
          swap(ans.vertices[i], ans.vertices[i + 1]); // 1点ずつ右にずらす
        }
      }
      else {
        for (int i = ra1; i > ra2; i--) {
          swap(ans.vertices[i], ans.vertices[i - 1]); // 1点ずつ左にずらす
        }
      }
      start_index = min(ra1, ra2) - 1;
      end_index = max(ra1, ra2) + 1;
    }

    // スコア計算
    double current_score = ans.score;
    ans.recalc_points(board, start_index, end_index); // 点の再計算
    double tmp_score = ans.calc_score(board);

    // 焼きなましで採用判定
    double diff_score = (tmp_score - current_score) * 12345.6;
    double prob = exp(diff_score / temp);
    if (prob > rand_01() || rand_xorshift() % 10000 == 0) {
      // 採用
      current_score = tmp_score;

      // ベスト更新
      if (current_score > best_ans.score) {
        best_ans = ans;
      }
    }
    else {
      // 元に戻す
      if (ra_exec_mode < 50) {
        // 近傍操作1 の巻き戻し
        swap(ans.vertices[ra1], ans.vertices[ra2]); // 2点の入れ替えを元に戻す
        ans.recalc_points(board); // 点の再計算
        ans.score = current_score; // スコアを元に戻す
      }
      else if (ra_exec_mode < 100) {
        // 近傍操作2 の巻き戻し
        reverse(ans.vertices.begin() + ra1, ans.vertices.begin() + ra2 + 1); // 区間reverseを元に戻す
        ans.recalc_points(board); // 点の再計算
        ans.score = current_score; // スコアを元に戻す
      }
      else if (ra_exec_mode < 300) {
        // 近傍操作3 の巻き戻し
        if (ra1 < ra2) {
          for (int i = ra2; i > ra1; i--) {
            swap(ans.vertices[i], ans.vertices[i - 1]); // 1点ずつ左にずらす
          }
        }
        else {
          for (int i = ra2; i < ra1; i++) {
            swap(ans.vertices[i], ans.vertices[i + 1]); // 1点ずつ右にずらす
          }
        }
        ans.recalc_points(board, start_index, end_index); // 点の再計算
        ans.score = current_score; // スコアを元に戻す
      }
    }
  }

  if (exec_mode != 3) {
    cerr << loop << endl;
  }

  ans = best_ans; // ベスト解を最終解に設定
}


ll solve_case(int case_num) {
  const double TIME_LIMIT = 2.9;
  start_timer();

  Board board = input_data(case_num);

  Answer ans;
  build_initial_path(board, ans);

  Answer initial_ans = ans;
  Answer best_ans = ans;
  const int SET_COUNT = 1;
  double start_time = get_elapsed_time();
  for (int i = 0; i < SET_COUNT; i++) {
    double time_limit = start_time + (TIME_LIMIT - start_time) / SET_COUNT * (i + 1);
    ans = initial_ans;
    run_simulated_annealing(time_limit, board, ans);
    cerr << "Set " << i + 1 << ": score = " << ans.score << ", time = " << get_elapsed_time() << endl;
    if (ans.score > best_ans.score) {
      best_ans = ans;
    }
  }
  ans = best_ans;

  output_data(case_num, board, ans);

  ll score = 0;
  if (exec_mode != 0) {
    score = ans.calc_score(board);
  }
  return score;
}

int main() {
  exec_mode = 2;

  if (exec_mode == 0) {
    solve_case(0);
  }
  else if (exec_mode <= 2) {
    ll sum_score = 0;
    for (int i = 0; i < 15; i++) {
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
