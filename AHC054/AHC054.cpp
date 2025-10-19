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

// 2次元キューのクラス
class Queue2D
{
private:
  static const int MAX_SIZE = 10000;
  int arr[MAX_SIZE][2];
  int head;
  int tail;

public:
  // コンストラクタ
  Queue2D() : head(0), tail(0) {}

  void clear_queue()
  {
    head = 0;
    tail = 0;
  }

  int front_x() const
  {
    return arr[head][0];
  }

  int front_y() const
  {
    return arr[head][1];
  }

  void push(int x, int y)
  {
    arr[tail][0] = x;
    arr[tail][1] = y;
    tail++;
  }

  void pop()
  {
    head++;
  }

  int size() const
  {
    return tail - head;
  }

  bool empty() const
  {
    return size() == 0;
  }
};

Queue2D que2D;

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

  std::random_device seed_gen;
  std::mt19937 engine(seed_gen());
}

// 上下左右
const int DX[4] = { -1, 1, 0, 0 };
const int DY[4] = { 0, 0, -1, 1 };

const double TIME_LIMIT = 1.8;
int exec_mode;
bool is_simulate;

const int MAX_N = 40;

int n;
int ti, tj;
int b[MAX_N][MAX_N];
int init_b[MAX_N][MAX_N];

int turn;
int confirmed[MAX_N][MAX_N];
int pi, pj;

// ローカルテスター用 入力
int current_q;
int qi[MAX_N * MAX_N], qj[MAX_N * MAX_N];
int init_qi[MAX_N * MAX_N], init_qj[MAX_N * MAX_N];
int generate_qi[10][MAX_N * MAX_N], generate_qj[10][MAX_N * MAX_N];
int ri, rj;

// 高速化用
int bfs_dp[MAX_N][MAX_N];
int bfs_visited[MAX_N][MAX_N];
int bfs_unknown_reachable[MAX_N][MAX_N];
int bfs_version;

void generate_q()
{
  vector<P> vp;
  int cnt = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      vp.push_back(P(i, j));
    }
  }

  // シャッフル
  std::shuffle(vp.begin(), vp.end(), engine);

  for (int i = 0; i < n * n; i++) {
    qi[i] = vp[i].first;
    qj[i] = vp[i].second;
  }
}

void generate_q10()
{
  for (int loop = 0; loop < 10; loop++) {
    vector<P> vp;
    int cnt = 0;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        vp.push_back(P(i, j));
      }
    }

    // シャッフル
    std::shuffle(vp.begin(), vp.end(), engine);

    for (int i = 0; i < n * n; i++) {
      generate_qi[loop][i] = vp[i].first;
      generate_qj[loop][i] = vp[i].second;
    }
  }
}

void setup_true_q()
{
  for (int i = 0; i < n * n - 1; i++) {
    qi[i] = init_qi[i];
    qj[i] = init_qj[i];
  }
}

void setup_generate_q10(int num)
{
  for (int i = 0; i < n * n - 1; i++) {
    qi[i] = generate_qi[num][i];
    qj[i] = generate_qj[num][i];
  }
}

void reset_board()
{
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      b[i][j] = init_b[i][j];
    }
  }
}

void initialize_state()
{
  current_q = 0;
  ri = -1;
  rj = -1;

  turn = 0;
  pi = 0;
  pj = n / 2;
  for (int i = 0; i < MAX_N; i++) {
    for (int j = 0; j < MAX_N; j++) {
      confirmed[i][j] = 0;
    }
  }

  bfs_version = 0;
  for (int i = 0; i < MAX_N; i++) {
    for (int j = 0; j < MAX_N; j++) {
      bfs_dp[i][j] = -1;
      bfs_visited[i][j] = 0;
    }
  }
}

void initialize_simulate()
{
  current_q = 0;
  ri = -1;
  rj = -1;

  turn = 0;
  pi = 0;
  pj = n / 2;
  for (int i = 0; i < MAX_N; i++) {
    for (int j = 0; j < MAX_N; j++) {
      confirmed[i][j] = 0;
    }
  }

  bfs_version = 0;
  for (int i = 0; i < MAX_N; i++) {
    for (int j = 0; j < MAX_N; j++) {
      bfs_dp[i][j] = -1;
      bfs_visited[i][j] = 0;
    }
  }

  reset_board();
}

void input_data(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    cin >> n;
    cin >> ti >> tj;
    for (int i = 0; i < n; i++) {
      string str;
      cin >> str;
      for (int j = 0; j < n; j++) {
        if (str[j] == '.') {
          b[i][j] = 0;
        }
        else {
          b[i][j] = 1;
        }
      }
    }
  }
  else {
    // ファイル入力
    ifs >> n;
    ifs >> ti >> tj;
    for (int i = 0; i < n; i++) {
      string str;
      ifs >> str;
      for (int j = 0; j < n; j++) {
        if (str[j] == '.') {
          b[i][j] = 0;
        }
        else {
          b[i][j] = 1;
        }
      }
    }

    for (int i = 0; i < n * n - 1; i++) {
      ifs >> qi[i] >> qj[i];
    }

    for (int i = 0; i < n * n - 1; i++) {
      init_qi[i] = qi[i];
      init_qj[i] = qj[i];
    }

    ifs.close();
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      init_b[i][j] = b[i][j];
    }
  }
}

void open_ofs(int case_num, ofstream& ofs)
{
  if (exec_mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofs.open(oss.str());
  }
}

// ゴールに到着可能か判定
bool can_reach_goal()
{
  bfs_version++;
  que2D.clear_queue();
  que2D.push(pi, pj);
  bfs_visited[pi][pj] = bfs_version;
  while (!que2D.empty()) {
    int i = que2D.front_x();
    int j = que2D.front_y();
    que2D.pop();
    if (i == ti && j == tj) {
      return true;
    }
    for (int d = 0; d < 4; d++) {
      int ni = i + DX[d];
      int nj = j + DY[d];
      if (ni < 0 || ni >= n || nj < 0 || nj >= n) {
        continue;
      }
      if (b[ni][nj] == 1) {
        continue;
      }
      if (bfs_visited[ni][nj] == bfs_version) {
        continue;
      }
      bfs_visited[ni][nj] = bfs_version;
      que2D.push(ni, nj);
    }
  }
  return false;
}

bool can_reach_all_cells(int margin)
{
  if (b[ti][tj] == 1) {
    return false;
  }

  bfs_version++;
  que2D.clear_queue();
  que2D.push(pi, pj);
  bfs_visited[pi][pj] = bfs_version;
  while (!que2D.empty()) {
    int i = que2D.front_x();
    int j = que2D.front_y();
    que2D.pop();
    for (int d = 0; d < 4; d++) {
      int ni = i + DX[d];
      int nj = j + DY[d];
      if (ni < 0 || ni >= n || nj < 0 || nj >= n) {
        continue;
      }
      if (b[ni][nj] == 1) {
        continue;
      }
      if (bfs_visited[ni][nj] == bfs_version) {
        continue;
      }
      bfs_visited[ni][nj] = bfs_version;
      que2D.push(ni, nj);
    }
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (b[i][j] == 0 && bfs_visited[i][j] != bfs_version) {
        margin--;
      }
    }
  }

  if (bfs_visited[ti][tj] != bfs_version) {
    margin = -1;
  }

  return margin >= 0;
}

// 未確認のマスをすべて空きマスであると仮定してBFS
int get_next_heros_move()
{
  bfs_version++;
  que2D.clear_queue();
  que2D.push(pi, pj);
  bfs_visited[pi][pj] = bfs_version;
  bfs_dp[pi][pj] = 0;
  while (!que2D.empty()) {
    int i = que2D.front_x();
    int j = que2D.front_y();
    que2D.pop();
    if (i == ri && j == rj) {
      break;
    }
    for (int d = 0; d < 4; d++) {
      int ni = i + DX[d];
      int nj = j + DY[d];
      if (ni < 0 || ni >= n || nj < 0 || nj >= n) {
        continue;
      }
      if (confirmed[ni][nj] == 1 && b[ni][nj] == 1) {
        continue;
      }
      int next_cost = bfs_dp[i][j] + 4;
      if (bfs_dp[i][j] == 0) {
        next_cost += d;
      }
      if (bfs_visited[ni][nj] == bfs_version && bfs_dp[ni][nj] <= next_cost) {
        continue;
      }
      bfs_visited[ni][nj] = bfs_version;
      bfs_dp[ni][nj] = next_cost;
      que2D.push(ni, nj);
    }
  }

  return bfs_dp[ri][rj] % 4;
}

bool is_reachable_unknown(int i, int j)
{
  return bfs_visited[i][j] == bfs_version && bfs_unknown_reachable[i][j] == 1;
}

// 未確認のマスをすべて空きマスであると仮定して到達可能なすべての未確認マスをBFSで求める
void update_unknown_reachability()
{
  bfs_version++;
  que2D.clear_queue();
  que2D.push(pi, pj);
  bfs_visited[pi][pj] = bfs_version;
  bfs_unknown_reachable[pi][pj] = 0;
  while (!que2D.empty()) {
    int i = que2D.front_x();
    int j = que2D.front_y();
    que2D.pop();
    for (int d = 0; d < 4; d++) {
      int ni = i + DX[d];
      int nj = j + DY[d];
      if (ni < 0 || ni >= n || nj < 0 || nj >= n) {
        continue;
      }
      if (confirmed[ni][nj] == 1 && b[ni][nj] == 1) {
        continue;
      }
      if (bfs_visited[ni][nj] == bfs_version) {
        continue;
      }
      bfs_visited[ni][nj] = bfs_version;
      bfs_unknown_reachable[ni][nj] = 1 - confirmed[ni][nj];
      que2D.push(ni, nj);
    }
  }
}

void update_confirmed_local_tester()
{
  if (turn == 0) {
    pi = 0;
    pj = n / 2;
    confirmed[pi][pj] = 1;
  }
  else {
    for (int d = 0; d < 4; d++) {
      int ni = pi;
      int nj = pj;
      while (true) {
        ni += DX[d];
        nj += DY[d];
        if (ni < 0 || ni >= n || nj < 0 || nj >= n) {
          break;
        }
        confirmed[ni][nj] = 1;
        if (b[ni][nj] == 1) {
          break;
        }
      }
    }

    if (ri == ti && rj == tj) {
      // 目的地イコールGoalになったら今後目的地が変更されることはない
    }
    else if (confirmed[ti][tj] == 1) {
      ri = ti;
      rj = tj;
    }
    else {
      update_unknown_reachability();
      if (ri != -1 && !is_reachable_unknown(ri, rj)) {
        ri = -1;
        rj = -1;
      }
      if (ri != -1 && confirmed[ri][rj] == 1) {
        ri = -1;
        rj = -1;
      }
      while (ri == -1) {
        int ni = qi[current_q];
        int nj = qj[current_q];
        current_q++;
        if (is_reachable_unknown(ni, nj)) {
          ri = ni;
          rj = nj;
          break;
        }
      }
    }

    int dir = get_next_heros_move();
    pi += DX[dir];
    pj += DY[dir];
  }
}

void update_confirmed()
{
  if (is_simulate) {
    update_confirmed_local_tester();
  }
  else if (exec_mode == 0) {
    // 標準入力
    cin >> pi >> pj;
    int x_size;
    cin >> x_size;
    vector<int> xs, ys;
    xs.resize(x_size);
    ys.resize(x_size);
    for (int i = 0; i < x_size; i++) {
      cin >> xs[i] >> ys[i];
    }
    for (int i = 0; i < x_size; i++) {
      confirmed[xs[i]][ys[i]] = 1;
    }
  }
  else {
    update_confirmed_local_tester();
  }
}

void put_each_turn_output(ofstream& ofs, const vector<P>& ps)
{
  if (is_simulate) {
    // 何もしない
  }
  else if (exec_mode == 0) {
    // 標準出力
    cout << ps.size();
    for (int i = 0; i < (int)ps.size(); i++) {
      cout << " " << ps[i].first << " " << ps[i].second;
    }
    cout << endl;
  }
  else {
    // ファイル出力
    ofs << ps.size();
    for (int i = 0; i < (int)ps.size(); i++) {
      ofs << " " << ps[i].first << " " << ps[i].second;
    }
    ofs << endl;
  }

  for (int i = 0; i < (int)ps.size(); i++) {
    b[ps[i].first][ps[i].second] = 1;
  }
}

struct SimulateParam
{
  int goal_guard_method = 0;
  int init_method = 0;
  int init_transpose = 0;
  int slide_i = 0;
  int slide_j = 0;
  vector<vector<int>> placed_init_trees;
};

struct SimulateResult
{
  int score = 0;
  vector<vector<int>> placed_init_trees;
};

void attempt(int i, int j, vector<P>& ps, int margin)
{
  if (i < 0 || i >= n || j < 0 || j >= n) {
    return;
  }

  if (b[i][j] == 1 || (i == ti && j == tj) || (i == pi && j == pj)) {
    return;
  }

  b[i][j] = 1;
  if (!can_reach_all_cells(margin)) {
    b[i][j] = 0;
  }
  else {
    ps.emplace_back(i, j);
  }
}

void init_make_goal_guard(vector<P>& ps, const SimulateParam& param)
{
  if (param.goal_guard_method == 0) {
    // 左上
    attempt(ti, tj + 1, ps, 5);
    for (int k = 0; k < 100; k++) {
      int i1 = ti - 1 - k;
      int j1 = tj + 1 - k;
      attempt(i1, j1, ps, 5);
      int i2 = ti + 1 - k;
      int j2 = tj - k;
      attempt(i2, j2, ps, 5);
    }
  }
  else if (param.goal_guard_method == 1) {
    // 左下
    attempt(ti - 1, tj, ps, 5);
    for (int k = 0; k < 100; k++) {
      int i1 = ti - 1 + k;
      int j1 = tj - 1 - k;
      attempt(i1, j1, ps, 5);
      int i2 = ti + k;
      int j2 = tj + 1 - k;
      attempt(i2, j2, ps, 5);
    }
  }
  else if (param.goal_guard_method == 2) {
    // 右下
    attempt(ti, tj - 1, ps, 5);
    for (int k = 0; k < 100; k++) {
      int i1 = ti + 1 + k;
      int j1 = tj - 1 + k;
      attempt(i1, j1, ps, 5);
      int i2 = ti - 1 + k;
      int j2 = tj + k;
      attempt(i2, j2, ps, 5);
    }
  }
  else if (param.goal_guard_method == 3) {
    // 右上
    attempt(ti + 1, tj, ps, 5);
    for (int k = 0; k < 100; k++) {
      int i1 = ti + 1 - k;
      int j1 = tj + 1 + k;
      attempt(i1, j1, ps, 5);
      int i2 = ti - k;
      int j2 = tj - 1 + k;
      attempt(i2, j2, ps, 5);
    }
  }
}

// helper: 配置判定ロジック
bool should_put(int ii, int jj, int n, int kind)
{
  switch (kind) {
    case 0: // init_0
      switch (ii % 6) {
        case 0: return jj % 3 == 2;
        case 1: return jj % 3 == 0;
        case 2: return false;
        case 3: return jj % 3 == 0;
        case 4: return jj % 3 == 2;
        case 5:
          if (ii % 12 == 5) return jj != n - 1;
          else return jj != 0;
      }
      break;
    case 1: // init_1
      switch (ii % 6) {
        case 0: return jj % 3 == 2;
        case 1: return jj % 3 == 0;
        case 2: return jj % 3 == 1;
        case 3: return jj % 3 == 2;
        case 4: return false;
        case 5:
          if (ii % 12 == 5) return jj != n - 1;
          else return jj != 0;
      }
      break;
    case 2: // init_2
      switch (ii % 4) {
        case 0: return jj % 3 == 2;
        case 1: return jj % 3 == 0;
        case 2: return false;
        case 3:
          if (ii % 8 == 3) return jj != n / 2;
          else return jj != 0;
      }
      break;
    default:
      cerr << "Error: unknown init_method " << kind << endl;
      break;
  }
  return false;
}

// 共通関数
vector<P> init_pattern(const SimulateParam& param)
{
  bool transpose = (param.init_transpose == 1);
  int kind = param.init_method;

  vector<P> ps;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (confirmed[i][j] == 1 || b[i][j] == 1 || (i == ti && j == tj) || (i == pi && j == pj)) {
        continue;
      }

      int ii = (i + param.slide_i + n) % n;
      int jj = (j + param.slide_j + n) % n;

      // 縦横を反転させる
      if (transpose) std::swap(ii, jj);

      if (should_put(ii, jj, n, kind)) {
        ps.emplace_back(i, j);
      }
    }
  }
  return ps;
}

int ng_count;
vector<P> init_pattern2(const SimulateParam& param)
{
  const int M = n / 4;
  vector<int> hs(M + 1), ws(M + 1);
  for (int i = 0; i <= M; i++) {
    hs[i] = i * 4;
    ws[i] = i * 4;
  }
  for (int _loop = 0; _loop < n % 4; _loop++) {
    int i = rand_xorshift() % M + 1;
    for (int j = i; j <= M; j++) {
      hs[j]++;
    }
    i = rand_xorshift() % M + 1;
    for (int j = i; j <= M; j++) {
      ws[j]++;
    }
  }

  vector<vector<int>> grid(M, vector<int>(M, -1));
  int si = -1;
  int sj = -1;
  int gi = -1;
  int gj = -1;
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < M; j++) {
      if (hs[i] <= 0 && 0 < hs[i + 1] && ws[j] <= n / 2 && n / 2 < ws[j + 1]) {
        si = i;
        sj = j;
      }
      if (hs[i] <= ti && ti < hs[i + 1] && ws[j] <= tj && tj < ws[j + 1]) {
        gi = i;
        gj = j;
      }
    }
  }

  // 方向をシャッフル
  vector<int> dx, dy;
  {
    vector<P> dirs = { P(-1,0), P(1,0), P(0,-1), P(0,1) };
    std::shuffle(dirs.begin(), dirs.end(), engine);
    for (auto d : dirs) {
      dx.push_back(d.first);
      dy.push_back(d.second);
    }
  }

  const int TOTAL = M * M;

  auto inside_blk = [&](int x, int y) -> bool {
    return (0 <= x && x < M && 0 <= y && y < M);
    };
  auto parity = [&](int x, int y) -> int { return (x + y) & 1; };

  // Hamiltonian のパリティ条件:
  // TOTAL が偶数 => 終端2点は異色, TOTAL が奇数 => 終端2点は同色
  const bool hamilton_ok =
    ((TOTAL % 2 == 0) ? (parity(si, sj) != parity(gi, gj))
      : (parity(si, sj) == parity(gi, gj)));

  const int targetLen = hamilton_ok ? TOTAL : (TOTAL - 1); // 塗る個数
  bool found = false;

  // 経路（ブロック座標の列）。grid[i][j] には訪問順を入れる（0..targetLen-1）
  vector<P> order; order.reserve(targetLen);

  double now_time = timer.get_elapsed_time();

  function<void(int, int, int)> dfs = [&](int x, int y, int step) {
    if (timer.get_elapsed_time() - now_time > 0.010) {
      return;
    }
    grid[x][y] = step;
    order.emplace_back(x, y);

    if (step == targetLen - 1) {
      // 最後のマスは gi,gj であることを要求
      if (x == gi && y == gj) {
        found = true;
        return;
      }
    }
    else if (step == targetLen - 2 && x == gi && y == gj) {
      found = true;
      return;
    }
    else {
      // 未探索マスが二つ以上に分かれてたら探索する意味がないので打ち切り
      int cnt = 0;
      bfs_version++;
      que2D.clear_queue();
      que2D.push(gi, gj);
      bfs_visited[gi][gj] = bfs_version;
      cnt++;
      while (!que2D.empty()) {
        int i = que2D.front_x();
        int j = que2D.front_y();
        que2D.pop();
        for (int d = 0; d < 4; d++) {
          int ni = i + dx[d];
          int nj = j + dy[d];
          if (!inside_blk(ni, nj) || grid[ni][nj] != -1) {
            continue;
          }
          if (bfs_visited[ni][nj] == bfs_version) {
            continue;
          }
          bfs_visited[ni][nj] = bfs_version;
          cnt++;
          que2D.push(ni, nj);
        }
      }
      if (cnt + step < targetLen - 1) {
        // 残りマス数が足りない
        grid[x][y] = -1;
        order.pop_back();
        return;
      }

      // dx,dy の優先順で貪欲に試す（基本は順序通り。解が出ないときだけバックトラック）
      for (int k = 0; k < 4 && !found; ++k) {
        int nx = x + dx[k], ny = y + dy[k];
        if (!inside_blk(nx, ny) || grid[nx][ny] != -1) {
          continue;
        }

        // ゴールには「最後の一手」でしか入らない
        if (nx == gi && ny == gj && step + 1 != targetLen - 1) {
          continue;
        }

        // 軽い安全策：残り 1 手未満で“袋小路”になりそうならスキップ
        int deg = 0;
        for (int t = 0; t < 4; ++t) {
          int ux = nx + dx[t], uy = ny + dy[t];
          if (inside_blk(ux, uy) && grid[ux][uy] == -1) {
            ++deg;
          }
        }
        if (deg == 0 && !(nx == gi && ny == gj && step + 1 == targetLen - 1)) {
          continue;
        }

        dfs(nx, ny, step + 1);
      }
    }

    if (found) {
      return;
    }
    // バックトラック
    grid[x][y] = -1;
    order.pop_back();
    };

  fill(grid.begin(), grid.end(), vector<int>(M, -1));
  order.clear();
  now_time = timer.get_elapsed_time();
  dfs(si, sj, 0);

  if (!found) {
    ng_count++;
    //cerr << "Warning: Hamiltonian path not found." << endl;
    return vector<P>();
  }

  //for (int i = 0; i < M; i++) {
  //  for (int j = 0; j < M; j++) {
  //    cout << setw(2) << grid[i][j] << " ";
  //  }
  //  cout << endl;
  //}

  // -1の処理
  vector<vector<bool>> isMinusOne(M, vector<bool>(M, false));
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < M; j++) {
      if (grid[i][j] == -1) {
        isMinusOne[i][j] = true;
        int min_num = 99999;
        for (int d = 0; d < 4; d++) {
          int ni = i + DX[d];
          int nj = j + DY[d];
          if (inside_blk(ni, nj) && grid[ni][nj] != -1) {
            min_num = min(min_num, grid[ni][nj]);
          }
        }
        if (min_num == 99999) {
          cerr << "Error: isolated cell in Hamiltonian path." << endl;
          continue;
        }
        order.insert(order.begin() + min_num, P(i, j));
        for (int k = min_num + 1; k < (int)order.size(); k++) {
          int x = order[k].first;
          int y = order[k].second;
          grid[x][y] = k;
        }
      }
    }
  }

  vector<P> ps;
  vector<P> ps_toge;

  // 迷路のように壁を作る
  for (int idx = 0; idx < (int)order.size(); idx++) {
    int bi = order[idx].first;
    int bj = order[idx].second;

    // 上下左右
    vector<P> walls[4];
    for (int j = ws[bj]; j < ws[bj + 1]; j++) {
      walls[0].emplace_back(hs[bi], j);
      walls[1].emplace_back(hs[bi + 1] - 1, j);
    }
    for (int i = hs[bi]; i < hs[bi + 1]; i++) {
      walls[2].emplace_back(i, ws[bj]);
      walls[3].emplace_back(i, ws[bj + 1] - 1);
    }

    int togeDir = -1;
    int min_toge = 99999;
    bool is_put_wall[4] = { false,false,false,false };

    if (isMinusOne[bi][bj]) {
      // 自分より小さい数で一番大きい数を探す
      int ref_dir = -1;
      int ref_max = -1;
      for (int k = 0; k < 4; k++) {
        int ni = bi + DX[k];
        int nj = bj + DY[k];
        if (inside_blk(ni, nj) && grid[ni][nj] != -1 && grid[ni][nj] < grid[bi][bj]) {
          if (ref_max < grid[ni][nj]) {
            ref_max = grid[ni][nj];
            ref_dir = k;
          }
        }
      }

      // ref_dir以外の方向には壁を作る
      for (int k = 0; k < 4; k++) {
        if (k == ref_dir) {
          continue;
        }

        // 外側には壁を作る必要はない
        int ni = bi + DX[k];
        int nj = bj + DY[k];
        if (!inside_blk(ni, nj)) {
          if (k < 2) {
            min_toge = -2;
            togeDir = k;
          }
          else {
            if (-1 < min_toge) {
              min_toge = -1;
              togeDir = k;
            }
          }
          continue;
        }

        if (grid[ni][nj] != -1 && grid[ni][nj] < min_toge) {
          min_toge = grid[ni][nj];
          togeDir = k;
        }

        for (auto p : walls[k]) {
          int i = p.first;
          int j = p.second;
          ps.emplace_back(i, j);
        }
        is_put_wall[k] = true;
      }
    }
    else {
      // 自分より大きい数で一番小さい数を探す
      int ref_dir = -1;
      int ref_min = 99999;
      for (int k = 0; k < 4; k++) {
        int ni = bi + DX[k];
        int nj = bj + DY[k];

        if (!inside_blk(ni, nj)) {
          continue;
        }

        if (grid[ni][nj] == -1 || grid[ni][nj] < grid[bi][bj]) {
          continue;
        }
        if (isMinusOne[ni][nj]) {
          continue;
        }

        if (grid[ni][nj] < ref_min) {
          ref_min = grid[ni][nj];
          ref_dir = k;
        }
      }

      // 進む方向以外の-1以外の自分より大きい数の方向には壁を作る
      for (int k = 0; k < 4; k++) {
        if (k == ref_dir) {
          continue;
        }
        // 外側には壁を作る必要はない
        int ni = bi + DX[k];
        int nj = bj + DY[k];
        if (!inside_blk(ni, nj)) {
          if (k < 2) {
            min_toge = -2;
            togeDir = k;
          }
          else {
            if (-1 < min_toge) {
              min_toge = -1;
              togeDir = k;
            }
          }
          continue;
        }

        if (grid[ni][nj] == -1) {
          continue;
        }
        if (isMinusOne[ni][nj]) {
          continue;
        }

        if (grid[ni][nj] < grid[bi][bj]) {
          if (grid[ni][nj] < min_toge) {
            min_toge = grid[ni][nj];
            togeDir = k;
          }
          continue;
        }

        for (auto p : walls[k]) {
          int i = p.first;
          int j = p.second;
          ps.emplace_back(i, j);
        }
        is_put_wall[k] = true;
      }
    }

    // 一方向にだけとげを付ける
    if (togeDir == 0) {
      for (int j = ws[bj]; j < ws[bj + 1]; j++) {
        int ii = hs[bi];
        int jj = j;
        if (is_put_wall[0]) {
          ii++;
        }
        if ((ii + jj) % 3 == 2) {
          ps_toge.emplace_back(ii, jj);
          ps_toge.emplace_back(ii + 1, jj - 1);
        }
      }
    }
    else if (togeDir == 1) {
      for (int j = ws[bj]; j < ws[bj + 1]; j++) {
        int ii = hs[bi + 1] - 1;
        int jj = j;
        if (is_put_wall[1]) {
          ii--;
        }
        if ((ii + jj) % 3 == 2) {
          ps_toge.emplace_back(ii, jj);
          ps_toge.emplace_back(ii - 1, jj + 1);
        }
      }
    }
    else if (togeDir == 2) {
      for (int i = hs[bi]; i < hs[bi + 1]; i++) {
        int ii = i;
        int jj = ws[bj];
        if (is_put_wall[2]) {
          jj++;
        }
        if ((ii + jj) % 3 == 2) {
          ps_toge.emplace_back(ii, jj);
          ps_toge.emplace_back(ii - 1, jj + 1);
        }
      }
    }
    else if (togeDir == 3) {
      for (int i = hs[bi]; i < hs[bi + 1]; i++) {
        int ii = i;
        int jj = ws[bj + 1] - 1;
        if (is_put_wall[3]) {
          jj--;
        }
        if ((ii + jj) % 3 == 2) {
          ps_toge.emplace_back(ii, jj);
          ps_toge.emplace_back(ii + 1, jj - 1);
        }
      }
    }
  }

  for (auto p : ps_toge) {
    ps.emplace_back(p);
  }

  return ps;
}

vector<P> init_common(const SimulateParam& param)
{
  vector<P> ps;

  if (is_simulate) {
    if (true) {
      init_make_goal_guard(ps, param);

      vector<P> tmp_ps = init_pattern(param);

      for (auto p : tmp_ps) {
        int i = p.first;
        int j = p.second;
        attempt(i, j, ps, 5);
      }
    }
    else {
      vector<P> tmp_ps;

      int try_count = 10;
      while (try_count && tmp_ps.empty()) {
        tmp_ps = init_pattern2(param);
        try_count--;
      }

      for (auto p : tmp_ps) {
        int i = p.first;
        int j = p.second;
        attempt(i, j, ps, 5);
      }

      init_make_goal_guard(ps, param);
    }
  }
  else {
    // paramのplaced_init_treesをそのまま使う
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (param.placed_init_trees[i][j] == 1) {
          attempt(i, j, ps, 5);
        }
      }
    }
  }

  return ps;
}

bool recreate_board()
{
  vector<P> vp;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (confirmed[i][j] == 0 && b[i][j] == 1 && init_b[i][j] == 0) {
        b[i][j] = 0;
        vp.emplace_back(i, j);
      }
    }
  }

  if (!can_reach_all_cells(99999)) {
    for (auto p : vp) {
      int i = p.first;
      int j = p.second;
      b[i][j] = 1;
    }
    return false;
  }

  vector<P> ps;
  for (auto p : vp) {
    int i = p.first;
    int j = p.second;
    attempt(i, j, ps, 5);
  }
  return true;
}

vector<P> put_this_turn()
{
  vector<P> res;

  // ゴールが見つかってしまうかどうかチェック
  if (confirmed[ti][tj] == 0) {
    int emergencyDir = -1;
    for (int d = 0; d < 4; d++) {
      int ni = pi;
      int nj = pj;
      while (true) {
        ni += DX[d];
        nj += DY[d];
        if (ni < 0 || ni >= n || nj < 0 || nj >= n) {
          break;
        }
        if (b[ni][nj] == 1) {
          break;
        }
        if (ni == ti && nj == tj) {
          emergencyDir = d;
          break;
        }
      }
      if (emergencyDir != -1) {
        break;
      }
    }

    if (emergencyDir != -1) {
      int d = emergencyDir;
      int ni = pi;
      int nj = pj;
      while (true) {
        ni += DX[d];
        nj += DY[d];
        if (ni == ti && nj == tj) {
          break;
        }

        // 確認済みの場合木は置けないのでcontinue
        if (confirmed[ni][nj] == 1) {
          continue;
        }

        // 木を置けるなら置く
        confirmed[ni][nj] = 1;
        b[ni][nj] = 1;
        update_unknown_reachability();
        if (is_reachable_unknown(ti, tj) && recreate_board()) {
          confirmed[ni][nj] = 0;
          break;
        }
        else {
          confirmed[ni][nj] = 0;
          b[ni][nj] = 0;
        }
      }
    }
  }

  for (int d = 0; d < 4; d++) {
    int ni = pi;
    int nj = pj;
    while (true) {
      ni += DX[d];
      nj += DY[d];
      if (ni < 0 || ni >= n || nj < 0 || nj >= n) {
        break;
      }
      if (b[ni][nj] == 1) {
        if (confirmed[ni][nj] == 0 && init_b[ni][nj] == 0) {
          res.emplace_back(ni, nj);
        }
        break;
      }
    }
  }

  return res;
}

SimulateResult method1(ofstream& ofs, const SimulateParam& param)
{
  SimulateResult result;
  result.placed_init_trees.resize(n, vector<int>(n, 0));

  // まずは適当に思いついたルールで配置を作ってみる
  init_common(param);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (b[i][j] == 1 && init_b[i][j] == 0) {
        result.placed_init_trees[i][j] = 1;
      }
    }
  }

  //for (int i = 0; i < n; i++) {
  //  for (int j = 0; j < n; j++) {
  //    if (b[i][j] == 1 || init_b[i][j] == 1) {
  //      cerr << "1";
  //    }
  //    else {
  //      if (i == ti && j == tj) {
  //        cerr << "G";
  //      }
  //      else if (i == 0 && j == n / 2) {
  //        cerr << "S";
  //      }
  //      else {
  //        cerr << "0";
  //      }
  //    }
  //  }
  //  cerr << endl;
  //}

  while (true) {
    update_confirmed();
    if (pi == ti && pj == tj) {
      break;
    }

    vector<P> ps;
    if (false) {
      ps = put_this_turn();
    }
    else {
      // 最初のターンに全部置く
      if (turn == 0) {
        for(int i = 0; i < n; i++) {
          for(int j = 0; j < n; j++) {
            if (confirmed[i][j] == 0 && b[i][j] == 1 && init_b[i][j] == 0) {
              ps.push_back(P(i, j));
            }
          }
        }
      }
    }

    put_each_turn_output(ofs, ps);

    turn++;
  }

  result.score = turn;
  return result;
}


int solve_case(int case_num)
{
  timer.start();

  initialize_state();

  input_data(case_num);

  ofstream ofs;
  open_ofs(case_num, ofs);

  SimulateParam bestParam;
  int bestScore = -1;

  is_simulate = true;
  generate_q10();

  int sim_loop = 0;
  ng_count = 0;
  while (true) {
    if (timer.get_elapsed_time() > TIME_LIMIT) {
      break;
    }
    sim_loop++;

    SimulateParam param;
    param.goal_guard_method = rand_xorshift() % 4;
    param.init_method = rand_xorshift() % 3;
    param.init_transpose = rand_xorshift() % 2;
    param.slide_i = rand_xorshift() % n;
    param.slide_j = rand_xorshift() % n;

    //cerr << "sim_loop = " << setw(4) << sim_loop << endl;
    //init_pattern2(param);
    //cerr << endl;
    //continue;

    int score_sum = 0;
    SimulateResult result;
    for (int loop = 0; loop < 4; loop++) {
      setup_generate_q10(loop);
      initialize_simulate();
      result = method1(ofs, param);
      score_sum += result.score;
    }
    if (score_sum > bestScore) {
      bestScore = score_sum;
      bestParam = param;
      bestParam.placed_init_trees = result.placed_init_trees;
    }
  }

  cerr << "ng_count = " << ng_count << endl;

  cerr
    << "sim_loop = " << setw(4) << sim_loop
    << ", bestScore = " << setw(5) << bestScore
    << ", init_method = " << bestParam.init_method
    << ", init_transpose = " << bestParam.init_transpose
    << ", goal_guard_method = " << bestParam.goal_guard_method
    << endl;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (bestParam.placed_init_trees[i][j] == 1 || init_b[i][j] == 1) {
        cerr << "1";
      }
      else {
        if (i == ti && j == tj) {
          cerr << "G";
        }
        else if (i == 0 && j == n / 2) {
          cerr << "S";
        }
        else {
          cerr << "0";
        }
      }
    }
    cerr << endl;
  }

  is_simulate = false;
  setup_true_q();

  initialize_simulate();
  auto result = method1(ofs, bestParam);

  if (ofs.is_open()) {
    ofs.close();
  }

  return result.score;
}

int main()
{
  exec_mode = 2;

  if (exec_mode == 0) {
    solve_case(0);
  }
  else if (exec_mode <= 2) {
    int sum_score = 0;
    for (int i = 0; i < 100; i++) {
      int score = solve_case(i);
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
