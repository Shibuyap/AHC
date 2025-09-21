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

  void shuffle_array(int* arr, int n)
  {
    for (int i = n - 1; i >= 0; i--) {
      int j = rand_xorshift() % (i + 1);
      int swa = arr[i];
      arr[i] = arr[j];
      arr[j] = swa;
    }
  }

  std::random_device seed_gen;
  std::mt19937 engine(seed_gen());
}

const int DX[4] = { -1, 1, 0, 0 };
const int DY[4] = { 0, 0, -1, 1 };

const double TIME_LIMIT = 1.9;
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
int ri, rj;

// 高速化用
int bfs_dp[MAX_N][MAX_N];
int bfs_visited[MAX_N][MAX_N];
int bfs_can_use[MAX_N][MAX_N];
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

void setup_true_q()
{
  for (int i = 0; i < n * n - 1; i++) {
    qi[i] = init_qi[i];
    qj[i] = init_qj[i];
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
  pi = -1;
  pj = -1;
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
  pi = -1;
  pj = -1;
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

int calculate_score()
{
  int res = 0;
  return res;
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

int get_next_heros_move()
{
  // 未確認のマスをすべて空きマスであると仮定してBFS
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

void get_can_use()
{
  // 未確認のマスをすべて空きマスであると仮定して到達可能なすべての未確認マスをBFSで求める
  bfs_version++;
  que2D.clear_queue();
  que2D.push(pi, pj);
  bfs_visited[pi][pj] = bfs_version;
  bfs_can_use[pi][pj] = 0;
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
      bfs_can_use[ni][nj] = 1 - confirmed[ni][nj];
      que2D.push(ni, nj);
    }
  }
}

// (pi,pj)から到達不可能なマスをすべて木で埋める
vector<P> fill_unreachable_with_trees()
{
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

  vector<P> res;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (b[i][j] == 0 && bfs_visited[i][j] != bfs_version) {
        b[i][j] = 1;
        res.push_back(P(i, j));
      }
    }
  }
  return res;
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

    if(ri == ti && rj == tj) {
      // 目的地イコールGoalになったら今後目的地が変更されることはない
    }
    else if (confirmed[ti][tj] == 1) {
      ri = ti;
      rj = tj;
    }
    else {
      get_can_use();
      if (ri != -1 && (bfs_visited[ri][rj] != bfs_version || bfs_can_use[ri][rj] == 0)) {
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
        if (bfs_visited[ni][nj] == bfs_version && bfs_can_use[ni][nj] == 1 && confirmed[ni][nj] == 0) {
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

struct SimulateParam
{
  int init_method = 0;
  int slide_i = 0;
  int slode_j = 0;
};

vector<P> init_0(const SimulateParam& param)
{
  vector<P> ps;

  attempt(ti, tj + 1, ps, 5);
  for (int k = 0; k < 100; k++) {
    int i1 = ti - 1 - k;
    int j1 = tj + 1 - k;
    attempt(i1, j1, ps, 5);
    int i2 = ti + 1 - k;
    int j2 = tj - k;
    attempt(i2, j2, ps, 5);
  }

  vector<P> tmp_ps;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (confirmed[i][j] == 1 || b[i][j] == 1 || (i == ti && j == tj) || (i == pi && j == pj)) {
        continue;
      }

      int ii = (i + param.slide_i + n) % n;
      int jj = (j + param.slode_j + n) % n;

      bool is_put = false;
      switch (ii % 6) {
      case 0:
        if (jj % 3 == 2) {
          is_put = true;
        }
        break;
      case 1:
        if (jj % 3 == 0) {
          is_put = true;
        }
        break;
      case 2:
        ;
        break;
      case 3:
        if (jj % 3 == 0) {
          is_put = true;
        }
        break;
      case 4:
        if (jj % 3 == 2) {
          is_put = true;
        }
        break;
      case 5:
        if (ii % 12 == 5) {
          if (jj != n - 1) {
            is_put = true;
          }
        }
        else {
          if (jj != 0) {
            is_put = true;
          }
        }
        break;
      default:
        assert(false);
        break;
      }

      if (is_put) {
        tmp_ps.push_back(P(i, j));
      }
    }
  }

  // シャッフル
  //std::shuffle(tmp_ps.begin(), tmp_ps.end(), engine);

  for (auto p : tmp_ps) {
    int i = p.first;
    int j = p.second;
    b[i][j] = 1;
    if (!can_reach_all_cells(5)) {
      b[i][j] = 0;
    }
    else {
      ps.emplace_back(i, j);
    }
  }

  return ps;
}

vector<P> init_1(const SimulateParam& param)
{
  vector<P> ps;

  attempt(ti, tj + 1, ps, 5);
  for (int k = 0; k < 100; k++) {
    int i1 = ti - 1 - k;
    int j1 = tj + 1 - k;
    attempt(i1, j1, ps, 5);
    int i2 = ti + 1 - k;
    int j2 = tj - k;
    attempt(i2, j2, ps, 5);
  }

  vector<P> tmp_ps;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (confirmed[i][j] == 1 || b[i][j] == 1 || (i == ti && j == tj) || (i == pi && j == pj)) {
        continue;
      }

      int ii = (i + param.slide_i + n) % n;
      int jj = (j + param.slode_j + n) % n;

      bool is_put = false;
      switch (ii % 6) {
      case 0:
        if (jj % 3 == 2) {
          is_put = true;
        }
        break;
      case 1:
        if (jj % 3 == 0) {
          is_put = true;
        }
        break;
      case 2:
        if (jj % 3 == 1) {
          is_put = true;
        }
        break;
      case 3:
        if (jj % 3 == 2) {
          is_put = true;
        }
        break;
      case 4:
        ;
        break;
      case 5:
        if (ii % 12 == 5) {
          if (jj != n - 1) {
            is_put = true;
          }
        }
        else {
          if (jj != 0) {
            is_put = true;
          }
        }
        break;
      default:
        assert(false);
        break;
      }

      if (is_put) {
        tmp_ps.push_back(P(i, j));
      }
    }
  }

  // シャッフル
  //std::shuffle(tmp_ps.begin(), tmp_ps.end(), engine);

  for (auto p : tmp_ps) {
    int i = p.first;
    int j = p.second;
    attempt(i, j, ps, 5);
  }

  return ps;
}

int method1(ofstream& ofs, const SimulateParam& param)
{
  while (true) {
    update_confirmed();
    if (pi == ti && pj == tj) {
      break;
    }

    vector<P> ps2;

    // まずは適当に思いついたルールで配置を作ってみる
    if (turn == 0) {
      if (param.init_method == 0) {
        ps2 = init_0(param);
      }
      else if (param.init_method == 1) {
        ps2 = init_1(param);
      }
      else {
        cerr << "Error: param.init_method = " << param.init_method << endl;
        ps2 = init_0(param);
      }
    }

    put_each_turn_output(ofs, ps2);

    turn++;
  }

  return turn;
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
  generate_q();

  int sim_loop = 0;
  while (true) {
    if (timer.get_elapsed_time() > TIME_LIMIT / 2) {
      break;
    }
    sim_loop++;

    SimulateParam param;
    param.init_method = rand_xorshift() % 2;
    param.slide_i = rand_xorshift() % 6;
    param.slode_j = rand_xorshift() % 6;

    initialize_simulate();
    int score = method1(ofs, param);
    if (score > bestScore) {
      bestScore = score;
      bestParam = param;
    }
  }

  cerr 
    << "sim_loop = " << sim_loop 
    << ", bestScore = " << bestScore 
    << ", init_method = " << bestParam.init_method
    << endl;

  is_simulate = false;
  setup_true_q();

  initialize_simulate();
  int score = method1(ofs, bestParam);

  if (ofs.is_open()) {
    ofs.close();
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
