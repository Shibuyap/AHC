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

const int DX[4] = { -1, 1, 0, 0 };
const int DY[4] = { 0, 0, -1, 1 };

const double TIME_LIMIT = 1.9;
int exec_mode;

const int MAX_N = 40;

int n;
int ti, tj;
int b[MAX_N][MAX_N];

int turn;
int confirmed[MAX_N][MAX_N];
int pi, pj;


// ローカルテスター用 入力
int current_q;
int qi[MAX_N * MAX_N - 1], qj[MAX_N * MAX_N - 1];
int ri, rj;

// 高速化用
int bfs_dp[MAX_N][MAX_N];
int bfs_visited[MAX_N][MAX_N];
int bfs_can_use[MAX_N][MAX_N];
int bfs_version;

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

    ifs.close();
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
  queue<P> que;
  que.push(P(pi, pj));
  bfs_visited[pi][pj] = bfs_version;
  while (!que.empty()) {
    P p = que.front();
    que.pop();
    int i = p.first;
    int j = p.second;
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
      que.push(P(ni, nj));
    }
  }
  return false;
}

bool can_reach_all_cells(int margin)
{
  bfs_version++;
  queue<P> que;
  que.push(P(pi, pj));
  bfs_visited[pi][pj] = bfs_version;
  while (!que.empty()) {
    P p = que.front();
    que.pop();
    int i = p.first;
    int j = p.second;
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
      que.push(P(ni, nj));
    }
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (b[i][j] == 0 && bfs_visited[i][j] != bfs_version) {
        margin--;
      }
    }
  }

  if (b[ti][tj] == 1) {
    margin = -1;
  }
  if (bfs_visited[ti][tj] != bfs_version) {
    margin = -1;
  }

  return margin >= 0;
}

int get_next_heros_move()
{
  //for(int i = 0; i < n; i++) {
  //  for(int j = 0; j < n; j++) {
  //    if(i == pi && j == pj) {
  //      cerr << 2;
  //    }
  //    else {
  //      cerr << b[i][j];
  //    }
  //  }
  //  cerr << endl;
  //}

  // 未確認のマスをすべて空きマスであると仮定してBFS
  bfs_version++;
  queue<P> que;
  que.push(P(pi, pj));
  bfs_visited[pi][pj] = bfs_version;
  bfs_dp[pi][pj] = 0;
  while (!que.empty()) {
    P p = que.front();
    que.pop();
    int i = p.first;
    int j = p.second;
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
      que.push(P(ni, nj));
    }
  }

  return bfs_dp[ri][rj] % 4;
}

void get_can_use()
{
  // 未確認のマスをすべて空きマスであると仮定して到達可能なすべての未確認マスをBFSで求める
  bfs_version++;
  queue<P> que;
  que.push(P(pi, pj));
  bfs_visited[pi][pj] = bfs_version;
  bfs_can_use[pi][pj] = 0;
  while (!que.empty()) {
    P p = que.front();
    que.pop();
    int i = p.first;
    int j = p.second;
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
      que.push(P(ni, nj));
    }
  }
}

void update_confirmed()
{
  if (exec_mode == 0) {
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
    // ローカルテスター用
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

      if (confirmed[ti][tj] == 1) {
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
}

void put_each_turn_output(ofstream& ofs, const vector<P>& ps)
{
  if (exec_mode == 0) {
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

vector<P> init_1()
{
  vector<P> ps;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i % 2 == 0 && j % 2 == 0 && confirmed[i][j] == 0 && b[i][j] == 0 && !(i == ti && j == tj)) {
        ps.emplace_back(i, j);
      }
    }
  }
  return ps;
}

vector<P> init_2()
{
  vector<P> ps;

  int loop_num = n * n / 4;
  for (int it = 0; it < loop_num; it++) {
    int i = rand_xorshift() % n;
    int j = rand_xorshift() % n;
    if (rand_xorshift() % 2 == 0) {
      i = rand_xorshift() % 3 + ti - 1;
      j = rand_xorshift() % 3 + tj - 1;
      i = max(i, 0);
      i = min(i, n - 1);
      j = max(j, 0);
      j = min(j, n - 1);
    }
    while (b[i][j] == 1 || (i == ti && j == tj) || (i == pi && j == pj)) {
      i = rand_xorshift() % n;
      j = rand_xorshift() % n;
    }

    b[i][j] = 1;
    if (!can_reach_goal()) {
      b[i][j] = 0;
    }
    else {
      ps.emplace_back(i, j);
    }
  }

  return ps;
}

vector<P> init_3()
{
  vector<P> ps;

  //int loop_num = n * n / 4;
  //for (int it = 0; it < loop_num; it++) {
  //  int i = rand_xorshift() % 5 + ti - 2;
  //  int j = rand_xorshift() % 5 + tj - 2;
  //  i = max(i, 0);
  //  i = min(i, n - 1);
  //  j = max(j, 0);
  //  j = min(j, n - 1);
  //  while (b[i][j] == 1 || (i == ti && j == tj) || (i == pi && j == pj)) {
  //    i = rand_xorshift() % 5 + ti - 2;
  //    j = rand_xorshift() % 5 + tj - 2;
  //    i = max(i, 0);
  //    i = min(i, n - 1);
  //    j = max(j, 0);
  //    j = min(j, n - 1);
  //  }

  //  b[i][j] = 1;
  //  if (!can_reach_all_cells()) {
  //    b[i][j] = 0;
  //  }
  //  else {
  //    ps.emplace_back(i, j);
  //  }
  //}

  //attempt(ti, tj + 1, ps);
  //for (int j = tj; j >= 0; j--) {
  //  attempt(ti - 1, j, ps);
  //  attempt(ti + 1, j, ps);
  //}

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
      bool is_put = false;
      if (i % 6 == 0) {
        if (j % 3 == 2) {
          is_put = true;
        }
      }
      else if (i % 6 == 1) {
        if (j % 3 == 0) {
          is_put = true;
        }
      }
      else if (i % 6 == 2) {
        ;
      }
      else if (i % 6 == 3) {
        if (j % 3 == 0) {
          is_put = true;
        }
      }
      else if (i % 6 == 4) {
        if (j % 3 == 2) {
          is_put = true;
        }
      }
      else if (i % 6 == 5) {
        if (i % 12 == 5) {
          if (j != n - 1) {
            is_put = true;
          }
        }
        else {
          if (j != 0) {
            is_put = true;
          }
        }
      }
      if (is_put) {
        tmp_ps.push_back(P(i, j));
      }
    }
  }

  // シャッフル
  int m = tmp_ps.size();
  int* arr = new int[m];
  for (int i = 0; i < m; i++) {
    arr[i] = i;
  }
  shuffle_array(arr, m);
  vector<P> shuffled_tmp_ps;
  for (int i = 0; i < m; i++) {
    shuffled_tmp_ps.push_back(tmp_ps[arr[i]]);
  }
  delete[] arr;
  tmp_ps = shuffled_tmp_ps;

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

vector<P> init_4()
{
  vector<P> ps;

  //int loop_num = n * n / 4;
  //int loop_num = 0;
  //for (int it = 0; it < loop_num; it++) {
  //  int i = rand_xorshift() % 5 + ti - 2;
  //  int j = rand_xorshift() % 5 + tj - 2;
  //  i = max(i, 0);
  //  i = min(i, n - 1);
  //  j = max(j, 0);
  //  j = min(j, n - 1);
  //  while (b[i][j] == 1 || (i == ti && j == tj) || (i == pi && j == pj)) {
  //    i = rand_xorshift() % 5 + ti - 2;
  //    j = rand_xorshift() % 5 + tj - 2;
  //    i = max(i, 0);
  //    i = min(i, n - 1);
  //    j = max(j, 0);
  //    j = min(j, n - 1);
  //  }

  //  attempt(i, j, ps, 5);
  //}

  attempt(ti, tj + 1, ps, 5);
  for (int k = 0; k < 100; k++) {
    int i1 = ti - 1 - k;
    int j1 = tj + 1 - k;
    attempt(i1, j1, ps, 5);
    int i2 = ti + 1 - k;
    int j2 = tj - k;
    attempt(i2, j2, ps, 5);
  }

  //attempt(ti, tj + 1, ps, 5);
  //for (int j = tj; j >= 0; j--) {
  //  attempt(ti - 1, j, ps, 5);
  //  attempt(ti + 1, j, ps, 5);
  //}

  vector<P> tmp_ps;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (confirmed[i][j] == 1 || b[i][j] == 1 || (i == ti && j == tj) || (i == pi && j == pj)) {
        continue;
      }
      bool is_put = false;
      if (i % 6 == 0) {
        if (j % 3 == 2) {
          is_put = true;
        }
      }
      else if (i % 6 == 1) {
        if (j % 3 == 0) {
          is_put = true;
        }
      }
      else if (i % 6 == 2) {
        if (j % 3 == 1) {
          is_put = true;
        }
      }
      else if (i % 6 == 3) {
        if (j % 3 == 2) {
          is_put = true;
        }
      }
      else if (i % 6 == 4) {
        ;
      }
      else if (i % 6 == 5) {
        if (i % 12 == 5) {
          if (j != n - 1) {
            is_put = true;
          }
        }
        else {
          if (j != 0) {
            is_put = true;
          }
        }
      }
      if (is_put) {
        tmp_ps.push_back(P(i, j));
      }
    }
  }

  // シャッフル
  int m = tmp_ps.size();
  int* arr = new int[m];
  for (int i = 0; i < m; i++) {
    arr[i] = i;
  }
  shuffle_array(arr, m);
  vector<P> shuffled_tmp_ps;
  for (int i = 0; i < m; i++) {
    shuffled_tmp_ps.push_back(tmp_ps[arr[i]]);
  }
  delete[] arr;
  //tmp_ps = shuffled_tmp_ps;

  for (auto p : tmp_ps) {
    int i = p.first;
    int j = p.second;
    attempt(i, j, ps, 5);
  }

  return ps;
}


int method1(ofstream& ofs)
{
  while (true) {
    update_confirmed();
    //cerr << "turn = " << turn << ", pi = " << pi << ", pj = " << pj << ", ri = " << ri << ", rj = " << rj << endl;
    if (pi == ti && pj == tj) {
      break;
    }

    vector<P> ps2;

    // まずは適当に思いついたルールで配置を作ってみる
    if (turn == 0) {
      //ps2 = init_1();
      //ps2 = init_2();
      ps2 = init_4();
      //ps2 = init_4();
    }

    put_each_turn_output(ofs, ps2);

    turn++;
  }

  return turn;
}

struct SimurateResult
{
  vector<P> ps;
  int score;
};

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
}

SimurateResult simurate(const vector<P>& ps)
{
  SimurateResult res;
  res.ps = ps;
  res.score = 0;
  return res;
}

int solve_case(int case_num)
{
  timer.start();

  initialize_state();

  input_data(case_num);

  ofstream ofs;
  open_ofs(case_num, ofs);

  int score = method1(ofs);

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
