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
#define rep(i, n) for (int i = 0; i < (n); ++i)
#define srep(i, s, t) for (int i = s; i < t; ++i)
#define drep(i, n) for (int i = (n) - 1; i >= 0; --i)
#define dsrep(i, s, t) for (int i = (t) - 1; i >= s; --i)

using namespace std;

typedef long long int ll;
typedef pair<int, int> P;
typedef pair<P, P> PP;



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

const ll INF = 1001001001001001001LL;
const int INT_INF = 1001001001;

const int DX[5] = { -1, 0, 1, 0,0 };
const int DY[5] = { 0, -1, 0, 1,0 };
const char DC[5] = { 'U', 'L', 'D', 'R','S' };

const double TIME_LIMIT = 1.9;
int exec_mode;

const int n = 30;
const int m = 10;
const int k = 10;

int sx[n], sy[n];
int v[n][n - 1];
int h[n - 1][n];

int c[k][m];
vector<int> ans;

int f[n][n];

int current_score;

int best_score;

void store_best_score()
{
  best_score = current_score;
}

void restore_best_score()
{
  current_score = best_score;
}

void initialize_state()
{
  ans.clear();
  current_score = 0;
}

void input_data(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  int _n, _m, _k;
  if (!ifs.is_open()) {
    // 標準入力
    cin >> _n >> _m >> _k;
    rep(i, m)
    {
      cin >> sx[i] >> sy[i];
    }
    rep(i, n)
    {
      string s;
      cin >> s;
      rep(j, n - 1)
      {
        if (s[j] == '1') {
          v[i][j] = 1;
        }
        else {
          v[i][j] = 0;
        }

      }
    }
    rep(i, n - 1)
    {
      string s;
      cin >> s;
      rep(j, n)
      {
        if (s[j] == '1') {
          h[i][j] = 1;
        }
        else {
          h[i][j] = 0;
        }
      }
    }
  }
  else {
    // ファイル入力
    ifs >> _n >> _m >> _k;
    rep(i, m)
    {
      ifs >> sx[i] >> sy[i];
    }
    rep(i, n)
    {
      string s;
      ifs >> s;
      rep(j, n - 1)
      {
        if (s[j] == '1') {
          v[i][j] = 1;
        }
        else {
          v[i][j] = 0;
        }
      }
    }
    rep(i, n - 1)
    {
      string s;
      ifs >> s;
      rep(j, n)
      {
        if (s[j] == '1') {
          h[i][j] = 1;
        }
        else {
          h[i][j] = 0;
        }
      }
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

void output_data(ofstream& ofs)
{
  if (exec_mode == 0) {
    // 標準出力
    rep(i, k)
    {
      rep(j, m)
      {
        cout << DC[c[i][j]] << ' ';
      }
      cout << endl;
    }
    for (auto ansx : ans) {
      cout << ansx << endl;
    }
  }
  else {
    // ファイル出力
    rep(i, k)
    {
      rep(j, m)
      {
        ofs << DC[c[i][j]] << ' ';
      }
      ofs << endl;
    }
    for (auto ansx : ans) {
      ofs << ansx << endl;
    }
  }
}

ll calculate_score()
{
  ll res = 3 * n * n - ans.size();
  return res;
}

bool is_ng(int x, int y, int dir)
{
  if (dir == 0) {
    if (x == 0 || h[x - 1][y]) {
      return true;
    }
  }
  if (dir == 1) {
    if (y == 0 || v[x][y - 1]) {
      return true;
    }
  }
  if (dir == 2) {
    if (x == n - 1 || h[x][y]) {
      return true;
    }
  }
  if (dir == 3) {
    if (y == n - 1 || v[x][y]) {
      return true;
    }
  }
  return false;
}

int bs[m];
int x[m], y[m];
int cnt;
void push_button(int z)
{
  int is_move = 0;
  rep(j, m)
  {
    int dir = c[z][j];
    if (dir == 4) continue;
    if (is_ng(x[j], y[j], dir)) continue;
    is_move = 1;
    x[j] += DX[dir];
    y[j] += DY[dir];
    if (f[x[j]][y[j]] == 0) {
      f[x[j]][y[j]] = 1;
      cnt++;
    }
  }
  if (is_move) {
    ans.push_back(z);
  }
}

int visited[n][n];
int visited_version = 0;
int dp[n][n];
int dp2[n][n];
int dp3[n][n];
int decide_direction()
{
  int min_dir = -1;
  int min_need = INT_INF;
  queue<P> q;

  while (!q.empty()) q.pop();
  visited_version++;
  rep(jjj, m)
  {
    if (visited[x[jjj]][y[jjj]] == visited_version) continue;
    q.push({ x[jjj], y[jjj] });
    visited[x[jjj]][y[jjj]] = visited_version;
    dp[x[jjj]][y[jjj]] = 0;
    dp3[x[jjj]][y[jjj]] = jjj;
  }
  //rep(i, n)
  //{
  //  rep(j, n)
  //  {
  //    dp[i][j] = INT_INF;
  //  }
  //}
  int gx = -1;
  int gy = -1;
  while (!q.empty()) {
    P p = q.front();
    q.pop();
    int cx = p.first;
    int cy = p.second;
    int jjj = dp3[cx][cy];
    int fin = 0;
    srep(dd, 4, 8)
    {
      int d = c[dd][jjj];
      if (is_ng(cx, cy, d)) continue;
      int nx = cx + DX[d];
      int ny = cy + DY[d];
      if (visited[nx][ny] == visited_version) continue;
      visited[nx][ny] = visited_version;
      dp[nx][ny] = dp[cx][cy] + 1;
      dp2[nx][ny] = dd;
      dp3[nx][ny] = jjj;
      if (f[nx][ny] == 0) {
        if (dp[nx][ny] < min_need) {
          min_need = dp[nx][ny];
          gx = nx;
          gy = ny;
        }

        fin = 1;
        break;
      }
      q.push({ nx, ny });
    }
    if (fin) break;
  }

  // 経路復元
  int dir = -1;
  int argdir = -1;
  int jjj = dp3[gx][gy];
  while (gx != x[jjj] || gy != y[jjj]) {
    jjj = dp3[gx][gy];
    argdir = dp2[gx][gy];
    dir = c[argdir][jjj];
    int nx = gx - DX[dir];
    int ny = gy - DY[dir];
    gx = nx;
    gy = ny;
  }
  min_dir = argdir;

  return argdir;
}

int simulate(vector<int> tmp, int breakpoint = 9999)
{
  rep(i, n)
  {
    rep(j, n)
    {
      f[i][j] = 0;
    }
  }

  rep(i, m)
  {
    x[i] = sx[i];
    y[i] = sy[i];
    f[x[i]][y[i]] = 1;
  }

  cnt = m;
  rep(i, tmp.size())
  {
    int z = tmp[i];
    push_button(z);
    if (cnt == n * n)break;
    if (i == breakpoint)break;
  }

  while (cnt < n * n) {
    int dir = decide_direction();
    if (dir == -1) break;
    push_button(dir);
    //cout << x[0] << ' ' << y[0] << ' ' << cnt << endl;
    if (calculate_score() < best_score) {
      return 0;
    }
  }
  return calculate_score();
}

int build_initial_solution(int num1, int num2, int n3, int n4, int make_c, int n5, int n6, int n7)
{
  if (make_c == 1) {
    rep(i, m)
    {
      bs[i] = (n4 >> i) & 1;
    }
    rep(j, m)
    {
      rep(i, k)
      {
        if (bs[j] == 0) {
          c[i][j] = i % 4;
        }
        else {
          c[i][j] = (i + 1) % 4;
        }
      }
      srep(i, 8, 10)
      {
        c[i][j] = rand_xorshift() % 4;
      }
    }
  }


  rep(i, n)
  {
    rep(j, n)
    {
      f[i][j] = 0;
    }
  }


  rep(i, m)
  {
    x[i] = sx[i];
    y[i] = sy[i];
    f[x[i]][y[i]] = 1;
  }

  cnt = m;

  rep(i, n7)
  {
    if (i % 2 == 0) {
      push_button(8);
    }
    else {
      push_button(9);
    }
  }

  vector<vector<int>> vs = {
    {0,1,3,2},
    {2,1,3,0},
    {1,0,2,3},
    {3,0,2,1},
  };

  while (cnt < n * n) {
    // 上に5回移動
    rep(i, num1)
    {
      push_button(vs[n3][0]);
    }
    // 左にn-1回移動
    rep(i, n5)
    {
      push_button(vs[n3][1]);
    }
    if (n6 == 1) {
      push_button(vs[n3][3]);
    }
    // 蛇腹に5列移動
    rep(i, num2)
    {
      if (i % 2 == 0) {
        rep(j, n5)
        {
          push_button(vs[n3][2]);
        }
      }
      else {
        rep(j, n5)
        {
          push_button(vs[n3][1]);
        }
      }
      push_button(vs[n3][3]);
    }
    break;
  }

  while (cnt < n * n) {
    int dir = decide_direction();
    if (dir == -1) break;
    push_button(dir);
    //cout << x[0] << ' ' << y[0] << ' ' << cnt << endl;
    if (calculate_score() < best_score) {
      return 0;
    }
  }
  if (cnt != n * n) {
    cerr << "cnt != n*n" << endl;
  }

  return calculate_score();
}

struct AnnealingParams
{
  double start_temperature[10];
  double end_temperature;
  double score_scale;
  int operation_thresholds[10];
};

void run_simulated_annealing(AnnealingParams annealingParams)
{
  store_best_score();

  double now_time = timer.get_elapsed_time();
  const double START_TEMP = annealingParams.start_temperature[0];
  const double END_TEMP = annealingParams.end_temperature;

  int loop = 0;
  while (true) {
    loop++;

    if (loop % 100 == 0) {
      now_time = timer.get_elapsed_time();
      if (now_time > TIME_LIMIT) { break; }
    }

    double progress_ratio = now_time / TIME_LIMIT;
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progress_ratio;

    // 近傍解作成
    int ra_exec_mode = rand_xorshift() % annealingParams.operation_thresholds[1];
    int ra1, ra2, ra3, ra4, ra5;
    int keep1, keep2, keep3, keep4, keep5;

    if (ra_exec_mode < annealingParams.operation_thresholds[0]) {
      // 近傍操作1
    }
    else if (ra_exec_mode < annealingParams.operation_thresholds[1]) {
      // 近傍操作2
    }

    // スコア計算
    double tmp_score = calculate_score();

    // 焼きなましで採用判定
    double diff_score = (tmp_score - current_score) * annealingParams.score_scale;
    double prob = exp(diff_score / temp);
    if (prob > rand_01()) {
      // 採用
      current_score = tmp_score;

      // ベスト更新
      if (current_score > best_score) {
        store_best_score();
      }
    }
    else {
      // 元に戻す
      if (ra_exec_mode < annealingParams.operation_thresholds[0]) {
        // 近傍操作1 の巻き戻し
      }
      else if (ra_exec_mode < annealingParams.operation_thresholds[1]) {
        // 近傍操作2 の巻き戻し
      }
    }
  }

  if (exec_mode != 0 && exec_mode != 3) {
    cerr << loop << endl;
  }

  restore_best_score();
}

ll solve_case(int case_num, AnnealingParams annealingParams)
{
  timer.start();

  initialize_state();

  input_data(case_num);

  ofstream ofs;
  open_ofs(case_num, ofs);

  vector<int> best_ans;
  int best_c[m][m];
  best_score = 0;
  int best_num1 = -1;
  int best_num2 = -1;
  int best_n3 = -1;
  int best_n4 = -1;
  int best_n5 = -1;
  int best_n6 = -1;
  int best_n7 = -1;
  int loop = 0;
  while (true) {
    if (timer.get_elapsed_time() > TIME_LIMIT / 2) break;
    loop++;
    int num1 = rand_range(0, 15);
    int num2 = rand_range(0, 15);
    int n3 = rand_range(0, 3);
    int n4 = rand_xorshift() % (1 << m);
    int n5 = (n - 1) - rand_xorshift() % 1;
    int n6 = rand_xorshift() % 2;
    int n7 = rand_xorshift() % 2;
    ans.clear();
    int current_score = build_initial_solution(num1, num2, n3, n4, 1, n5, n6, n7);
    if (current_score > best_score) {
      best_score = current_score;
      rep(i, k)
      {
        rep(j, m)
        {
          best_c[i][j] = c[i][j];
        }
      }
      best_ans = ans;
      best_num1 = num1;
      best_num2 = num2;
      best_n3 = n3;
      best_n4 = n4;
      best_n5 = n5;
      best_n6 = n6;
      best_n7 = n7;
    }
  }
  //rep(num1, 10)
  //{
  //  rep(num2, 7)
  //  {
  //    ans.clear();
  //    current_score = build_initial_solution(num1,num2);
  //    if( current_score > best_score ) {
  //      best_score = current_score;
  //      best_ans = ans;
  //      best_num1 = num1;
  //      best_num2 = num2;
  //    }
  //  }
  //}
  ans = best_ans;
  rep(i, k)
  {
    rep(j, m)
    {
      c[i][j] = best_c[i][j];
    }
  }
  cerr << "loop = " << loop << endl;
  //cout << "best_num1 = " << best_num1 << ", best_num2 = " << best_num2 << ", best_score = " << best_score << endl;

  while (true) {
    if (timer.get_elapsed_time() > TIME_LIMIT) break;
    loop++;
    int ra = rand_xorshift() % 100;
    if (ra < 20) {
      int row1 = rand_range(4, 7);
      int row2 = rand_range(4, 7);
      while (row1 == row2) {
        row2 = rand_range(4, 7);
      }
      int col = rand_range(0, m - 1);
      int val = rand_xorshift() % 4;
      swap(c[row1][col], c[row2][col]);
      ans.clear();
      int current_score = build_initial_solution(best_num1, best_num2, best_n3, best_n4, 0, best_n5, best_n6, best_n7);
      if (current_score > best_score) {
        cerr << "loop = " << loop << ", best_score = " << current_score << endl;
        best_score = current_score;
        rep(i, k)
        {
          rep(j, m)
          {
            best_c[i][j] = c[i][j];
          }
        }
        best_ans = ans;
      }
      else {
        swap(c[row1][col], c[row2][col]);
      }
    }
    else if (ra < 40) {
      int num = rand_xorshift() % best_ans.size();
      vector<int> tmp = best_ans;
      tmp.erase(tmp.begin() + num);
      ans.clear();
      int current_score = simulate(tmp);
      //cerr << current_score << endl;
      if (current_score > best_score) {
        cerr << "loop = " << loop << ", best_score AAA = " << current_score << endl;
        best_score = current_score;
        rep(i, k)
        {
          rep(j, m)
          {
            best_c[i][j] = c[i][j];
          }
        }
        best_ans = ans;
      }
    }
    else if (ra < 100) {
      int num = rand_xorshift() % best_ans.size();
      vector<int> tmp = best_ans;
      int ope = rand_xorshift() % 8;
      tmp.insert(tmp.begin() + num, ope);
      ans.clear();
      int current_score = simulate(tmp, num);
      //cerr << current_score << endl;
      if (current_score > best_score) {
        cerr << "loop = " << loop << ", best_score BBB = " << current_score << endl;
        best_score = current_score;
        rep(i, k)
        {
          rep(j, m)
          {
            best_c[i][j] = c[i][j];
          }
        }
        best_ans = ans;
      }
    }
  }
  ans = best_ans;
  rep(i, k)
  {
    rep(j, m)
    {
      c[i][j] = best_c[i][j];
    }
  }
  current_score = best_score;

  cerr << best_num1 << ' ' << best_num2 << ' ' << best_n3 << ' ' << best_n4 << ' ' << best_n5 << ' ' << best_n6 << ' ' << best_n7 << endl;

  // 解答を出力
  output_data(ofs);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (exec_mode != 0) {
    score = calculate_score();
  }
  return score;
}

int main()
{
  exec_mode = 2;

  AnnealingParams annealingParams;
  annealingParams.start_temperature[0] = 2048.0;
  annealingParams.start_temperature[1] = 2048.0;
  annealingParams.start_temperature[2] = 2048.0;
  annealingParams.start_temperature[3] = 2048.0;
  annealingParams.start_temperature[4] = 2048.0;
  annealingParams.start_temperature[5] = 2048.0;
  annealingParams.start_temperature[6] = 2048.0;
  annealingParams.start_temperature[7] = 2048.0;
  annealingParams.start_temperature[8] = 2048.0;
  annealingParams.start_temperature[9] = 2048.0;
  annealingParams.end_temperature = 0.0;
  annealingParams.score_scale = 12345.0;
  annealingParams.operation_thresholds[0] = 100;
  annealingParams.operation_thresholds[1] = 200;
  annealingParams.operation_thresholds[2] = 300;
  annealingParams.operation_thresholds[3] = 400;
  annealingParams.operation_thresholds[4] = 500;
  annealingParams.operation_thresholds[5] = 600;
  annealingParams.operation_thresholds[6] = 700;
  annealingParams.operation_thresholds[7] = 800;
  annealingParams.operation_thresholds[8] = 900;
  annealingParams.operation_thresholds[9] = 1000;

  if (exec_mode == 0) {
    solve_case(0, annealingParams);
  }
  else if (exec_mode <= 2) {
    ll sum_score = 0;
    for (int i = 0; i < 15; ++i) {
      ll score = solve_case(i, annealingParams);
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
  else if (exec_mode == 3) {
    int loop_count = 0;
    AnnealingParams best_annealingParams;
    ll best_sum_score = 0;

    while (true) {
      AnnealingParams new_annealingParams;
      new_annealingParams.start_temperature[0] = pow(2.0, rand_01() * 20);
      new_annealingParams.end_temperature = 0.0;
      new_annealingParams.score_scale = pow(2.0, rand_01() * 20);
      new_annealingParams.operation_thresholds[0] = rand() % 101;

      ll sum_score = 0;
      for (int i = 0; i < 15; ++i) {
        ll score = solve_case(i, new_annealingParams);
        sum_score += score;

        // シード0が悪ければ打ち切り
        if (i == 0 && score < 0) {
          break;
        }
      }

      cerr << "loop_count = " << loop_count
        << ", sum_score = " << sum_score
        << ", start_temperature = " << new_annealingParams.start_temperature[0]
        << ", end_temperature = " << new_annealingParams.end_temperature
        << ", score_scale = " << new_annealingParams.score_scale
        << ", operation_thresholds = " << new_annealingParams.operation_thresholds[0]
        << endl;

      if (sum_score > best_sum_score) {
        best_sum_score = sum_score;
        best_annealingParams = new_annealingParams;
      }

      loop_count++;
    }
  }

  return 0;
}
