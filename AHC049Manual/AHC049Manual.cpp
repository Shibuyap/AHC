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

using namespace std;

typedef long long int ll;
typedef pair<int, int> P;
typedef pair<P, P> PP;

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

const ll INF = 1001001001001001001LL;
const int INT_INF = 1001001001;

const int DX[4] = { -1, 0, 1, 0 };
const int DY[4] = { 0, -1, 0, 1 };
const char DIR_CHAR[4] = { 'U', 'L', 'D', 'R' };

const double TIME_LIMIT = 1.8;
int exec_mode;

const int MAX_N = 20;

struct Box
{
  int w;  // 重さ
  int d;  // 耐久力
  int initial_d;  // 初期耐久力
  bool exist;  // 箱がまだ存在するか

  Box() : w(0), d(0), initial_d(0), exist(false) {}
  Box(int weight, int durability) : w(weight), d(durability), initial_d(durability), exist(true) {}
};

struct Input
{
  int n;
  Box grid[MAX_N][MAX_N];  // グリッド上の箱

  void read(istream& is)
  {
    is >> n;

    // 重さを読み込み
    rep(i, n)
    {
      rep(j, n)
      {
        int w;
        is >> w;
        grid[i][j].w = w;
        grid[i][j].exist = (w > 0);
      }
    }

    // 耐久力を読み込み
    rep(i, n)
    {
      rep(j, n)
      {
        int d;
        is >> d;
        grid[i][j].d = d;
        grid[i][j].initial_d = d;
      }
    }
  }
};

struct State
{
  int x, y;  // 現在位置
  vector<Box> hand;  // 手に持っている箱（下から順）
  Box grid[MAX_N][MAX_N];  // グリッド上の箱
  int remaining_boxes;  // 残っている箱の数
  int move_count;  // 移動回数
  vector<char> actions;  // 実行した操作の履歴

  State() : x(0), y(0), remaining_boxes(0), move_count(0) {}

  void init(const Input& input)
  {
    x = 0;
    y = 0;
    hand.clear();
    remaining_boxes = 0;
    move_count = 0;
    actions.clear();

    rep(i, input.n)
    {
      rep(j, input.n)
      {
        grid[i][j] = input.grid[i][j];
        if (grid[i][j].exist && !(i == 0 && j == 0)) {
          remaining_boxes++;
        }
      }
    }
  }

  bool pickup()
  {
    if (!grid[x][y].exist) return false;
    if (x == 0 && y == 0) return false;  // 出入り口では拾えない

    hand.push_back(grid[x][y]);
    grid[x][y].exist = false;
    remaining_boxes--;
    actions.push_back('1');
    return true;
  }

  bool put()
  {
    if (hand.empty()) return false;
    if (grid[x][y].exist) return false;
    if (x == 0 && y == 0) return false;  // 出入り口には置けない

    grid[x][y] = hand.back();
    hand.pop_back();
    remaining_boxes++;
    actions.push_back('2');
    return true;
  }

  bool move(int dir)
  {
    int nx = x + DX[dir];
    int ny = y + DY[dir];

    if (nx < 0 || nx >= MAX_N || ny < 0 || ny >= MAX_N) return false;

    x = nx;
    y = ny;
    move_count++;

    // 手に持っている箱の耐久力を減少
    int total_weight = 0;
    for (int i = hand.size() - 1; i >= 0; i--) {
      hand[i].d -= total_weight;
      total_weight += hand[i].w;
    }

    actions.push_back(DIR_CHAR[dir]);

    // 出入り口に到着したら箱を全て運び出す
    if (x == 0 && y == 0 && !hand.empty()) {
      hand.clear();
    }

    return true;
  }

  ll calc_score(int n) const
  {
    if (remaining_boxes > 0) {
      return n * n - remaining_boxes;
    }
    else {
      return n * n + 2 * n * n * n - move_count;
    }
  }
};

Input input;
State current_state;
State best_state;

void store_best_state()
{
  best_state = current_state;
}

void restore_best_state()
{
  current_state = best_state;
}

bool is_out_of_range(int x, int y)
{
  if (x < 0 || MAX_N <= x || y < 0 || MAX_N <= y) return true;
  return false;
}

void initialize_state()
{
  current_state.init(input);
}

void input_data(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    input.read(cin);
  }
  else {
    // ファイル入力
    input.read(ifs);
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
  const vector<char>& actions = best_state.actions;

  if (exec_mode == 0) {
    // 標準出力
    for (char c : actions) {
      cout << c << endl;
    }
  }
  else {
    // ファイル出力
    for (char c : actions) {
      ofs << c << endl;
    }
  }
}

ll calculate_score()
{
  return best_state.calc_score(input.n);
}

void build_initial_solution()
{
  current_state.init(input);

  // 複数箱を同時運搬する貪欲解法
  while (current_state.remaining_boxes > 0) {
    // 一番重い箱
    int max_w = -1;
    int target_x = -1, target_y = -1;

    for (int i = 0; i < input.n; i++) {
      for (int j = 0; j < input.n; j++) {
        if (current_state.grid[i][j].exist) {
          if (current_state.grid[i][j].w > max_w) {
            max_w = current_state.grid[i][j].w;
            target_x = i;
            target_y = j;
          }
        }
      }
    }

    if (target_x == -1) break;

    // 目標位置への移動経路上で箱を拾う
    // まずx方向に移動
    while (current_state.x != target_x) {
      if (current_state.x < target_x) {
        current_state.move(2);  // Down
      }
      else {
        current_state.move(0);  // Up
      }
    }

    // 次にy方向に移動
    while (current_state.y != target_y) {
      if (current_state.y < target_y) {
        current_state.move(3);  // Right
      }
      else {
        current_state.move(1);  // Left
      }
    }

    // 目標位置の箱を拾う
    if (current_state.grid[current_state.x][current_state.y].exist) {
      current_state.pickup();
    }

    // 出入り口へ戻る（途中で箱を拾いながら）
    while (current_state.x != 0) {
      current_state.move(0);  // Up

      // 移動後の位置に箱があれば拾う
      if (current_state.grid[current_state.x][current_state.y].exist) {
        int total_weight = 0;
        for (const auto& box : current_state.hand) {
          total_weight += box.w;
        }

        int dist_to_home = current_state.x + current_state.y;
        if (dist_to_home * total_weight < current_state.grid[current_state.x][current_state.y].d) {
          current_state.pickup();
        }
      }
    }

    while (current_state.y != 0) {
      current_state.move(1);  // Left

      // 移動後の位置に箱があれば拾う
      if (current_state.grid[current_state.x][current_state.y].exist) {
        int total_weight = 0;
        for (const auto& box : current_state.hand) {
          total_weight += box.w;
        }

        int dist_to_home = current_state.x + current_state.y;
        if (dist_to_home * total_weight < current_state.grid[current_state.x][current_state.y].d) {
          current_state.pickup();
        }
      }
    }
  }

  best_state = current_state;
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
  store_best_state();

  double now_time = get_elapsed_time();
  const double START_TEMP = annealingParams.start_temperature[0];
  const double END_TEMP = annealingParams.end_temperature;

  int loop = 0;
  while (true) {
    loop++;

    if (loop % 100 == 0) {
      now_time = get_elapsed_time();
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
    ll tmp_score = current_state.calc_score(input.n);
    ll old_score = best_state.calc_score(input.n);

    // 焼きなましで採用判定
    double diff_score = (tmp_score - old_score) * annealingParams.score_scale;
    double prob = exp(diff_score / temp);
    if (prob > rand_01()) {
      // 採用

      // ベスト更新
      if (tmp_score > best_state.calc_score(input.n)) {
        store_best_state();
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

  restore_best_state();
}

ll solve_case(int case_num, AnnealingParams annealingParams)
{
  start_timer();

  initialize_state();

  input_data(case_num);

  ofstream ofs;
  open_ofs(case_num, ofs);

  build_initial_solution();

  // 焼きなまし実行
  // run_simulated_annealing(annealingParams);

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
  exec_mode = 2;  // 0: 標準入出力, 1-2: ファイル入出力

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
    for (int i = 0; i < 150; ++i) {
      ll score = solve_case(i, annealingParams);
      sum_score += score;
      if (exec_mode == 1) {
        cerr << score << endl;
      }
      else {
        cerr << "case = " << setw(3) << i << ", "
          << "score = " << setw(6) << score << ", "
          << "sum = " << setw(8) << sum_score << ", "
          << "time = " << setw(5) << get_elapsed_time() << ", "
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
      for (int i = 0; i < 150; ++i) {
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
