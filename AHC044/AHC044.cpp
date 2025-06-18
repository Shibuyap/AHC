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

typedef long long int ll;
typedef pair<int, int> P;
typedef pair<P, P> PP;

static uint32_t rand_u32()
{
  static uint32_t x = 123456789;
  static uint32_t y = 362436069;
  static uint32_t z = 521288629;
  static uint32_t w = 88675123;
  uint32_t t;

  t = x ^ (x << 11);
  x = y;
  y = z;
  z = w;
  return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
}

static double rand_unit()
{
  return (rand_u32() + 0.5) * (1.0 / UINT_MAX);
}

constexpr ll INF = 1001001001001001001;
constexpr int INT_INF = 1001001001;

constexpr int DX[4] = { -1, 0, 1, 0 };
constexpr int DY[4] = { 0, -1, 0, 1 };

constexpr double TL = 1.8;
int mode;

std::chrono::steady_clock::time_point startTimeClock;

static void reset_timer()
{
  startTimeClock = std::chrono::steady_clock::now();
}

static double get_elapsed_time()
{
  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - startTimeClock;
  return elapsed.count();
}

const int N = 100;
const int L = 500000;

std::array<int, N> t;
std::array<int, N> orig_index;

class State
{
public:
  std::array<int, N> a, b;
  std::array<int, N> best_a, best_b;
  std::array<double, N> est_counts, est_counts_cur;
  int current_score = 0, best_score = 0;

  void copy_to_best()
  {
    best_score = current_score;
    best_a = a;
    best_b = b;
  }

  void copy_from_best()
  {
    current_score = best_score;
    a = best_a;
    b = best_b;
  }

  void init_state()
  {
    current_score = 0;
    best_score = 0;
  }
};

// 入力を受け取る関数
void read_input(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    int _n, _l;
    cin >> _n >> _l;
    for (int i = 0; i < N; ++i)cin >> t[i];
  }
  else {
    // ファイル入力
    int _n, _l;
    ifs >> _n >> _l;
    for (int i = 0; i < N; ++i)ifs >> t[i];
  }

  vector<P> vp;
  for (int i = 0; i < N; ++i) {
    vp.push_back(P(t[i], i));
  }
  sort(vp.begin(), vp.end());
  for (int i = 0; i < N; ++i) {
    t[i] = vp[i].first;
    orig_index[i] = vp[i].second;
  }
}

// 出力ファイルストリームを開く関数
static void open_output_file(int case_num, ofstream& ofs)
{
  if (mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofs.open(oss.str());
  }
}

// スコアを計算する関数
static int calc_score(const State& st)
{
  int visit_cnt[N] = {};
  int cur_v = 0;
  for (int i = 0; i < L; ++i) {
    visit_cnt[cur_v]++;
    if (visit_cnt[cur_v] % 2 == 1) {
      cur_v = st.a[cur_v];
    }
    else {
      cur_v = st.b[cur_v];
    }
  }

  int res = 1000000;
  for (int i = 0; i < N; ++i) {
    res -= abs(visit_cnt[i] - t[i]);
  }
  return res;
}

static int calc_score_sample(int sample_steps, const State& st)
{
  int visit_cnt[N] = {};
  int cur_v = 0;
  for (int i = 0; i < sample_steps; ++i) {
    visit_cnt[cur_v]++;
    if (visit_cnt[cur_v] % 2 == 1) {
      cur_v = st.a[cur_v];
    }
    else {
      cur_v = st.b[cur_v];
    }
  }

  double res = 1000000;
  double scale_to_full = (double)L / sample_steps;
  for (int i = 0; i < N; ++i) {
    res -= abs(visit_cnt[i] * scale_to_full - t[i]);
  }
  return max(0, (int)round(res));
}

std::array<int, N> sorted_a;
std::array<int, N> sorted_b;
vector<P> cnt_sorted_vec(N);
std::array<int, N> index_map;
static int sort_and_estimate(int sample_steps, State& st)
{
  int visit_cnt[N] = {};
  int cur_v = 0;
  for (int i = 0; i < sample_steps; ++i) {
    visit_cnt[cur_v]++;
    if (visit_cnt[cur_v] % 2 == 1) {
      cur_v = st.a[cur_v];
    }
    else {
      cur_v = st.b[cur_v];
    }
    constexpr int EARLY_ABORT_STEP = 10000;
    if (i == EARLY_ABORT_STEP) {
      for (int j = 0; j < N; ++j) {
        if (t[j] > 100 && visit_cnt[j] == 0) {
          return 0;
        }
      }
    }
  }

  for (int i = 0; i < N; ++i) {
    cnt_sorted_vec[i].first = visit_cnt[i];
    cnt_sorted_vec[i].second = i;
  }
  sort(cnt_sorted_vec.begin(), cnt_sorted_vec.end());
  for (int i = 0; i < N; ++i) {
    index_map[cnt_sorted_vec[i].second] = i;
  }

  for (int i = 0; i < N; ++i) {
    sorted_a[i] = index_map[st.a[cnt_sorted_vec[i].second]];
    sorted_b[i] = index_map[st.b[cnt_sorted_vec[i].second]];
  }

  for (int i = 0; i < N; ++i) {
    if (cnt_sorted_vec[i].first > 0) {
      break;
    }
    if (t[i] > 100) {
      return 0;
    }
  }

  double scale_to_full = (double)L / sample_steps;
  for (int i = 0; i < N; ++i) {
    st.est_counts[i] = cnt_sorted_vec[i].first * scale_to_full;
  }

  double res = 1000000;
  for (int i = 0; i < N; ++i) {
    res -= abs(st.est_counts[i] - t[i]);
  }
  return max(0, (int)round(res));
}

std::array<double, N> est_counts_buf;
int dfs_stack[1000];
int dfs_tmp[1000];
static void reset_est_counts(const State& st)
{
  for (int i = 0; i < N; ++i) {
    est_counts_buf[i] = st.est_counts_cur[i];
  }
}

static void update_est_counts(int s, double diff, const State& st)
{
  dfs_stack[0] = s;
  int stack_tail = 1;
  for (int dfs = 0; dfs < 6; ++dfs) {
    int next_tail = 0;
    for (int j = 0; j < stack_tail; ++j) {
      int num = dfs_stack[j];
      est_counts_buf[num] += diff;
      dfs_tmp[next_tail] = st.a[num];
      next_tail++;
      dfs_tmp[next_tail] = st.b[num];
      next_tail++;
    }
    for (int j = 0; j < next_tail; ++j) {
      dfs_stack[j] = dfs_tmp[j];
    }
    stack_tail = next_tail;
    diff /= 2;
  }
}

static double calc_est_score()
{
  double res = 1000000;
  for (int i = 0; i < N; ++i) {
    res -= abs(est_counts_buf[i] - t[i]);
  }
  return max(0, (int)round(res));
}

// 解答を出力する関数
static void write_output(ofstream& ofs, const State& st)
{
  int a_out[N] = {};
  int b_out[N] = {};
  for (int i = 0; i < N; ++i) {
    a_out[orig_index[i]] = orig_index[st.a[i]];
    b_out[orig_index[i]] = orig_index[st.b[i]];
  }

  if (mode == 0) {
    // 標準出力
    for (int i = 0; i < N; ++i)cout << a_out[i] << ' ' << b_out[i] << endl;
  }
  else {
    // ファイル出力
    for (int i = 0; i < N; ++i)ofs << a_out[i] << ' ' << b_out[i] << endl;
  }
}

static void build_initial_solution(State& st)
{
  // ランダムに初期解作成
  double now_time = get_elapsed_time();
  int loop_cnt_init = 0;
  while (true) {
    loop_cnt_init++;

    if (loop_cnt_init % 100 == 0) {
      now_time = get_elapsed_time();
      if (now_time > TL / 10) { break; }
    }

    if (rand_u32() % 2 == 0) {
      for (int i = 0; i < N; ++i) {
        st.a[i] = rand_u32() % N;
        st.b[i] = rand_u32() % N;
      }
    }
    else {
      for (int i = 0; i < N; ++i) {
        if (rand_u32() % 2 == 0) {
          st.a[i] = rand_u32() % N;
          st.b[i] = (i + 1) % N;
        }
        else {
          st.a[i] = (i + 1) % N;
          st.b[i] = rand_u32() % N;
        }
      }
    }

    int tmp_score = sort_and_estimate(10000, st);
    if (tmp_score > st.current_score) {
      st.current_score = tmp_score;
      st.a = sorted_a;
      st.b = sorted_b;
      st.copy_to_best();
    }
  }

  st.copy_from_best();
  st.current_score = sort_and_estimate(25000, st);
  st.a = sorted_a;
  st.b = sorted_b;
  st.est_counts_cur = st.est_counts;
  st.copy_to_best();

  if (mode != 0 && mode != 3) {
    cout << loop_cnt_init << endl;
  }
}

// ハイパーパラメータ
struct Params
{
  double start_temp;
  double end_temp;
  double multiple_value;
  int operation_thresholds[10];
};

static int get_weighted_rand_node(const State& st)
{
  int cum_diff[N] = {};
  for (int i = 0; i < N; ++i) {
    if (i > 0) {
      cum_diff[i] = cum_diff[i - 1];
    }
    cum_diff[i] += abs(t[i] - st.est_counts_cur[i]);
  }
  int rand_pick = rand_u32() % cum_diff[N - 1];
  int node_id = lower_bound(cum_diff, cum_diff + N, rand_pick) - cum_diff;
  return node_id;
}

static void simulated_annealing(const Params& params, State& st)
{
  st.copy_to_best();

  build_initial_solution(st);

  int accept_stat[10][2];
  for (int i = 0; i < 10; ++i) {
    for (int j = 0; j < 2; ++j) {
      accept_stat[i][j] = 0;
    }
  }

  double now_time = get_elapsed_time();
  const double START_TEMP = params.start_temp;
  const double END_TEMP = params.end_temp;
  int iter = 0;
  while (true) {
    iter++;

    if (iter % 100 == 0) {
      now_time = get_elapsed_time();
      if (now_time > TL) { break; }
    }

    int accept_move = 1;

    const double progress_ratio = now_time / TL;
    const double temp = START_TEMP + (END_TEMP - START_TEMP) * progress_ratio;
    const int NEAR = 5;
    int r_op = rand_u32() % params.operation_thresholds[9];
    int i1, i2, choice_flag;
    int backup;
    if (r_op < params.operation_thresholds[0]) {
      accept_stat[0][1]++;

      i1 = get_weighted_rand_node(st);

      if (rand_u32() % 2 == 0) {
        i2 = get_weighted_rand_node(st);
      }
      else {
        i2 = i1 + rand_u32() % (NEAR * 2 + 1) - NEAR;
        while (i2 < 0 || i2 >= N) {
          i2 = i1 + rand_u32() % (NEAR * 2 + 1) - NEAR;
        }
      }
      choice_flag = rand_u32() % 2;

      reset_est_counts(st);
      if (choice_flag == 0) {
        update_est_counts(st.a[i1], -st.est_counts_cur[i1] / 2.0, st);
        backup = st.a[i1];
        st.a[i1] = i2;
        update_est_counts(i2, st.est_counts_cur[i1] / 2.0, st);
      }
      else {
        update_est_counts(st.b[i1], -st.est_counts_cur[i1] / 2.0, st);
        backup = st.b[i1];
        st.b[i1] = i2;
        update_est_counts(i2, st.est_counts_cur[i1] / 2.0, st);
      }

      double early_est_score = calc_est_score();
      if (early_est_score < st.current_score - 10000) {
        accept_move = 0;
      }
    }
    else if (r_op < params.operation_thresholds[2]) {
      accept_stat[2][1]++;
      i1 = get_weighted_rand_node(st);
      if (rand_u32() % 10 == 0) {
        i2 = get_weighted_rand_node(st);
      }
      else {
        i2 = i1 + rand_u32() % (NEAR * 2 + 1) - NEAR;
        while (i2 < 0 || i2 >= N) {
          i2 = i1 + rand_u32() % (NEAR * 2 + 1) - NEAR;
        }
      }
      choice_flag = rand_u32() % 2;

      reset_est_counts(st);
      if (choice_flag == 0) {
        update_est_counts(st.a[i1], -st.est_counts_cur[i1] / 2.0, st);
        update_est_counts(st.a[i2], -st.est_counts_cur[i2] / 2.0, st);
        swap(st.a[i1], st.a[i2]);
        update_est_counts(st.a[i1], st.est_counts_cur[i1] / 2.0, st);
        update_est_counts(st.a[i2], st.est_counts_cur[i2] / 2.0, st);
      }
      else {
        update_est_counts(st.b[i1], -st.est_counts_cur[i1] / 2.0, st);
        update_est_counts(st.b[i2], -st.est_counts_cur[i2] / 2.0, st);
        swap(st.b[i1], st.b[i2]);
        update_est_counts(st.b[i1], st.est_counts_cur[i1] / 2.0, st);
        update_est_counts(st.b[i2], st.est_counts_cur[i2] / 2.0, st);
      }

      double early_est_score = calc_est_score();
      if (early_est_score < st.current_score - 10000) {
        accept_move = 0;
      }
    }
    else if (r_op < params.operation_thresholds[4]) {
      accept_stat[4][1]++;
      i1 = get_weighted_rand_node(st);
      swap(st.a[i1], st.b[i1]);
    }
    else if (r_op < params.operation_thresholds[5]) {
      accept_stat[5][1]++;
      i1 = get_weighted_rand_node(st);
      if (rand_u32() % 10 == 0) {
        i2 = get_weighted_rand_node(st);
      }
      else {
        i2 = i1 + rand_u32() % (NEAR * 2 + 1) - NEAR;
        while (i2 < 0 || i2 >= N) {
          i2 = i1 + rand_u32() % (NEAR * 2 + 1) - NEAR;
        }
      }

      reset_est_counts(st);
      update_est_counts(st.a[i1], -st.est_counts_cur[i1] / 2.0, st);
      update_est_counts(st.b[i2], -st.est_counts_cur[i2] / 2.0, st);
      swap(st.a[i1], st.b[i2]);
      update_est_counts(st.a[i1], st.est_counts_cur[i1] / 2.0, st);
      update_est_counts(st.b[i2], st.est_counts_cur[i2] / 2.0, st);
      double early_est_score = calc_est_score();
      if (early_est_score < st.current_score - 10000) {
        accept_move = 0;
      }
    }

    // スコア計算
    double tmp_score = 0;
    if (accept_move) {
      tmp_score = sort_and_estimate(25000, st);

      double diff_score = (tmp_score - st.current_score) * params.multiple_value;
      double prob = exp(diff_score / temp);
      accept_move = prob > rand_unit();
    }

    if (accept_move) {
      // 採用
      st.current_score = tmp_score;
      for (int i = 0; i < N; ++i) {
        st.a[i] = sorted_a[i];
        st.b[i] = sorted_b[i];
        st.est_counts_cur[i] = st.est_counts[i];
      }

      // Best解よりもいいか
      if (st.current_score > st.best_score) {
        st.copy_to_best();
      }
    }
    else {
      // 元に戻す
      if (r_op < params.operation_thresholds[0]) {
        accept_stat[0][0]++;
        if (choice_flag == 0) {
          st.a[i1] = backup;
        }
        else {
          st.b[i1] = backup;
        }
      }
      else if (r_op < params.operation_thresholds[2]) {
        accept_stat[2][0]++;
        if (choice_flag == 0) {
          swap(st.a[i1], st.a[i2]);
        }
        else {
          swap(st.b[i1], st.b[i2]);
        }
      }
      else if (r_op < params.operation_thresholds[4]) {
        accept_stat[4][0]++;
        swap(st.a[i1], st.b[i1]);
      }
      else if (r_op < params.operation_thresholds[5]) {
        accept_stat[5][0]++;
        swap(st.a[i1], st.b[i2]);
      }
    }
  }

  if (mode != 0 && mode != 3) {
    cout << iter << endl;
    for (int i = 0; i < 10; ++i) {
      cout << accept_stat[i][1] - accept_stat[i][0] << " / " << accept_stat[i][1] << endl;
    }
  }

  st.copy_from_best();
}

// 問題を解く関数
static ll solve(int problem_num, Params params)
{
  reset_timer();

  State st;
  st.init_state();

  // 入力受け取り
  read_input(problem_num);

  // 出力ファイルストリームオープン
  ofstream ofs;
  open_output_file(problem_num, ofs);

  // 焼きなまし
  simulated_annealing(params, st);

  // 解答を出力
  write_output(ofs, st);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (mode != 0) {
    score = calc_score(st);
  }
  return score;
}

int main()
{
  mode = 2;

  Params params;
  params.start_temp = 2000000.0;
  params.end_temp = 0.0;
  params.multiple_value = 12345.0;
  params.operation_thresholds[0] = 200;
  params.operation_thresholds[1] = 200;
  params.operation_thresholds[2] = 400;
  params.operation_thresholds[3] = 400;
  params.operation_thresholds[4] = 440;
  params.operation_thresholds[5] = 700;
  params.operation_thresholds[6] = 700;
  params.operation_thresholds[7] = 700;
  params.operation_thresholds[8] = 700;
  params.operation_thresholds[9] = 700;

  if (mode == 0) {
    solve(0, params);
  }
  else if (mode <= 2) {
    ll sum = 0;
    for (int i = 0; i < 15; ++i) {
      ll score = solve(i, params);
      sum += score;
      if (mode == 1) {
        cout << score << endl;
      }
      else {
        cout << "num = " << setw(2) << i << ", ";
        cout << "score = " << setw(4) << score << ", ";
        cout << "sum = " << setw(5) << sum << ", ";
        cout << "time = " << setw(5) << get_elapsed_time() << ", ";
        cout << endl;
      }
    }
  }

  return 0;
}
