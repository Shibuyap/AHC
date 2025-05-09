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

constexpr double TIME_LIMIT = 1.8;
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

const int n = 100;
const int L = 500000;

std::array<int, n> t;
std::array<int, n> orig_index;

int current_score;
std::array<int, n> a, b;

std::array<double, n> est_counts;
std::array<double, n> est_counts_cur;

int best_score;
std::array<int, n> best_a, best_b;

static void copy_to_best()
{
  best_score = current_score;
  best_a = a;
  best_b = b;
}

static void copy_from_best()
{
  current_score = best_score;
  a = best_a;
  b = best_b;
}

// 複数のケースを処理する際に、内部状態を初期化する関数
void init_state()
{
  current_score = 0;
  best_score = 0;
}

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
    for (int i = 0; i < n; ++i)cin >> t[i];
  }
  else {
    // ファイル入力
    int _n, _l;
    ifs >> _n >> _l;
    for (int i = 0; i < n; ++i)ifs >> t[i];
  }

  vector<P> vp;
  for (int i = 0; i < n; ++i) {
    vp.push_back(P(t[i], i));
  }
  sort(vp.begin(), vp.end());
  for (int i = 0; i < n; ++i) {
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
static int calc_score_exact()
{
  int visit_cnt[n] = {};
  int cur_v = 0;
  for (int i = 0; i < L; ++i) {
    visit_cnt[cur_v]++;
    if (visit_cnt[cur_v] % 2 == 1) {
      cur_v = a[cur_v];
    }
    else {
      cur_v = b[cur_v];
    }
  }

  int res = 1000000;
  for (int i = 0; i < n; ++i) {
    res -= abs(visit_cnt[i] - t[i]);
  }
  return res;
}

static int calc_score_sample(int sample_steps)
{
  int visit_cnt[n] = {};
  int cur_v = 0;
  for (int i = 0; i < sample_steps; ++i) {
    visit_cnt[cur_v]++;
    if (visit_cnt[cur_v] % 2 == 1) {
      cur_v = a[cur_v];
    }
    else {
      cur_v = b[cur_v];
    }
  }

  double res = 1000000;
  double scale_to_full = (double)L / sample_steps;
  for (int i = 0; i < n; ++i) {
    res -= abs(visit_cnt[i] * scale_to_full - t[i]);
  }
  return max(0, (int)round(res));
}

std::array<int, n> sorted_a;
std::array<int, n> sorted_b;
vector<P> cnt_sorted_vec(n);
std::array<int, n> index_map;
static int sort_and_estimate(int sample_steps)
{
  int visit_cnt[n] = {};
  int cur_v = 0;
  for (int i = 0; i < sample_steps; ++i) {
    visit_cnt[cur_v]++;
    if (visit_cnt[cur_v] % 2 == 1) {
      cur_v = a[cur_v];
    }
    else {
      cur_v = b[cur_v];
    }
    constexpr int EARLY_ABORT_STEP = 10000;
    if (i == EARLY_ABORT_STEP) {
      for (int j = 0; j < n; ++j) {
        if (t[j] > 100 && visit_cnt[j] == 0) {
          return 0;
        }
      }
    }
  }

  for (int i = 0; i < n; ++i) {
    cnt_sorted_vec[i].first = visit_cnt[i];
    cnt_sorted_vec[i].second = i;
  }
  sort(cnt_sorted_vec.begin(), cnt_sorted_vec.end());
  for (int i = 0; i < n; ++i) {
    index_map[cnt_sorted_vec[i].second] = i;
  }

  for (int i = 0; i < n; ++i) {
    sorted_a[i] = index_map[a[cnt_sorted_vec[i].second]];
    sorted_b[i] = index_map[b[cnt_sorted_vec[i].second]];
  }

  for (int i = 0; i < n; ++i) {
    if (cnt_sorted_vec[i].first > 0) {
      break;
    }
    if (t[i] > 100) {
      return 0;
    }
  }

  double scale_to_full = (double)L / sample_steps;
  for (int i = 0; i < n; ++i) {
    est_counts[i] = cnt_sorted_vec[i].first * scale_to_full;
  }

  double res = 1000000;
  for (int i = 0; i < n; ++i) {
    res -= abs(est_counts[i] - t[i]);
  }
  return max(0, (int)round(res));
}

std::array<double, n> est_counts_buf;
int dfs_stack[1000];
int dfs_tmp[1000];
static void reset_estimated_counts()
{
  for (int i = 0; i < n; ++i) {
    est_counts_buf[i] = est_counts_cur[i];
  }
}

static void update_estimated_counts(int s, double diff)
{
  dfs_stack[0] = s;
  int stack_tail = 1;
  for (int dfs = 0; dfs < 6; ++dfs) {
    int next_tail = 0;
    for (int j = 0; j < stack_tail; ++j) {
      int num = dfs_stack[j];
      est_counts_buf[num] += diff;
      dfs_tmp[next_tail] = a[num];
      next_tail++;
      dfs_tmp[next_tail] = b[num];
      next_tail++;
    }
    for (int j = 0; j < next_tail; ++j) {
      dfs_stack[j] = dfs_tmp[j];
    }
    stack_tail = next_tail;
    diff /= 2;
  }
}

static double calc_estimated_score()
{
  double res = 1000000;
  for (int i = 0; i < n; ++i) {
    res -= abs(est_counts_buf[i] - t[i]);
  }
  return max(0, (int)round(res));
}

// 解答を出力する関数
static void write_output(ofstream& ofs)
{
  int aaaa[n] = {};
  int bbbb[n] = {};
  for (int i = 0; i < n; ++i) {
    aaaa[orig_index[i]] = orig_index[a[i]];
    bbbb[orig_index[i]] = orig_index[b[i]];
  }

  if (mode == 0) {
    // 標準出力
    for (int i = 0; i < n; ++i)cout << aaaa[i] << ' ' << bbbb[i] << endl;
  }
  else {
    // ファイル出力
    for (int i = 0; i < n; ++i)ofs << aaaa[i] << ' ' << bbbb[i] << endl;
  }
}

static void build_initial_solution()
{
  // ランダムに初期解作成
  double now_time = get_elapsed_time();
  int loop_cnt_init = 0;
  while (true) {
    loop_cnt_init++;

    if (loop_cnt_init % 100 == 0) {
      now_time = get_elapsed_time();
      if (now_time > TIME_LIMIT / 10) break;
    }

    if (rand_u32() % 2 == 0) {
      for (int i = 0; i < n; ++i) {
        a[i] = rand_u32() % n;
        b[i] = rand_u32() % n;
      }
    }
    else {
      for (int i = 0; i < n; ++i) {
        if (rand_u32() % 2 == 0) {
          a[i] = rand_u32() % n;
          b[i] = (i + 1) % n;
        }
        else {
          a[i] = (i + 1) % n;
          b[i] = rand_u32() % n;
        }
      }
    }

    int tmpScore = sort_and_estimate(10000);
    if (tmpScore > current_score) {
      current_score = tmpScore;
      a = sorted_a;
      b = sorted_b;
      copy_to_best();
    }
  }

  copy_from_best();
  current_score = sort_and_estimate(25000);
  a = sorted_a;
  b = sorted_b;
  est_counts_cur = est_counts;
  copy_to_best();

  if (mode != 0 && mode != 3) {
    cout << loop_cnt_init << endl;
  }
}

// ハイパーパラメータ
struct AnnealingParams
{
  double start_temp;
  double end_temp;
  double multiple_value;
  int operation_thresholds[10];
};

static int get_omomi_rand_n()
{
  int cum_abs_diff[n] = {};
  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      cum_abs_diff[i] = cum_abs_diff[i - 1];
    }
    cum_abs_diff[i] += abs(t[i] - est_counts_cur[i]);
  }
  int rand_pick = rand_u32() % cum_abs_diff[n - 1];
  int node_id = lower_bound(cum_abs_diff, cum_abs_diff + n, rand_pick) - cum_abs_diff;
  return node_id;
}

static void simulated_annealing(AnnealingParams params)
{
  copy_to_best();

  build_initial_solution();

  int accept_stat[10][2];
  for (int i = 0; i < 10; ++i) {
    for (int j = 0; j < 2; ++j) {
      accept_stat[i][j] = 0;
    }
  }

  double now_time = get_elapsed_time();
  const double START_TEMP = params.start_temp;
  const double END_TEMP = params.end_temp;
  int anneal_iter = 0;
  while (true) {
    anneal_iter++;

    if (anneal_iter % 100 == 0) {
      now_time = get_elapsed_time();
      if (now_time > TIME_LIMIT) break;
    }

    int accept_move = 1;

    const double progress_ratio = now_time / TIME_LIMIT;
    const double temp = START_TEMP + (END_TEMP - START_TEMP) * progress_ratio;
    const int NEAR = 5;
    int op_rand = rand_u32() % params.operation_thresholds[9];
    int idx1, idx2, choice_flag;
    int backup_dst;
    if (op_rand < params.operation_thresholds[0]) {
      accept_stat[0][1]++;

      idx1 = get_omomi_rand_n();

      if (rand_u32() % 2 == 0) {
        idx2 = get_omomi_rand_n();
      }
      else {
        idx2 = idx1 + rand_u32() % (NEAR * 2 + 1) - NEAR;
        while (idx2 < 0 || idx2 >= n) {
          idx2 = idx1 + rand_u32() % (NEAR * 2 + 1) - NEAR;
        }
      }
      choice_flag = rand_u32() % 2;

      reset_estimated_counts();
      if (choice_flag == 0) {
        update_estimated_counts(a[idx1], -est_counts_cur[idx1] / 2.0);
        backup_dst = a[idx1];
        a[idx1] = idx2;
        update_estimated_counts(idx2, est_counts_cur[idx1] / 2.0);
      }
      else {
        update_estimated_counts(b[idx1], -est_counts_cur[idx1] / 2.0);
        backup_dst = b[idx1];
        b[idx1] = idx2;
        update_estimated_counts(idx2, est_counts_cur[idx1] / 2.0);
      }

      double earlyScore = calc_estimated_score();
      if (earlyScore < current_score - 10000) {
        accept_move = 0;
      }
    }
    else if (op_rand < params.operation_thresholds[2]) {
      accept_stat[2][1]++;
      idx1 = get_omomi_rand_n();
      if (rand_u32() % 10 == 0) {
        idx2 = get_omomi_rand_n();
      }
      else {
        idx2 = idx1 + rand_u32() % (NEAR * 2 + 1) - NEAR;
        while (idx2 < 0 || idx2 >= n) {
          idx2 = idx1 + rand_u32() % (NEAR * 2 + 1) - NEAR;
        }
      }
      choice_flag = rand_u32() % 2;

      reset_estimated_counts();
      if (choice_flag == 0) {
        update_estimated_counts(a[idx1], -est_counts_cur[idx1] / 2.0);
        update_estimated_counts(a[idx2], -est_counts_cur[idx2] / 2.0);
        swap(a[idx1], a[idx2]);
        update_estimated_counts(a[idx1], est_counts_cur[idx1] / 2.0);
        update_estimated_counts(a[idx2], est_counts_cur[idx2] / 2.0);
      }
      else {
        update_estimated_counts(b[idx1], -est_counts_cur[idx1] / 2.0);
        update_estimated_counts(b[idx2], -est_counts_cur[idx2] / 2.0);
        swap(b[idx1], b[idx2]);
        update_estimated_counts(b[idx1], est_counts_cur[idx1] / 2.0);
        update_estimated_counts(b[idx2], est_counts_cur[idx2] / 2.0);
      }

      double earlyScore = calc_estimated_score();
      if (earlyScore < current_score - 10000) {
        accept_move = 0;
      }
    }
    else if (op_rand < params.operation_thresholds[4]) {
      accept_stat[4][1]++;
      idx1 = get_omomi_rand_n();
      swap(a[idx1], b[idx1]);
    }
    else if (op_rand < params.operation_thresholds[5]) {
      accept_stat[5][1]++;
      idx1 = get_omomi_rand_n();
      if (rand_u32() % 10 == 0) {
        idx2 = get_omomi_rand_n();
      }
      else {
        idx2 = idx1 + rand_u32() % (NEAR * 2 + 1) - NEAR;
        while (idx2 < 0 || idx2 >= n) {
          idx2 = idx1 + rand_u32() % (NEAR * 2 + 1) - NEAR;
        }
      }

      reset_estimated_counts();
      update_estimated_counts(a[idx1], -est_counts_cur[idx1] / 2.0);
      update_estimated_counts(b[idx2], -est_counts_cur[idx2] / 2.0);
      swap(a[idx1], b[idx2]);
      update_estimated_counts(a[idx1], est_counts_cur[idx1] / 2.0);
      update_estimated_counts(b[idx2], est_counts_cur[idx2] / 2.0);
      double early_score = calc_estimated_score();
      if (early_score < current_score - 10000) {
        accept_move = 0;
      }
    }

    // スコア計算
    double trial_score = 0;
    if (accept_move) {
      trial_score = sort_and_estimate(25000);

      double diff_score_2 = (trial_score - current_score) * params.multiple_value;
      double prob_2 = exp(diff_score_2 / temp);
      accept_move = prob_2 > rand_unit();
    }

    if (accept_move) {
      // 採用
      current_score = trial_score;
      for (int i = 0; i < n; ++i) {
        a[i] = sorted_a[i];
        b[i] = sorted_b[i];
        est_counts_cur[i] = est_counts[i];
      }

      // Best解よりもいいか
      if (current_score > best_score) {
        copy_to_best();
      }
    }
    else {
      // 元に戻す
      if (op_rand < params.operation_thresholds[0]) {
        accept_stat[0][0]++;
        if (choice_flag == 0) {
          a[idx1] = backup_dst;
        }
        else {
          b[idx1] = backup_dst;
        }
      }
      else if (op_rand < params.operation_thresholds[2]) {
        accept_stat[2][0]++;
        if (choice_flag == 0) {
          swap(a[idx1], a[idx2]);
        }
        else {
          swap(b[idx1], b[idx2]);
        }
      }
      else if (op_rand < params.operation_thresholds[4]) {
        accept_stat[4][0]++;
        swap(a[idx1], b[idx1]);
      }
      else if (op_rand < params.operation_thresholds[5]) {
        accept_stat[5][0]++;
        swap(a[idx1], b[idx2]);
      }
    }
  }

  if (mode != 0 && mode != 3) {
    cout << anneal_iter << endl;
    for (int i = 0; i < 10; ++i) {
      cout << accept_stat[i][1] - accept_stat[i][0] << " / " << accept_stat[i][1] << endl;
    }
  }

  copy_from_best();
}

// 問題を解く関数
static ll solve_case(int problem_num, AnnealingParams hypers)
{
  reset_timer();

  // 複数ケース回すときに内部状態を初期値に戻す
  init_state();

  // 入力受け取り
  read_input(problem_num);

  // 出力ファイルストリームオープン
  ofstream ofs;
  open_output_file(problem_num, ofs);

  // 焼きなまし
  simulated_annealing(hypers);

  // 解答を出力
  write_output(ofs);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (mode != 0) {
    score = calc_score_exact();
  }
  return score;
}

int main()
{
  mode = 2;

  AnnealingParams params;
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
    solve_case(0, params);
  }
  else if (mode <= 2) {
    ll sum = 0;
    for (int i = 0; i < 15; ++i) {
      ll score = solve_case(i, params);
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
