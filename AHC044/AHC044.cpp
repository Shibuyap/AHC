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
#define drep(i, n) for (int i = (n)-1; i >= 0; --i)
#define dsrep(i, s, t) for (int i = (t)-1; i >= s; --i)

using namespace std;

typedef long long int ll;
typedef pair<int, int> P;
typedef pair<P, P> PP;

static uint32_t rand_u32() {
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

static double rand_unit() {
  return (rand_u32() + 0.5) * (1.0 / UINT_MAX);
}

constexpr ll INF = 1001001001001001001;
constexpr int INT_INF = 1001001001;

constexpr int DX[4] = { -1, 0, 1, 0 };
constexpr int DY[4] = { 0, -1, 0, 1 };

constexpr double TIME_LIMIT = 1.8;
int mode;

std::chrono::steady_clock::time_point startTimeClock;

static void reset_timer() {
  startTimeClock = std::chrono::steady_clock::now();
}

static double get_elapsed_time() {
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

static void copy_to_best() {
  best_score = current_score;
  best_a = a;
  best_b = b;
}

static void copy_from_best() {
  current_score = best_score;
  a = best_a;
  b = best_b;
}

// 複数のケースを処理する際に、内部状態を初期化する関数
void init_state() {
  current_score = 0;
  best_score = 0;
}

// 入力を受け取る関数
void read_input(int problemNum) {
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    int nnn, lll;
    cin >> nnn >> lll;
    rep(i, n)cin >> t[i];
  }
  else {
    // ファイル入力
    int nnn, lll;
    ifs >> nnn >> lll;
    rep(i, n)ifs >> t[i];
  }

  vector<P> vp;
  rep(i, n) {
    vp.push_back(P(t[i], i));
  }
  sort(vp.begin(), vp.end());
  rep(i, n) {
    t[i] = vp[i].first;
    orig_index[i] = vp[i].second;
  }
}

// 出力ファイルストリームを開く関数
static void open_output_file(int probNum, ofstream& ofs) {
  if (mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << probNum << ".txt";
    ofs.open(oss.str());
  }
}

// スコアを計算する関数
static int calc_score_exact() {
  int cnt[n] = {};
  int now = 0;
  rep(i, L) {
    cnt[now]++;
    if (cnt[now] % 2 == 1) {
      now = a[now];
    }
    else {
      now = b[now];
    }
  }

  int res = 1000000;
  rep(i, n) {
    res -= abs(cnt[i] - t[i]);
  }
  return res;
}

static int calc_score_sample(int k) {
  int cnt[n] = {};
  int now = 0;
  rep(i, k) {
    cnt[now]++;
    if (cnt[now] % 2 == 1) {
      now = a[now];
    }
    else {
      now = b[now];
    }
  }

  double res = 1000000;
  double mul = (double)L / k;
  rep(i, n) {
    res -= abs(cnt[i] * mul - t[i]);
  }
  return max(0, (int)round(res));
}

std::array<int, n> sorted_a;
std::array<int, n> sorted_b;
vector<P> cnt_sorted_vec(n);
std::array<int, n> index_map;
static int sort_and_estimate(int k) {
  int cnt[n] = {};
  int now = 0;
  rep(i, k) {
    cnt[now]++;
    if (cnt[now] % 2 == 1) {
      now = a[now];
    }
    else {
      now = b[now];
    }
    if (i == 10000) {
      rep(j, n) {
        if (t[j] > 100 && cnt[j] == 0) {
          return 0;
        }
      }
    }
  }

  rep(i, n) {
    cnt_sorted_vec[i].first = cnt[i];
    cnt_sorted_vec[i].second = i;
  }
  sort(cnt_sorted_vec.begin(), cnt_sorted_vec.end());
  rep(i, n) {
    index_map[cnt_sorted_vec[i].second] = i;
  }

  rep(i, n) {
    sorted_a[i] = index_map[a[cnt_sorted_vec[i].second]];
    sorted_b[i] = index_map[b[cnt_sorted_vec[i].second]];
  }

  rep(i, n) {
    if (cnt_sorted_vec[i].first > 0) {
      break;
    }
    if (t[i] > 100) {
      return 0;
    }
  }

  double mul = (double)L / k;
  rep(i, n) {
    est_counts[i] = cnt_sorted_vec[i].first * mul;
  }

  double res = 1000000;
  rep(i, n) {
    res -= abs(est_counts[i] - t[i]);
  }
  return max(0, (int)round(res));
}

std::array<double, n> est_counts_buf;
int dfs_stack[1000];
int dfs_tmp[1000];
static void reset_estimated_counts() {
  rep(i, n) {
    est_counts_buf[i] = est_counts_cur[i];
  }
}

static void update_estimated_counts(int s, double diff) {
  dfs_stack[0] = s;
  int tail = 1;
  rep(dfs, 6) {
    int newTail = 0;
    rep(j, tail) {
      int num = dfs_stack[j];
      est_counts_buf[num] += diff;
      dfs_tmp[newTail] = a[num];
      newTail++;
      dfs_tmp[newTail] = b[num];
      newTail++;
    }
    rep(j, newTail) {
      dfs_stack[j] = dfs_tmp[j];
    }
    tail = newTail;
    diff /= 2;
  }
}

static double calc_estimated_score() {
  double res = 1000000;
  rep(i, n) {
    res -= abs(est_counts_buf[i] - t[i]);
  }
  return max(0, (int)round(res));
}

// 解答を出力する関数
static void write_output(ofstream& ofs) {
  int aaaa[n] = {};
  int bbbb[n] = {};
  rep(i, n) {
    aaaa[orig_index[i]] = orig_index[a[i]];
    bbbb[orig_index[i]] = orig_index[b[i]];
  }

  if (mode == 0) {
    // 標準出力
    rep(i, n)cout << aaaa[i] << ' ' << bbbb[i] << endl;
  }
  else {
    // ファイル出力
    rep(i, n)ofs << aaaa[i] << ' ' << bbbb[i] << endl;
  }
}

static void build_initial_solution() {
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
      rep(i, n) {
        a[i] = rand_u32() % n;
        b[i] = rand_u32() % n;
      }
    }
    else {
      rep(i, n) {
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
  int diff_sum[n] = {};
  rep(i, n)
  {
    if (i > 0) {
      diff_sum[i] = diff_sum[i - 1];
    }
    diff_sum[i] += abs(t[i] - est_counts_cur[i]);
  }
  int ra = rand_u32() % diff_sum[n - 1];
  int res = lower_bound(diff_sum, diff_sum + n, ra) - diff_sum;
  return res;
}

static void simulated_annealing(AnnealingParams hypers) {
  copy_to_best();

  build_initial_solution();

  int saitakuCount[10][2];
  rep(i, 10) {
    rep(j, 2) {
      saitakuCount[i][j] = 0;
    }
  }

  double now_time = get_elapsed_time();
  const double START_TEMP = hypers.start_temp;
  const double END_TEMP = hypers.end_temp;
  int loop_cnt_sa = 0;
  while (true) {
    loop_cnt_sa++;

    if (loop_cnt_sa % 100 == 0) {
      now_time = get_elapsed_time();
      if (now_time > TIME_LIMIT) break;
    }

    int ok = 1;

    double progressRatio = now_time / TIME_LIMIT;
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;
    int NEAR = 5;
    int raMode = rand_u32() % hypers.operation_thresholds[9];
    int ra1, ra2, ra3;
    int keep1;
    if (raMode < hypers.operation_thresholds[0]) {
      saitakuCount[0][1]++;

      //ra1 = rand_u32() % n;
      ra1 = get_omomi_rand_n();

      if (rand_u32() % 2 == 0) {
        //ra2 = rand_u32() % n;
        ra2 = get_omomi_rand_n();
      }
      else {
        ra2 = ra1 + rand_u32() % (NEAR * 2 + 1) - NEAR;
        while (ra2 < 0 || ra2 >= n) {
          ra2 = ra1 + rand_u32() % (NEAR * 2 + 1) - NEAR;
        }
      }
      ra3 = rand_u32() % 2;

      reset_estimated_counts();
      if (ra3 == 0) {
        update_estimated_counts(a[ra1], -est_counts_cur[ra1] / 2.0);
        keep1 = a[ra1];
        a[ra1] = ra2;
        update_estimated_counts(ra2, est_counts_cur[ra1] / 2.0);
      }
      else {
        update_estimated_counts(b[ra1], -est_counts_cur[ra1] / 2.0);
        keep1 = b[ra1];
        b[ra1] = ra2;
        update_estimated_counts(ra2, est_counts_cur[ra1] / 2.0);
      }

      double earlyScore = calc_estimated_score();
      if (earlyScore < current_score - 10000) {
        ok = 0;
      }
    }
    else if (raMode < hypers.operation_thresholds[2]) {
      saitakuCount[2][1]++;
      //ra1 = rand_u32() % n;
      ra1 = get_omomi_rand_n();
      if (rand_u32() % 10 == 0) {
        //ra2 = rand_u32() % n;
        ra2 = get_omomi_rand_n();
      }
      else {
        ra2 = ra1 + rand_u32() % (NEAR * 2 + 1) - NEAR;
        while (ra2 < 0 || ra2 >= n) {
          ra2 = ra1 + rand_u32() % (NEAR * 2 + 1) - NEAR;
        }
      }
      ra3 = rand_u32() % 2;

      reset_estimated_counts();
      if (ra3 == 0) {
        update_estimated_counts(a[ra1], -est_counts_cur[ra1] / 2.0);
        update_estimated_counts(a[ra2], -est_counts_cur[ra2] / 2.0);
        swap(a[ra1], a[ra2]);
        update_estimated_counts(a[ra1], est_counts_cur[ra1] / 2.0);
        update_estimated_counts(a[ra2], est_counts_cur[ra2] / 2.0);
      }
      else {
        update_estimated_counts(b[ra1], -est_counts_cur[ra1] / 2.0);
        update_estimated_counts(b[ra2], -est_counts_cur[ra2] / 2.0);
        swap(b[ra1], b[ra2]);
        update_estimated_counts(b[ra1], est_counts_cur[ra1] / 2.0);
        update_estimated_counts(b[ra2], est_counts_cur[ra2] / 2.0);
      }

      double earlyScore = calc_estimated_score();
      if (earlyScore < current_score - 10000) {
        ok = 0;
      }
    }
    else if (raMode < hypers.operation_thresholds[4]) {
      saitakuCount[4][1]++;
      //ra1 = rand_u32() % n;
      ra1 = get_omomi_rand_n();
      swap(a[ra1], b[ra1]);
    }
    else if (raMode < hypers.operation_thresholds[5]) {
      saitakuCount[5][1]++;
      //ra1 = rand_u32() % n;
      ra1 = get_omomi_rand_n();
      if (rand_u32() % 10 == 0) {
        //ra2 = rand_u32() % n;
        ra2 = get_omomi_rand_n();
      }
      else {
        ra2 = ra1 + rand_u32() % (NEAR * 2 + 1) - NEAR;
        while (ra2 < 0 || ra2 >= n) {
          ra2 = ra1 + rand_u32() % (NEAR * 2 + 1) - NEAR;
        }
      }

      reset_estimated_counts();
      update_estimated_counts(a[ra1], -est_counts_cur[ra1] / 2.0);
      update_estimated_counts(b[ra2], -est_counts_cur[ra2] / 2.0);
      swap(a[ra1], b[ra2]);
      update_estimated_counts(a[ra1], est_counts_cur[ra1] / 2.0);
      update_estimated_counts(b[ra2], est_counts_cur[ra2] / 2.0);
      double earlyScore = calc_estimated_score();
      if (earlyScore < current_score - 10000) {
        ok = 0;
      }
    }

    // スコア計算
    double tmpScore2 = 0;
    if (ok) {
      tmpScore2 = sort_and_estimate(25000);

      double diffScore2 = (tmpScore2 - current_score) * hypers.multiple_value;
      double prob2 = exp(diffScore2 / temp);
      ok = prob2 > rand_unit();
    }

    if (ok) {
      // 採用
      current_score = tmpScore2;
      rep(i, n) {
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
      if (raMode < hypers.operation_thresholds[0]) {
        saitakuCount[0][0]++;
        if (ra3 == 0) {
          a[ra1] = keep1;
        }
        else {
          b[ra1] = keep1;
        }
      }
      else if (raMode < hypers.operation_thresholds[2]) {
        saitakuCount[2][0]++;
        if (ra3 == 0) {
          swap(a[ra1], a[ra2]);
        }
        else {
          swap(b[ra1], b[ra2]);
        }
      }
      else if (raMode < hypers.operation_thresholds[4]) {
        saitakuCount[4][0]++;
        swap(a[ra1], b[ra1]);
      }
      else if (raMode < hypers.operation_thresholds[5]) {
        saitakuCount[5][0]++;
        swap(a[ra1], b[ra2]);
      }
    }
  }

  if (mode != 0 && mode != 3) {
    cout << loop_cnt_sa << endl;
    rep(i, 10) {
      cout << saitakuCount[i][1] - saitakuCount[i][0] << " / " << saitakuCount[i][1] << endl;
    }
  }

  copy_from_best();
}

// 問題を解く関数
static ll solve_case(int problem_num, AnnealingParams hypers) {
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

int main() {
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    rand_u32();
  }

  mode = 2;

  AnnealingParams HYPERS;
  HYPERS.start_temp = 2000000.0;
  HYPERS.end_temp = 0.0;
  HYPERS.multiple_value = 12345.0;
  HYPERS.operation_thresholds[0] = 200;
  HYPERS.operation_thresholds[1] = 200;
  HYPERS.operation_thresholds[2] = 400;
  HYPERS.operation_thresholds[3] = 400;
  HYPERS.operation_thresholds[4] = 440;
  HYPERS.operation_thresholds[5] = 700;
  HYPERS.operation_thresholds[6] = 700;
  HYPERS.operation_thresholds[7] = 700;
  HYPERS.operation_thresholds[8] = 700;
  HYPERS.operation_thresholds[9] = 700;

  if (mode == 0) {
    solve_case(0, HYPERS);
  }
  else if (mode <= 2) {
    ll sum = 0;
    srep(i, 0, 15) {
      ll score = solve_case(i, HYPERS);
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
