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

const double TIME_LIMIT = 1.8;
int exec_mode;

const int n = 36;
const int m = 12;
const int L = 1000000;

string target_strings[n];
vector<int> target_sequences[n];
int preferences[n];

int state_chars[m];

int transition_probs[m][m];
int current_score;

int best_transition_probs[m][m];
int best_score;

void store_best_score()
{
  best_score = current_score;
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m; ++j) {
      best_transition_probs[i][j] = transition_probs[i][j];
    }
  }
}

void restore_best_score()
{
  current_score = best_score;
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m; ++j) {
      transition_probs[i][j] = best_transition_probs[i][j];
    }
  }
}

void initialize_state()
{
  current_score = 0;
}

void input_data(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (exec_mode == 0 || !ifs.is_open()) {
    // 標準入力
    int _n, _m, _L;
    cin >> _n >> _m >> _L;
    for (int i = 0; i < n; ++i) {
      cin >> target_strings[i];
      cin >> preferences[i];
    }
  }
  else {
    // ファイル入力
    int _n, _m, _L;
    ifs >> _n >> _m >> _L;
    for (int i = 0; i < n; ++i) {
      ifs >> target_strings[i];
      ifs >> preferences[i];
    }
  }

  for (int i = 0; i < n; ++i) {
    vector<int> v;
    for (int j = 0; j < target_strings[i].size(); ++j) {
      if (target_strings[i][j] == ' ') { continue; }
      v.push_back(target_strings[i][j] - 'a');
    }
    target_sequences[i] = v;
  }

  vector<pair<int, vector<int>>> sorted_preferences;
  for (int i = 0; i < n; ++i) {
    sorted_preferences.push_back(make_pair(preferences[i], target_sequences[i]));
  }
  sort(sorted_preferences.begin(), sorted_preferences.end());
  for (int i = 0; i < n; ++i) {
    preferences[i] = sorted_preferences[i].first;
    target_sequences[i] = sorted_preferences[i].second;
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

// a[][] を 0〜100 の整数として持っている前提
inline void build_transition(double P[12][12])
{
  constexpr double kInv = 1.0 / 100.0;
  for (int i = 0; i < 12; ++i)
    for (int j = 0; j < 12; ++j)
      P[i][j] = transition_probs[i][j] * kInv;        // row-stochastic
}

inline void stationary_dist(const double P[12][12],
  double pi[12],          // ←出力 (正規化済み)
  int max_iter = 20)
{
  // 初期ベクトルは一様で十分
  double v[12], nv[12];
  for (int i = 0; i < 12; ++i) v[i] = 1.0 / 12.0;

  for (int it = 0; it < max_iter; ++it) {
    // nv = v * P
    for (int j = 0; j < 12; ++j) nv[j] = 0.0;
    for (int i = 0; i < 12; ++i) {
      const double vi = v[i];
      for (int j = 0; j < 12; ++j) nv[j] += vi * P[i][j];
    }
    // 収束判定 (L1 ノルム)
    double diff = 0.0;
    for (int j = 0; j < 12; ++j) diff += fabs(nv[j] - v[j]);
    if (diff < 1e-12) { break; }
    memcpy(v, nv, sizeof(v));
  }
  memcpy(pi, v, sizeof(double) * 12);
}

// 旧コード互換：visited[i] = π_i × 1e6
inline void build_visited(int visited[12], const double pi[12])
{
  constexpr double SCALE = 1'000'000.0;   // 旧コードと同じスケール
  for (int i = 0; i < 12; ++i)
    visited[i] = static_cast<int>(pi[i] * SCALE + 0.5);
}

inline long double prob_one_or_more(long double k)
{
  constexpr long double ONE_MILLION = 1'000'000.0L;
  if (k <= 0.0L)                  return 0.0L;   // 起こらない
  if (k >= ONE_MILLION - 1)       return 1.0L;   // ほぼ確実に起こる

  long double p = k / ONE_MILLION;               // 1 試行あたり成功確率
  // Q = 1 - (1-p)^1e6  = 1 - exp(1e6 * log(1-p))
  long double ln_term = ONE_MILLION * log1pl(-p);
  return 1.0L - expl(ln_term);
}

int evaluate_score(bool all = false)
{
  double P[12][12];
  double pi[12];
  int    visited[12];

  build_transition(P);
  stationary_dist(P, pi);     // ←高速・決定的
  build_visited(visited, pi); // ←従来と同じスケール

  double res = 0;

  for (int i = 0; i < n; ++i) {
    double dp[2] = {};
    double dp2[2] = {};

    for (int j = 0; j < target_sequences[i].size(); ++j) {
      int num = target_sequences[i][j];
      if (j == 0) {
        dp[0] = visited[num];
        dp[1] = visited[num + 6];
      }
      else {
        dp2[0] = 0;
        dp2[1] = 0;

        dp2[0] += dp[0] * transition_probs[target_sequences[i][j - 1]][target_sequences[i][j]] / 100.0;
        dp2[1] += dp[0] * transition_probs[target_sequences[i][j - 1]][target_sequences[i][j] + 6] / 100.0;
        dp2[0] += dp[1] * transition_probs[target_sequences[i][j - 1] + 6][target_sequences[i][j]] / 100.0;
        dp2[1] += dp[1] * transition_probs[target_sequences[i][j - 1] + 6][target_sequences[i][j] + 6] / 100.0;

        dp[0] = dp2[0];
        dp[1] = dp2[1];
      }
    }

    double appearance_prob = prob_one_or_more(dp[0] + dp[1]);
    res += appearance_prob * preferences[i];
    if (all) {
      cout << i << " " << dp[0] + dp[1] << " " << appearance_prob << " " << preferences[i] << endl;
    }
  }

  return round(res);
}

void output_data(ofstream& ofs)
{
  if (exec_mode == 0) {
    // 標準出力
    for (int i = 0; i < m; ++i) {
      cout << (char)('a' + state_chars[i]);
      for (int j = 0; j < m; ++j) {
        cout << ' ' << transition_probs[i][j];
      }
      cout << endl;
    }
  }
  else {
    // ファイル出力
    for (int i = 0; i < m; ++i) {
      ofs << (char)('a' + state_chars[i]);
      for (int j = 0; j < m; ++j) {
        ofs << ' ' << transition_probs[i][j];
      }
      ofs << endl;
    }
  }
}

vector<int> normalize_to_100(const vector<double>& v)
{
  double sum = 0;
  for (int i = 0; i < 12; ++i) {
    sum += v[i];
  }
  vector<int> res(12);
  if (sum < 0.1) {
    for (int i = 0; i < 100; ++i) {
      res[i % 12]++;
    }
    return res;
  }

  for (int i = 0; i < 12; ++i) {
    res[i] = round(v[i] / sum * 100.0);
  }
  int sum2 = 0;
  for (int i = 0; i < 12; ++i) {
    sum2 += res[i];
  }
  if (sum2 > 100) {
    int diff = sum2 - 100;
    while (diff > 0) {
      int x = rand_xorshift() % 12;
      if (res[x] > 0) {
        res[x]--;
        diff--;
      }
    }
  }
  else if (sum2 < 100) {
    int diff = 100 - sum2;
    while (diff > 0) {
      int x = rand_xorshift() % 12;
      res[x]++;
      diff--;
    }
  }
  return res;
}

void build_seed_solution(double time_limit)
{
  int max_transition_probs[m][m];
  int max_score = -1;

  int selected_indices[n];
  int state_sequence[20];

  vector<double> cnt_sum[12];
  for (int i = 0; i < 12; ++i) {
    cnt_sum[i] = vector<double>(12);
  }

  for (int i = 0; i < m; ++i) {
    state_chars[i] = i % 6;
  }

  int iter = 0;
  while (true) {
    if (get_elapsed_time() > time_limit) {
      break;
    }
    iter++;

    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < m; ++j) {
        transition_probs[i][j] = 0;
      }
    }

    int num_selected = rand_xorshift() % 10 + 2;
    int sz = 0;
    for (int i = 0; i < num_selected; ++i) {
      int string_index = n - 1 - i;
      if (rand_xorshift() % 100 < 90) {
        selected_indices[sz] = string_index;
        sz++;
      }
    }

    for (int l = 0; l < sz; ++l) {
      int string_index = selected_indices[l];
      for (int i = 0; i < target_sequences[string_index].size(); ++i) {
        state_sequence[i] = target_sequences[string_index][i] + rand_xorshift() % 2 * 6;
      }
      for (int i = 0; i < target_sequences[string_index].size() - 1; ++i) {
        int x = state_sequence[i];
        int y = state_sequence[i + 1];
        cnt_sum[x][y] += preferences[string_index];
      }
    }

    for (int i = 0; i < 12; ++i) {
      auto v = normalize_to_100(cnt_sum[i]);
      for (int j = 0; j < 12; ++j) {
        transition_probs[i][j] = v[j];
      }
      for (int j = 0; j < 12; ++j) {
        cnt_sum[i][j] = 0;
      }
    }

    int score = evaluate_score();
    if (score > max_score) {
      max_score = score;
      for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
          max_transition_probs[i][j] = transition_probs[i][j];
        }
      }
    }
  }

  cerr << iter << ' ' << max_score << ' ' << get_elapsed_time() << endl;

  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m; ++j) {
      transition_probs[i][j] = max_transition_probs[i][j];
    }
  }
}

struct AnnealingParams
{
  double start_temperature[10];
  double end_temperature;
  double score_scale;
  int operation_thresholds[10];
};

void run_simulated_annealing(AnnealingParams params, int mode, double time_limit)
{
  current_score = evaluate_score();
  store_best_score();

  double start_time = get_elapsed_time();
  double now_time = get_elapsed_time();
  const double START_TEMP = params.start_temperature[0];
  const double END_TEMP = params.end_temperature;

  int loop = 0;
  while (true) {
    loop++;

    if (loop % 100 == 0) {
      now_time = get_elapsed_time();
      if (now_time > time_limit) { break; }
    }

    //if (current_score * 1.2 < best_score) {
    //  restore_best_score();
    //  continue;
    //}

    double progress_ratio = (now_time - start_time) / (time_limit - start_time);
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progress_ratio;

    // 近傍解作成
    int operation_type = rand_xorshift() % params.operation_thresholds[2];
    int ra1, ra2, ra3, ra4, ra5;
    int keep1, keep2, keep3, keep4, keep5;

    if (operation_type < params.operation_thresholds[0]) {
      // 近傍操作1
      ra5 = rand_xorshift() % 15 + 1;
      ra1 = rand_xorshift() % 2;
      ra2 = rand_xorshift() % 6;
      ra3 = rand_xorshift() % 6;
      int ng = 0;
      while (transition_probs[ra2 + ra1 * 6][ra3 + ra1 * 6] < ra5) {
        ra2 = rand_xorshift() % 6;
        ra3 = rand_xorshift() % 6;
        ng++;
        if (ng > 100)break;
      }
      if (ng > 100)continue;
      ra4 = rand_xorshift() % 6;
      while (ra3 == ra4) {
        ra4 = rand_xorshift() % 6;
      }

      transition_probs[ra2 + ra1 * 6][ra3 + ra1 * 6] -= ra5;
      transition_probs[ra2 + ra1 * 6][ra4 + ra1 * 6] += ra5;
    }
    else if (operation_type < params.operation_thresholds[1]) {
      // 近傍操作2
      ra1 = 0;
      ra5 = rand_xorshift() % 15 + 1;
      ra2 = rand_xorshift() % 12;
      ra3 = rand_xorshift() % 12;
      int ng = 0;
      while (transition_probs[ra2][ra3] < ra5) {
        ra2 = rand_xorshift() % 12;
        ra3 = rand_xorshift() % 12;
        ng++;
        if (ng > 100)break;
      }
      if (ng > 100)continue;
      ra4 = rand_xorshift() % 12;
      while (ra3 == ra4) {
        ra4 = rand_xorshift() % 12;
      }

      transition_probs[ra2][ra3] -= ra5;
      transition_probs[ra2][ra4] += ra5;
    }
    else if (operation_type < params.operation_thresholds[2]) {
      // 近傍操作2
      ra1 = rand_xorshift() % 6;
      ra2 = rand_xorshift() % 6;
      swap(transition_probs[ra1][ra2], transition_probs[ra1][ra2 + 6]);
    }
    else if (operation_type < params.operation_thresholds[3]) {
      // 近傍操作2
      ra1 = rand_xorshift() % 6;
      for (int j = 0; j < 12; ++j) {
        swap(transition_probs[ra1][j], transition_probs[ra1 + 6][j]);
      }
    }

    // スコア計算
    double tmp_score = evaluate_score(false);

    // 焼きなましで採用判定
    double diff_score = (tmp_score - current_score) * params.score_scale;
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
      if (operation_type < params.operation_thresholds[0]) {
        // 近傍操作1 の巻き戻し
        transition_probs[ra2 + ra1 * 6][ra3 + ra1 * 6] += ra5;
        transition_probs[ra2 + ra1 * 6][ra4 + ra1 * 6] -= ra5;
      }
      else if (operation_type < params.operation_thresholds[1]) {
        // 近傍操作2 の巻き戻し
        transition_probs[ra2][ra3] += ra5;
        transition_probs[ra2][ra4] -= ra5;
      }
      else if (operation_type < params.operation_thresholds[2]) {
        // 近傍操作2 の巻き戻し
        swap(transition_probs[ra1][ra2], transition_probs[ra1][ra2 + 6]);
      }
      else if (operation_type < params.operation_thresholds[3]) {
        // 近傍操作2 の巻き戻し
        for (int j = 0; j < 12; ++j) {
          swap(transition_probs[ra1][j], transition_probs[ra1 + 6][j]);
        }
      }
    }
  }

  cerr << loop << endl;

  restore_best_score();
}

ll solve_case(int case_num, AnnealingParams annealingParams)
{
  start_timer();

  initialize_state();

  input_data(case_num);

  ofstream ofs;
  open_ofs(case_num, ofs);

  build_seed_solution(1.0);

  // 焼きなまし実行
  run_simulated_annealing(annealingParams, 0, 1.5);

  annealingParams.start_temperature[0] = 5000048.0;
  run_simulated_annealing(annealingParams, 1, 1.95);

  // 解答を出力
  output_data(ofs);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (true) {
    score = evaluate_score(false);
  }
  return score;
}

int main()
{
  exec_mode = 2;

  AnnealingParams annealingParams;
  annealingParams.start_temperature[0] = 2000048.0;
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
    cerr << solve_case(0, annealingParams) << endl;
  }
  else if (exec_mode <= 2) {
    ll sum_score = 0;
    for (int i = 0; i < 15; ++i)
    {
      ll score = solve_case(i, annealingParams);
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
