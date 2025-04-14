#include <algorithm>
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
// std::shuffle(v.begin(), v.end(), engine);

const ll INF = 1001001001001001001LL;
const int INT_INF = 1001001001;

const int DX[4] = { -1, 0, 1, 0 };
const int DY[4] = { 0, -1, 0, 1 };


double time_limit = 2.9;
int exec_mode;

static const ll PERFECT_SCORE = 100000000;
static const int N = 20;
static const int MAX_PATTERNS = 1000;
static const int MIN_PATTERN_LENGTH = 2;
static const int MAX_PATTERN_LENGTH = 20;

class Patterns
{
private:
  const int MERGE_REQUIRED_LENGTH = 5;
public:
  vector<vector<vector<int>>> pattern;
  int pattern_count;

  void initialize(vector<string> strs) {
    pattern.clear();
    rep(i, MAX_PATTERN_LENGTH + 1) {
      pattern.push_back(vector<vector<int>>());
    }
    for (int i = 0; i < strs.size(); i++) {
      string s = strs[i];
      int len = s.size();
      vector<int> tmp(len);
      for (int j = 0; j < len; j++) {
        tmp[j] = s[j] - 'A' + 1;
      }
      pattern[len].push_back(tmp);
    }
    for (int i = MIN_PATTERN_LENGTH; i <= MAX_PATTERN_LENGTH; i++) {
      sort(pattern[i].begin(), pattern[i].end());
    }
    pattern_count = 0;
    for (int i = MIN_PATTERN_LENGTH; i <= MAX_PATTERN_LENGTH; i++) {
      pattern_count += pattern[i].size();
    }
  }

  void clear() {
    pattern_count = 0;
    for (int i = MIN_PATTERN_LENGTH; i <= MAX_PATTERN_LENGTH; i++) {
      pattern[i].clear();
    }
  }

  void merge() {
    vector<vector<int>> patterns_tmp;
    vector<vector<int>> patterns_keep;
    for (int i = MIN_PATTERN_LENGTH; i <= MAX_PATTERN_LENGTH; i++) {
      for (auto& vec : pattern[i]) {
        if (vec.size() < MERGE_REQUIRED_LENGTH) {
          patterns_keep.push_back(vec);
        }
        else {
          patterns_tmp.push_back(vec);
        }
      }
    }

    while (true) {
      int idx1 = -1;
      int idx2 = -1;
      int diff = 0;
      int max_len = 0;
      for (int i = 0; i < patterns_tmp.size(); i++) {
        for (int j = 0; j < patterns_tmp.size(); j++) {
          if (i == j) {
            continue;
          }
          int szi = patterns_tmp[i].size();
          int szj = patterns_tmp[j].size();
          rep(k, szi - MERGE_REQUIRED_LENGTH + 1) {
            int new_len = min(szi - k, szj);
            if (new_len < max_len && szi + szj < MAX_PATTERN_LENGTH) {
              break;
            }
            bool ok = true;
            for (int l = 0; l < new_len; l++) {
              if (patterns_tmp[i][l + k] != patterns_tmp[j][l]) {
                ok = false;
                break;
              }
            }

            if (ok && max(szi, szj + k) > MAX_PATTERN_LENGTH) {
              for (int l = 0; l < szj; l++) {
                if (l + k >= MAX_PATTERN_LENGTH) {
                  if (patterns_tmp[i][l + k - MAX_PATTERN_LENGTH] != patterns_tmp[j][l]) {
                    ok = false;
                    break;
                  }
                  else {
                    new_len++;
                  }
                }
              }
            }

            if (ok) {
              if (new_len > max_len) {
                max_len = new_len;
                idx1 = i;
                idx2 = j;
                diff = k;
              }
            }
          }
        }
      }

      if (idx1 == -1) {
        break;
      }

      vector<int> new_vec = patterns_tmp[idx1];
      int szj = patterns_tmp[idx2].size();
      rep(k, szj) {
        if (new_vec.size() >= MAX_PATTERN_LENGTH) {
          break;
        }
        if (k + diff < new_vec.size()) {
          continue;
        }
        else {
          new_vec.push_back(patterns_tmp[idx2][k]);
        }
      }

      if (idx1 < idx2) {
        patterns_tmp.erase(patterns_tmp.begin() + idx2);
        patterns_tmp.erase(patterns_tmp.begin() + idx1);
      }
      else {
        patterns_tmp.erase(patterns_tmp.begin() + idx1);
        patterns_tmp.erase(patterns_tmp.begin() + idx2);
      }
      patterns_tmp.push_back(new_vec);
    }

    clear();
    for (auto& vec : patterns_keep) {
      pattern[vec.size()].push_back(vec);
    }
    for (auto& vec : patterns_tmp) {
      pattern[vec.size()].push_back(vec);
    }
    for (int i = MIN_PATTERN_LENGTH; i <= MAX_PATTERN_LENGTH; i++) {
      sort(pattern[i].begin(), pattern[i].end());
    }

    pattern_count = 0;
    for (int i = MIN_PATTERN_LENGTH; i <= MAX_PATTERN_LENGTH; i++) {
      pattern_count += pattern[i].size();
    }
  }
};

class MatchedFlags
{
private:
  int count;
  bool row_flags[N][N];
  bool col_flags[N][N];

public:
  MatchedFlags() {
    clear();
  }

  void clear() {
    count = 0;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        row_flags[i][j] = false;
        col_flags[i][j] = false;
      }
    }
  }

  void set_flag(int row, int col, int dir, bool b) {
    if (dir == 0) {
      count -= row_flags[row][col];
      row_flags[row][col] = b;
      count += row_flags[row][col];
    }
    else {
      count -= col_flags[row][col];
      col_flags[row][col] = b;
      count += col_flags[row][col];
    }
  }

  int get_count() {
    return count;
  }
};

Patterns patterns;

class State
{
public:
  Patterns& patterns;

  int matched_count;
  int grid[N][N];
  vector<vector<MatchedFlags>> matched_flags;

  State() = delete;
  State(Patterns& p) : patterns(p) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        grid[i][j] = 0;
      }
    }

    matched_flags.resize(MAX_PATTERN_LENGTH + 1);
    for (int i = MIN_PATTERN_LENGTH; i <= MAX_PATTERN_LENGTH; i++) {
      matched_flags[i].resize(patterns.pattern[i].size());
      for (int j = 0; j < patterns.pattern[i].size(); j++) {
        matched_flags[i][j].clear();
      }
    }

    matched_count = 0;
  }

  void generate_random() {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        grid[i][j] = rand_xorshift() % 8 + 1;
      }
    }
  }

  void generate_greedy() {
    for (auto& v : patterns.pattern) {
      for (auto& vec : v) {
        int row = rand_xorshift() % N;
        int col = rand_xorshift() % N;
        for (int i = 0; i < vec.size(); i++) {
          grid[row][(col + i) % N] = vec[i];
        }
      }
    }
  }

  void recalc_all() {
    for (int i = MIN_PATTERN_LENGTH; i <= MAX_PATTERN_LENGTH; i++) {
      for (int j = 0; j < patterns.pattern[i].size(); j++) {
        matched_flags[i][j].clear();
        rep(k, N) {
          rep(l, N) {
            matched_flags[i][j].set_flag(k, l, 0, is_matched(k, l, patterns.pattern[i][j], 0));
            matched_flags[i][j].set_flag(k, l, 1, is_matched(k, l, patterns.pattern[i][j], 1));
          }
        }
      }
    }

    matched_count = 0;
    for (int i = MIN_PATTERN_LENGTH; i <= MAX_PATTERN_LENGTH; i++) {
      for (int j = 0; j < patterns.pattern[i].size(); j++) {
        matched_count += matched_flags[i][j].get_count() > 0;
      }
    }
  }

  void update_one_point(int row, int col) {
    for (int i = MIN_PATTERN_LENGTH; i <= MAX_PATTERN_LENGTH; i++) {
      for (int j = 0; j < patterns.pattern[i].size(); j++) {
        matched_count -= matched_flags[i][j].get_count() > 0;
        rep(k, i) {
          matched_flags[i][j].set_flag(row, (col + N - k) % N, 0, is_matched(row, (col + N - k) % N, patterns.pattern[i][j], 0));
          matched_flags[i][j].set_flag((row + N - k) % N, col, 1, is_matched((row + N - k) % N, col, patterns.pattern[i][j], 1));
        }
        matched_count += matched_flags[i][j].get_count() > 0;
      }
    }
  }

  ll get_score() {
    return PERFECT_SCORE * matched_count / patterns.pattern_count;
  }

private:
  bool is_matched(int row, int col, const vector<int>& vec, int dir) {
    if (dir == 0) {
      for (int i = 0; i < vec.size(); i++) {
        if (grid[row][(col + i) % N] != vec[i]) {
          return false;
        }
      }
    }
    else {
      for (int i = 0; i < vec.size(); i++) {
        if (grid[(row + i) % N][col] != vec[i]) {
          return false;
        }
      }
    }
    return true;
  }
};


void store_best_score()
{
}

void restore_best_score()
{
}

bool is_out_of_range(int x, int y)
{
  //if (x < 0 || n <= x || y < 0 || n <= y) return true;
  return false;
}

void initialize_state()
{
}

void input_data(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());


  int _n, _m;
  vector<string> vs;
  if (!ifs.is_open()) {
    // 標準入力
    cin >> _n >> _m;
    rep(i, _m)
    {
      string s;
      cin >> s;
      vs.push_back(s);
    }
  }
  else {
    // ファイル入力
    ifs >> _n >> _m;
    rep(i, _m)
    {
      string s;
      ifs >> s;
      vs.push_back(s);
    }
  }
  patterns.initialize(vs);
}

void open_ofs(int case_num, ofstream& ofs)
{
  if (exec_mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofs.open(oss.str());
  }
}

void output_data(ofstream& ofs, State& state)
{
  if (exec_mode == 0) {
    // 標準出力
    rep(i, N)
    {
      rep(j, N)
      {
        if (state.grid[i][j] == 0) {
          cout << '.';
        }
        else {
          cout << (char)(state.grid[i][j] + 'A' - 1);
        }
      }
      cout << endl;
    }
  }
  else {
    // ファイル出力
    rep(i, N)
    {
      rep(j, N)
      {
        if (state.grid[i][j] == 0) {
          ofs << '.';
        }
        else {
          ofs << (char)(state.grid[i][j] + 'A' - 1);
        }
      }
      ofs << endl;
    }
  }
}

void build_initial_solution()
{
}

struct AnnealingParams
{
  double start_temperature[10];
  double end_temperature;
  double score_scale;
  int operation_threshold[10];
};

void run_simulated_annealing(AnnealingParams annealingParams, State& state)
{
  store_best_score();

  double now_time = get_elapsed_time();
  const double START_TEMP = annealingParams.start_temperature[1];
  const double END_TEMP = annealingParams.end_temperature;

  int iteration_count = 0;
  while (true) {
    iteration_count++;

    if (iteration_count % 100 == 0) {
      now_time = get_elapsed_time();
      if (now_time > time_limit) break;
    }

    double progress_ratio = now_time / time_limit;
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progress_ratio;

    // 近傍解作成
    int ra_exec_mode = rand_xorshift() % annealingParams.operation_threshold[1];
    int ra1, ra2, ra3, ra4, ra5;
    int keep1, keep2, keep3, keep4, keep5;

    ll current_score = state.get_score();

    int row = rand_xorshift() % N;
    int col = rand_xorshift() % N;

    int candidate_value = rand_xorshift() % 8 + 1;
    int old_value = state.grid[row][col];

    int len = rand_xorshift() % MAX_PATTERN_LENGTH + 1;
    while (patterns.pattern[len].size() == 0) {
      len = rand_xorshift() % MAX_PATTERN_LENGTH + 1;
    }
    int index = rand_xorshift() % patterns.pattern[len].size();
    int dir = rand_xorshift() % 2;
    //dir = 0;
    vector<int> old_values(len);
    for (int i = 0; i < len; i++) {
      if (dir == 0) {
        old_values[i] = state.grid[row][(col + i) % N];
      }
      else {
        old_values[i] = state.grid[(row + i) % N][col];
      }
    }

    if (ra_exec_mode < annealingParams.operation_threshold[0]) {
      // 近傍操作1
      state.grid[row][col] = candidate_value;
      state.update_one_point(row, col);
    }
    else if (ra_exec_mode < annealingParams.operation_threshold[1]) {
      // 近傍操作2
      for (int i = 0; i < len; i++) {
        if (dir == 0) {
          state.grid[row][(col + i) % N] = patterns.pattern[len][index][i];
          state.update_one_point(row, (col + i) % N);
        }
        else {
          state.grid[(row + i) % N][col] = patterns.pattern[len][index][i];;
          state.update_one_point((row + i) % N, col);
        }
      }
    }

    // スコア計算
    double tmp_score = state.get_score();

    // 焼きなましで採用判定
    double diff_score = (tmp_score - current_score) * annealingParams.score_scale;
    double prob = exp(diff_score / temp);
    if (prob > rand_01()) {
      // 採用
    }
    else {
      // 元に戻す
      if (ra_exec_mode < annealingParams.operation_threshold[0]) {
        // 近傍操作1 の巻き戻し
        state.grid[row][col] = old_value;
        state.update_one_point(row, col);
      }
      else if (ra_exec_mode < annealingParams.operation_threshold[1]) {
        // 近傍操作2 の巻き戻し
        for (int i = 0; i < len; i++) {
          if (dir == 0) {
            state.grid[row][(col + i) % N] = old_values[i];
            state.update_one_point(row, (col + i) % N);
          }
          else {
            state.grid[(row + i) % N][col] = old_values[i];
            state.update_one_point((row + i) % N, col);
          }
        }
      }
    }
  }

  if (exec_mode != 0 && exec_mode != 3) {
    cout << iteration_count << endl;
  }

  restore_best_score();
}

ll solve_case(int case_num, AnnealingParams annealingParams)
{
  start_timer();

  initialize_state();

  input_data(case_num);

  patterns.merge();

  State state(patterns);

  ofstream ofs;
  open_ofs(case_num, ofs);

  state.generate_random();
  state.generate_greedy();
  state.recalc_all();

  // 焼きなまし実行
  run_simulated_annealing(annealingParams, state);

  // 解答を出力
  output_data(ofs, state);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (exec_mode != 0) {
    score = state.get_score();
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
  annealingParams.score_scale = 1234567.0;
  annealingParams.operation_threshold[0] = 190;
  annealingParams.operation_threshold[1] = 200;
  annealingParams.operation_threshold[2] = 300;
  annealingParams.operation_threshold[3] = 400;
  annealingParams.operation_threshold[4] = 500;
  annealingParams.operation_threshold[5] = 600;
  annealingParams.operation_threshold[6] = 700;
  annealingParams.operation_threshold[7] = 800;
  annealingParams.operation_threshold[8] = 900;
  annealingParams.operation_threshold[9] = 1000;

  if (exec_mode == 0) {
    solve_case(0, annealingParams);
  }
  else if (exec_mode <= 2) {
    ll sum_score = 0;
    srep(i, 0, 15)
    {
      ll score = solve_case(i, annealingParams);
      sum_score += score;
      if (exec_mode == 1) {
        cout << score << endl;
      }
      else {
        cout << "case = " << setw(2) << i << ", "
          << "score = " << setw(4) << score << ", "
          << "sum = " << setw(5) << sum_score << ", "
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
      new_annealingParams.operation_threshold[0] = rand() % 101;

      ll sum_score = 0;
      srep(i, 0, 15)
      {
        ll score = solve_case(i, new_annealingParams);
        sum_score += score;

        // シード0が悪ければ打ち切り
        if (i == 0 && score < 0) {
          break;
        }
      }

      cout << "loop_count = " << loop_count
        << ", sum_score = " << sum_score
        << ", start_temperature = " << new_annealingParams.start_temperature[0]
        << ", end_temperature = " << new_annealingParams.end_temperature
        << ", score_scale = " << new_annealingParams.score_scale
        << ", operation_threshold = " << new_annealingParams.operation_threshold[0]
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
