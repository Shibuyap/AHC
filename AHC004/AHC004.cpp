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

double TIME_LIMIT = 2.9;
int exec_mode;

static const ll PERFECT_SCORE = 100000000;
static const int N = 20;
static const int MAX_PATTERNS = 1000;
static const int MIN_PATTERN_LENGTH = 2;
static const int MAX_PATTERN_LENGTH = 20;
static const int CHARACTER_SIZE = 8;

struct Pattern
{
public:
  vector<int> pattern;
  int merged_count;

  Pattern() = delete;
  Pattern(vector<int> p) : pattern(p), merged_count(1) {}
};

struct Patterns
{
public:
  vector<vector<Pattern>> vv_patterns;
  int pattern_count = 0;
};

class PatternsManager
{
public:
  Patterns initial_patterns;
  Patterns merged_patterns;

private:
  const int MERGE_UNUSE_LENGTH = 1;

public:
  void initialize(vector<string> strs) {
    initial_patterns.vv_patterns.clear();
    initial_patterns.vv_patterns.resize(MAX_PATTERN_LENGTH + 1);
    for (int i = 0; i < strs.size(); i++) {
      string s = strs[i];
      int len = s.size();
      vector<int> tmp(len);
      for (int j = 0; j < len; j++) {
        tmp[j] = s[j] - 'A';
      }
      initial_patterns.vv_patterns[len].emplace_back(tmp);
    }
    initial_patterns.pattern_count = 0;
    for (auto& v_patterns : initial_patterns.vv_patterns) {
      for (auto& pattern : v_patterns) {
        initial_patterns.pattern_count += pattern.merged_count;
      }
    }
  }

  void build_merge_patterns(double time_limit) {
    vector<Pattern> patterns_tmp;
    vector<Pattern> patterns_keep;
    for (auto& v_patterns : initial_patterns.vv_patterns) {
      for (auto& pattern : v_patterns) {
        if (pattern.pattern.size() < MERGE_UNUSE_LENGTH) {
          patterns_keep.emplace_back(pattern);
        }
        else {
          patterns_tmp.emplace_back(pattern);
        }
      }
    }

    while (true) {
      if (get_elapsed_time() > time_limit) {
        break;
      }

      int idx1 = -1;
      int idx2 = -1;
      int diff = 0;
      int max_len = 0;
      for (int i = patterns_tmp.size() - 1; i >= 0; i--) {
        for (int j = 0; j < patterns_tmp.size(); j++) {
          if (i == j) {
            continue;
          }
          int size_i = patterns_tmp[i].pattern.size();
          int size_j = patterns_tmp[j].pattern.size();
          rep(k, size_i) {
            int new_len = min(size_i - k, size_j);
            if (new_len < max_len && size_i + size_j < MAX_PATTERN_LENGTH) {
              break;
            }
            bool ok = true;
            for (int l = 0; l < new_len; l++) {
              if (patterns_tmp[i].pattern[l + k] != patterns_tmp[j].pattern[l]) {
                ok = false;
                break;
              }
            }

            if (ok && max(size_i, size_j + k) > MAX_PATTERN_LENGTH) {
              for (int l = 0; l < size_j; l++) {
                if (l + k >= MAX_PATTERN_LENGTH) {
                  if (patterns_tmp[i].pattern[l + k - MAX_PATTERN_LENGTH] != patterns_tmp[j].pattern[l]) {
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
        if (max_len >= 10) {
          break;
        }
      }

      if (idx1 == -1) {
        break;
      }

      Pattern merged_pattern = patterns_tmp[idx1];
      int size_j = patterns_tmp[idx2].pattern.size();
      rep(k, size_j) {
        if (merged_pattern.pattern.size() >= MAX_PATTERN_LENGTH) {
          break;
        }
        if (k + diff < merged_pattern.pattern.size()) {
          continue;
        }
        else {
          merged_pattern.pattern.push_back(patterns_tmp[idx2].pattern[k]);
        }
      }
      merged_pattern.merged_count += patterns_tmp[idx2].merged_count;

      if (idx1 < idx2) {
        patterns_tmp.erase(patterns_tmp.begin() + idx2);
        patterns_tmp.erase(patterns_tmp.begin() + idx1);
      }
      else {
        patterns_tmp.erase(patterns_tmp.begin() + idx1);
        patterns_tmp.erase(patterns_tmp.begin() + idx2);
      }
      patterns_tmp.push_back(merged_pattern);
    }

    merged_patterns.vv_patterns.clear();
    merged_patterns.vv_patterns.resize(MAX_PATTERN_LENGTH + 1);
    for (auto& pattern : patterns_keep) {
      merged_patterns.vv_patterns[pattern.pattern.size()].push_back(pattern);
    }
    for (auto& pattern : patterns_tmp) {
      merged_patterns.vv_patterns[pattern.pattern.size()].push_back(pattern);
    }

    merged_patterns.pattern_count = 0;
    for (auto& v_patterns : merged_patterns.vv_patterns) {
      for (auto& pattern : v_patterns) {
        merged_patterns.pattern_count += pattern.merged_count;
      }
    }
  }
};

class MatchedFlag
{
private:
  int count;
  bool row_flags[N][N];
  bool col_flags[N][N];

public:
  MatchedFlag() {
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

class State
{
public:
  PatternsManager& patterns;

  int matched_count;
  int grid[N][N];
  vector<vector<MatchedFlag>> matched_flags;

  State() = delete;
  State(PatternsManager& p) : patterns(p) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        grid[i][j] = 0;
      }
    }

    reset_matched_flags(patterns.initial_patterns);
  }

  void generate_random() {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        grid[i][j] = rand_xorshift() % CHARACTER_SIZE;
      }
    }
  }

  void generate_random_empty() {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        if (grid[i][j] == CHARACTER_SIZE) {
          grid[i][j] = rand_xorshift() % CHARACTER_SIZE;
        }
      }
    }
  }

  void greedy_after_assemble(const Patterns& patterns) {
    recalc_all(patterns);

    // 置けてないパターンを長い順に置けるところに置いていく
    for (int len = MAX_PATTERN_LENGTH; len >= MIN_PATTERN_LENGTH; len--) {
      for (int pat_index = 0; pat_index < patterns.vv_patterns[len].size(); pat_index++) {
        if (matched_flags[len][pat_index].get_count() > 0) {
          continue;
        }

        int ok = 0;
        rep(i, N) {
          rep(j, N) {
            // 行に置けるか
            {
              ok = 1;
              rep(k, len) {
                if (patterns.vv_patterns[len][pat_index].pattern[k] != grid[i][(j + k) % N] && grid[i][(j + k) % N] != CHARACTER_SIZE) {
                  ok = 0;
                  break;
                }
              }
              if (ok) {
                rep(k, len) {
                  grid[i][(j + k) % N] = patterns.vv_patterns[len][pat_index].pattern[k];
                }
                break;
              }
            }
            // 列に置けるか
            {
              ok = 1;
              rep(k, len) {
                if (patterns.vv_patterns[len][pat_index].pattern[k] != grid[(i + k) % N][j] && grid[(i + k) % N][j] != CHARACTER_SIZE) {
                  ok = 0;
                  break;
                }
              }
              if (ok) {
                rep(k, len) {
                  grid[(i + k) % N][j] = patterns.vv_patterns[len][pat_index].pattern[k];
                }
                break;
              }
            }
          }
          if (ok) {
            break;
          }
        }
      }
    }
  }

  void reset_matched_flags(const Patterns& patterns) {
    matched_flags.resize(MAX_PATTERN_LENGTH + 1);
    for (int i = MIN_PATTERN_LENGTH; i <= MAX_PATTERN_LENGTH; i++) {
      matched_flags[i].resize(patterns.vv_patterns[i].size());
      for (int j = 0; j < patterns.vv_patterns[i].size(); j++) {
        matched_flags[i][j].clear();
      }
    }

    matched_count = 0;
  }

  void recalc_all(const Patterns& patterns) {
    reset_matched_flags(patterns);

    for (int i = MIN_PATTERN_LENGTH; i <= MAX_PATTERN_LENGTH; i++) {
      for (int j = 0; j < patterns.vv_patterns[i].size(); j++) {
        matched_flags[i][j].clear();
        rep(k, N) {
          rep(l, N) {
            matched_flags[i][j].set_flag(k, l, 0, is_matched(k, l, patterns.vv_patterns[i][j].pattern, 0));
            matched_flags[i][j].set_flag(k, l, 1, is_matched(k, l, patterns.vv_patterns[i][j].pattern, 1));
          }
        }
      }
    }

    matched_count = 0;
    for (int i = MIN_PATTERN_LENGTH; i <= MAX_PATTERN_LENGTH; i++) {
      for (int j = 0; j < patterns.vv_patterns[i].size(); j++) {
        if (matched_flags[i][j].get_count() > 0) {
          matched_count += patterns.vv_patterns[i][j].merged_count;
        }
      }
    }
  }

  void update_one_point(int row, int col, const Patterns& patterns) {
    for (int i = MIN_PATTERN_LENGTH; i <= MAX_PATTERN_LENGTH; i++) {
      for (int j = 0; j < patterns.vv_patterns[i].size(); j++) {
        if (matched_flags[i][j].get_count() > 0) {
          matched_count -= patterns.vv_patterns[i][j].merged_count;
        }
        rep(k, i) {
          matched_flags[i][j].set_flag(row, (col + N - k) % N, 0, is_matched(row, (col + N - k) % N, patterns.vv_patterns[i][j].pattern, 0));
          matched_flags[i][j].set_flag((row + N - k) % N, col, 1, is_matched((row + N - k) % N, col, patterns.vv_patterns[i][j].pattern, 1));
        }
        if (matched_flags[i][j].get_count() > 0) {
          matched_count += patterns.vv_patterns[i][j].merged_count;
        }
      }
    }
  }

  ll get_score(const Patterns& patterns) {
    return PERFECT_SCORE * matched_count / patterns.pattern_count;
  }

  void assemble(double time_limit, const Patterns& patterns) {
    vector<vector<int>> vv;
    drep(i, MAX_PATTERN_LENGTH + 1) {
      for (auto& pattern : patterns.vv_patterns[i]) {
        vv.push_back(pattern.pattern);
        if (vv.size() == N * 2) {
          break;
        }
      }
      if (vv.size() == N * 2) {
        break;
      }
    }

    vector<vector<vector<P>>> prefixMaps;
    rep(i, 6) {
      prefixMaps.push_back(vector<vector<P>>(1 << (i * 3)));
      if (i == 0) {
        continue;
      }
      srep(j, 0, vv.size()) {
        rep(k, vv[j].size()) {
          int num = 0;
          rep(l, i) {
            num *= CHARACTER_SIZE;
            num += vv[j][(k + l) % vv[j].size()];
          }
          prefixMaps[i][num].emplace_back(j, k);
        }
      }
    }

    int best_g[N][N];
    int best_score = -1;

    int g[N][N];
    rep(i, N) {
      rep(j, N) {
        g[i][j] = CHARACTER_SIZE;
      }
    }

    int used[N * 2] = {};
    int used_version = 0;
    int decided_col[N] = {};
    int decided_col_2[N] = {};
    int decided_col_count = 0;
    int decided_row[N] = {};
    int decided_row_2[N] = {};
    int decided_row_count = 0;

    // 1行目は一番長いので固定

    srep(ii, 0, vv.size()) {
      rep(j, vv[ii].size()) {
        g[0][j] = vv[ii][j];
      }
      srep(i, 0, vv.size()) {
        srep(j, 0, vv.size()) {
          if (i == j || ii == i || ii == j) {
            continue;
          }
          rep(k, N) {
            rep(l, N) {
              // 2行目と3行目をセット
              rep(m, vv[i].size()) {
                g[1][(k + m) % N] = vv[i][m];
              }
              rep(m, vv[j].size()) {
                g[2][(l + m) % N] = vv[j][m];
              }

              // 各列を決定していく
              rep(m, N) {
                decided_col[m] = -1;
              }
              decided_col_count = 0;
              used_version++;
              used[ii] = used_version;
              used[i] = used_version;
              used[j] = used_version;
              int ng_count = 0;
              while (true) {
                bool has_change = false;
                rep(m, N) {
                  if (decided_col[m] != -1) {
                    continue;
                  }
                  if (g[0][m] == CHARACTER_SIZE || g[1][m] == CHARACTER_SIZE || g[2][m] == CHARACTER_SIZE) {
                    continue;
                  }
                  int num = g[0][m] * CHARACTER_SIZE * CHARACTER_SIZE + g[1][m] * CHARACTER_SIZE + g[2][m];

                  int idx1 = -1;
                  int idx2 = -1;
                  for (auto& pre : prefixMaps[3][num]) {
                    if (used[pre.first] == used_version) {
                      continue;
                    }
                    else {
                      if (idx1 == -1) {
                        idx1 = pre.first;
                        idx2 = pre.second;
                      }
                      else {
                        idx1 = -2;
                        break;
                      }
                    }
                  }

                  if (idx1 == -1) {
                    ng_count++;
                    if (ng_count >= 1) {
                      break;
                    }
                  }
                  else if (idx1 >= 0) {
                    decided_col[m] = idx1;
                    decided_col_2[m] = idx2;
                    decided_col_count++;
                    used[idx1] = used_version;
                    has_change = true;
                  }
                }

                if (ng_count >= 1) {
                  break;
                }
                if (!has_change) {
                  break;
                }
              }

              // 失敗したら元に戻す
              if (ng_count >= 1) {
                rep(m, vv[i].size()) {
                  g[1][(k + m) % N] = CHARACTER_SIZE;
                }
                rep(m, vv[j].size()) {
                  g[2][(l + m) % N] = CHARACTER_SIZE;
                }
                continue;
              }

              //cout << i << " " << j << " " << k << " " << l << " decided_col_count: " << decided_col_count << endl;

              // 行決定に使う列を決定する
              const int COL_LENGTH = 1;
              int start_col = -1;
              rep(m, N) {
                bool ok = true;
                rep(n, COL_LENGTH) {
                  if (decided_col[(m + n) % N] == -1) {
                    ok = false;
                    break;
                  }
                }
                if (ok) {
                  start_col = m;
                  break;
                }
              }
              if (start_col != -1) {
                rep(mm, COL_LENGTH) {
                  int m = (start_col + mm) % N;
                  rep(n, N) {
                    int idx = (decided_col_2[m] + n) % N;
                    g[n][m] = idx < vv[decided_col[m]].size() ? vv[decided_col[m]][idx] : CHARACTER_SIZE;
                  }
                }
              }

              // 各行を決定していく
              rep(m, N) {
                decided_row[m] = -1;
              }
              decided_row_count = 0;
              ng_count = 0;
              while (start_col != -1) {
                bool has_change = false;
                srep(m, 3, N) {
                  if (decided_row[m] != -1) {
                    continue;
                  }

                  bool is_blank = false;
                  int num = 0;
                  rep(n, COL_LENGTH) {
                    if (g[m][(start_col + n) % N] == CHARACTER_SIZE) {
                      is_blank = true;
                      break;
                    }
                    num *= CHARACTER_SIZE;
                    num += g[m][(start_col + n) % N];
                  }
                  if (is_blank) {
                    continue;
                  }

                  int idx1 = -1;
                  int idx2 = -1;
                  for (auto& pre : prefixMaps[COL_LENGTH][num]) {
                    if (used[pre.first] == used_version) {
                      continue;
                    }
                    else {
                      bool is_ok = true;

                      rep(n, N) {
                        if (decided_col[n] != -1) {
                          int idx_col = (decided_col_2[n] + m) % N;
                          int idx_row = (pre.second + n) % N;
                          if (idx_col < vv[decided_col[n]].size() && idx_row < vv[pre.first].size()) {
                            if (vv[pre.first][idx_row] != vv[decided_col[n]][idx_col]) {
                              is_ok = false;
                              break;
                            }
                          }
                        }
                      }

                      if (is_ok) {
                        if (idx1 == -1) {
                          idx1 = pre.first;
                          idx2 = pre.second;
                        }
                        else {
                          idx1 = -2;
                          break;
                        }
                      }
                      else {
                        continue;
                      }
                    }
                  }

                  if (idx1 == -1) {
                    ng_count++;
                    //break;
                  }
                  else if (idx1 >= 0) {
                    decided_row[m] = idx1;
                    decided_row_2[m] = idx2;
                    decided_row_count++;
                    used[idx1] = used_version;
                    has_change = true;
                  }
                }

                if (!has_change) {
                  break;
                }
              }

              if (ng_count <= 10000 && decided_col_count + decided_row_count * 100 > best_score) {
                best_score = decided_col_count + decided_row_count * 100;
                rep(n, N) {
                  if (n < 3) {
                    rep(m, N) {
                      best_g[n][m] = g[n][m];
                    }
                  }
                  else {
                    rep(m, N) {
                      if (decided_col[m] != -1) {
                        int idx = (decided_col_2[m] + n) % N;
                        best_g[n][m] = idx < vv[decided_col[m]].size() ? vv[decided_col[m]][idx] : CHARACTER_SIZE;
                      }
                      else {
                        best_g[n][m] = CHARACTER_SIZE;
                      }
                    }
                  }
                }

                srep(n, 3, N) {
                  if (decided_row[n] != -1) {
                    rep(m, N) {
                      int idx = (decided_row_2[n] + m) % N;
                      if (idx < vv[decided_row[n]].size()) {
                        best_g[n][m] = vv[decided_row[n]][idx];
                      }
                    }
                  }
                }
              }

              // 後片付け
              rep(m, vv[i].size()) {
                g[1][(k + m) % N] = CHARACTER_SIZE;
              }
              rep(m, vv[j].size()) {
                g[2][(l + m) % N] = CHARACTER_SIZE;
              }
            }
          }
        }
      }

      if (best_score >= (N - 3) * 100) {
        break;
      }
      if (get_elapsed_time() > time_limit) {
        break;
      }
    }
    if (best_score > 0) {
      rep(i, N) {
        rep(j, N) {
          grid[i][j] = best_g[i][j];
        }
      }
    }
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

void input_data(int case_num, PatternsManager& patterns_manager)
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
  patterns_manager.initialize(vs);
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
        if (state.grid[i][j] == CHARACTER_SIZE) {
          cout << '.';
        }
        else {
          cout << (char)(state.grid[i][j] + 'A');
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
        if (state.grid[i][j] == CHARACTER_SIZE) {
          ofs << '.';
        }
        else {
          ofs << (char)(state.grid[i][j] + 'A');
        }
      }
      ofs << endl;
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

void run_simulated_annealing(AnnealingParams annealingParams, State& state, const Patterns& patterns)
{
  double now_time = get_elapsed_time();
  const double START_TEMP = annealingParams.start_temperature[0];
  const double END_TEMP = annealingParams.end_temperature;

  int iteration_count = 0;
  while (true) {
    iteration_count++;

    if (iteration_count % 10 == 0) {
      now_time = get_elapsed_time();
      if (now_time > TIME_LIMIT) break;
    }

    double progress_ratio = now_time / TIME_LIMIT;
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progress_ratio;

    // 近傍解作成
    int ra_exec_mode = rand_xorshift() % annealingParams.operation_thresholds[1];
    int ra1, ra2, ra3, ra4, ra5;
    int keep1, keep2, keep3, keep4, keep5;

    ll current_score = state.get_score(patterns);

    int row = rand_xorshift() % N;
    int col = rand_xorshift() % N;

    int candidate_value = rand_xorshift() % CHARACTER_SIZE;
    int old_value = state.grid[row][col];

    int len = rand_xorshift() % MAX_PATTERN_LENGTH + 1;
    while (patterns.vv_patterns[len].size() == 0) {
      len = rand_xorshift() % MAX_PATTERN_LENGTH + 1;
    }
    int index = rand_xorshift() % patterns.vv_patterns[len].size();
    int dir = rand_xorshift() % 2;

    int row2 = rand_xorshift() % N;
    while (row2 == row) {
      row2 = rand_xorshift() % N;
    }

    vector<int> old_values(len);
    for (int i = 0; i < len; i++) {
      if (dir == 0) {
        old_values[i] = state.grid[row][(col + i) % N];
      }
      else {
        old_values[i] = state.grid[(row + i) % N][col];
      }
    }

    if (ra_exec_mode < annealingParams.operation_thresholds[0]) {
      // 近傍操作1
      state.grid[row][col] = candidate_value;
      state.update_one_point(row, col, patterns);
    }
    else if (ra_exec_mode < annealingParams.operation_thresholds[1]) {
      // 近傍操作2
      for (int i = 0; i < len; i++) {
        if (dir == 0) {
          state.grid[row][(col + i) % N] = patterns.vv_patterns[len][index].pattern[i];
          state.update_one_point(row, (col + i) % N, patterns);
        }
        else {
          state.grid[(row + i) % N][col] = patterns.vv_patterns[len][index].pattern[i];;
          state.update_one_point((row + i) % N, col, patterns);
        }
      }
    }
    else if (ra_exec_mode < annealingParams.operation_thresholds[2]) {
      // 近傍操作3
      int keep[N];
      for (int i = 0; i < N; i++) {
        keep[i] = state.grid[row][i];
      }
      for (int i = 0; i < N; i++) {
        state.grid[row][i] = state.grid[row2][i];
        state.update_one_point(row, i, patterns);
      }
      for (int i = 0; i < N; i++) {
        state.grid[row2][i] = keep[i];
        state.update_one_point(row2, i, patterns);
      }
    }

    // スコア計算
    double tmp_score = state.get_score(patterns);

    // 焼きなましで採用判定
    double diff_score = (tmp_score - current_score) * annealingParams.score_scale;
    double prob = exp(diff_score / temp);
    if (prob > rand_01()) {
      // 採用
    }
    else {
      // 元に戻す
      if (ra_exec_mode < annealingParams.operation_thresholds[0]) {
        // 近傍操作1 の巻き戻し
        state.grid[row][col] = old_value;
        state.update_one_point(row, col, patterns);
      }
      else if (ra_exec_mode < annealingParams.operation_thresholds[1]) {
        // 近傍操作2 の巻き戻し
        for (int i = 0; i < len; i++) {
          if (dir == 0) {
            state.grid[row][(col + i) % N] = old_values[i];
            state.update_one_point(row, (col + i) % N, patterns);
          }
          else {
            state.grid[(row + i) % N][col] = old_values[i];
            state.update_one_point((row + i) % N, col, patterns);
          }
        }
      }
      else if (ra_exec_mode < annealingParams.operation_thresholds[2]) {
        // 近傍操作3 の巻き戻し
        int keep[N];
        for (int i = 0; i < N; i++) {
          keep[i] = state.grid[row][i];
        }
        for (int i = 0; i < N; i++) {
          state.grid[row][i] = state.grid[row2][i];
          state.update_one_point(row, i, patterns);
        }
        for (int i = 0; i < N; i++) {
          state.grid[row2][i] = keep[i];
          state.update_one_point(row2, i, patterns);
        }
      }
    }
  }

  if (exec_mode >= 3) {
    cerr << "iteration_count = " << iteration_count << endl;
  }
  if (exec_mode >= 2) {
    cerr << "焼きなまし処理終了 : " << get_elapsed_time() << " sec" << endl;
  }
}

ll solve_case(int case_num, AnnealingParams annealingParams)
{
  start_timer();

  PatternsManager patterns_manager;

  input_data(case_num, patterns_manager);

  patterns_manager.build_merge_patterns(TIME_LIMIT * 0.6);

  if (exec_mode >= 3) {
    for (int i = 0; i < patterns_manager.merged_patterns.vv_patterns.size(); i++) {
      cerr << setw(3) << i << " ";
    }
    cerr << endl;
    for (int i = 0; i < patterns_manager.merged_patterns.vv_patterns.size(); i++) {
      cerr << setw(3) << patterns_manager.merged_patterns.vv_patterns[i].size() << " ";
    }
    cerr << endl;
  }
  if (exec_mode >= 2) {
    cerr << "パターンマージ処理終了 : " << get_elapsed_time() << " sec" << endl;
  }

  ofstream ofs;
  open_ofs(case_num, ofs);

  State state(patterns_manager);

  state.assemble(TIME_LIMIT * 0.7, patterns_manager.merged_patterns);

  if (exec_mode >= 3) {
    rep(i, N)
    {
      rep(j, N)
      {
        if (state.grid[i][j] == CHARACTER_SIZE) {
          cerr << '.';
        }
        else {
          cerr << (char)(state.grid[i][j] + 'A');
        }
      }
      cerr << endl;
    }
  }
  if (exec_mode >= 2) {
    cerr << "完全復元処理終了 : " << get_elapsed_time() << " sec" << endl;
  }

  state.recalc_all(patterns_manager.merged_patterns);
  if (state.get_score(patterns_manager.merged_patterns) < PERFECT_SCORE) {
    //state.recalc_all(patterns_manager.initial_patterns);
    //cerr << state.get_score(patterns_manager.initial_patterns) << endl;

    state.greedy_after_assemble(patterns_manager.merged_patterns);
    //state.recalc_all(patterns_manager.initial_patterns);
    //cerr << state.get_score(patterns_manager.initial_patterns) << endl;

    state.greedy_after_assemble(patterns_manager.initial_patterns);
    //state.recalc_all(patterns_manager.initial_patterns);
    //cerr << state.get_score(patterns_manager.initial_patterns) << endl;

    state.generate_random_empty();
    state.recalc_all(patterns_manager.initial_patterns);
    if (state.get_score(patterns_manager.initial_patterns) < PERFECT_SCORE) {
      // 焼きなまし実行
      run_simulated_annealing(annealingParams, state, patterns_manager.initial_patterns);
    }
  }

  // 解答を出力
  output_data(ofs, state);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (exec_mode != 0) {
    state.recalc_all(patterns_manager.initial_patterns);
    score = state.get_score(patterns_manager.initial_patterns);
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
  else if (exec_mode < 100) {
    ll sum_score = 0;
    srep(i, 0, 100)
    {
      if (exec_mode == 1) {
      }
      else {
        cerr << endl;
      }
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
  else if (exec_mode == 100) {
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
      srep(i, 0, 15)
      {
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
