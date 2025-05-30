#include <chrono>
#include <climits>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <utility>
#include <vector>

using namespace std;

typedef long long int ll;
typedef pair<int, int> P;
typedef pair<P, P> PP;

// タイマー
namespace {
  std::chrono::steady_clock::time_point start_time_clock;

  void start_timer() {
    start_time_clock = std::chrono::steady_clock::now();
  }

  double get_elapsed_time() {
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time_clock;
    return elapsed.count();
  }
}

// 乱数
namespace {
  static uint32_t rand_xorshift() {
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

  static double rand_01() {
    return (rand_xorshift() + 0.5) * (1.0 / UINT_MAX);
  }

  static double rand_range(double l, double r) {
    return l + (r - l) * rand_01();
  }

  static uint32_t rand_range(uint32_t l, uint32_t r) {
    return l + rand_xorshift() % (r - l + 1); // [l, r]
  }

  void shuffle_array(int* arr, int n) {
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
int L;

char true_genome[N][N];
vector<string> true_genome_strs;
void generate_true_genome(int l, int m) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      true_genome[i][j] = 'A' + rand_xorshift() % CHARACTER_SIZE;
    }
  }
  true_genome_strs.clear();
  for (int _ = 0; _ < m; ++_) {
    string s;

    int i = rand_xorshift() % N;
    int j = rand_xorshift() % N;
    int dir = rand_xorshift() % 2;
    int len = rand_range((uint32_t)l - 2, (uint32_t)l + 2);
    for (int k = 0; k < len; ++k) {
      if (dir == 0) {
        s += true_genome[i][(j + k) % N];
      }
      else {
        s += true_genome[(i + k) % N][j];
      }
    }
    true_genome_strs.push_back(s);
  }
}

struct Pattern {
public:
  vector<int> pattern;
  int merged_count;

  Pattern() = delete;
  Pattern(vector<int> p) : pattern(p), merged_count(1) {}
};

struct Patterns {
public:
  vector<vector<Pattern>> vv_patterns;
  int pattern_count = 0;
};

class PatternsManager {
public:
  Patterns initial_patterns;
  Patterns merged_patterns;

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

  void build_merge_patterns(double time_limit, int need_length) {
    vector<Pattern> patterns_tmp;
    for (auto& v_patterns : initial_patterns.vv_patterns) {
      for (auto& pattern : v_patterns) {
        patterns_tmp.emplace_back(pattern);
      }
    }

    {
      vector<Pattern> patterns_tmp_2;
      for (int i = 0; i < patterns_tmp.size(); ++i) {
        bool skip = false;
        for (int j = 0; j < patterns_tmp.size(); ++j) {
          if (i == j) {
            continue;
          }
          const auto& p1 = patterns_tmp[i].pattern;
          const auto& p2 = patterns_tmp[j].pattern;
          if (p1.size() > p2.size()) {
            continue;
          }
          for (int k = 0; k < p2.size() - p1.size() + 1; ++k) {
            bool ok = true;
            for (int l = 0; l < p1.size(); l++) {
              if (p1[l] != p2[l + k]) {
                ok = false;
                break;
              }
            }
            if (ok) {
              skip = true;
              break;
            }
          }
          if (skip) {
            break;
          }
        }
        if (!skip) {
          patterns_tmp_2.push_back(patterns_tmp[i]);
        }
      }

      patterns_tmp = patterns_tmp_2;
    }

    while (true) {
      if (get_elapsed_time() > time_limit) {
        break;
      }

      int idx1 = -1;
      int idx2 = -1;
      int diff = 0;
      int max_len = 0;
      int max_min_len = 9999;
      for (int i = patterns_tmp.size() - 1; i >= 0; --i) {
        int size_i = patterns_tmp[i].pattern.size();
        if (size_i < max_len) {
          continue;
        }
        for (int j = patterns_tmp.size() - 1; j >= 0; --j) {
          if (i == j) {
            continue;
          }
          int size_j = patterns_tmp[j].pattern.size();
          for (int k = 0; k < size_i; ++k) {
            int new_len = min(size_i - k, size_j);
            int new_len2 = max(0, size_j + k - N);
            if (new_len + new_len2 <= max_len) {
              break;
            }

            if (new_len <= need_length && new_len < size_j) {
              continue;
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
                }
              }
            }

            if (ok) {
              if (new_len + new_len2 > max_len) {
                max_len = new_len + new_len2;
                max_min_len = max(size_i, size_j);
                idx1 = i;
                idx2 = j;
                diff = k;
              }
              else if (new_len + new_len2 == max_len && max(size_j, size_j) < max_min_len) {
                max_len = new_len + new_len2;
                max_min_len = max(size_i, size_j);
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
      for (int k = 0; k < size_j; ++k) {
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
      if (patterns_tmp.size() == N * 2) {
        break;
      }
    }

    merged_patterns.vv_patterns.clear();
    merged_patterns.vv_patterns.resize(MAX_PATTERN_LENGTH + 1);
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

class MatchedFlag {
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

class State {
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
        for (int i = 0; i < N; ++i) {
          for (int j = 0; j < N; ++j) {
            // 行に置けるか
            {
              ok = 1;
              for (int k = 0; k < len; ++k) {
                if (patterns.vv_patterns[len][pat_index].pattern[k] != grid[i][(j + k) % N] && grid[i][(j + k) % N] != CHARACTER_SIZE) {
                  ok = 0;
                  break;
                }
              }
              if (ok) {
                for (int k = 0; k < len; ++k) {
                  grid[i][(j + k) % N] = patterns.vv_patterns[len][pat_index].pattern[k];
                }
                break;
              }
            }
            // 列に置けるか
            {
              ok = 1;
              for (int k = 0; k < len; ++k) {
                if (patterns.vv_patterns[len][pat_index].pattern[k] != grid[(i + k) % N][j] && grid[(i + k) % N][j] != CHARACTER_SIZE) {
                  ok = 0;
                  break;
                }
              }
              if (ok) {
                for (int k = 0; k < len; ++k) {
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

  void recalc_all(const Patterns& patterns, bool rough = false) {
    reset_matched_flags(patterns);

    for (int i = MIN_PATTERN_LENGTH; i <= MAX_PATTERN_LENGTH; i++) {
      for (int j = 0; j < patterns.vv_patterns[i].size(); j++) {
        matched_flags[i][j].clear();
        for (int k = 0; k < N; ++k) {
          for (int l = 0; l < N; ++l) {
            matched_flags[i][j].set_flag(k, l, 0, is_matched(k, l, patterns.vv_patterns[i][j].pattern, 0));
            matched_flags[i][j].set_flag(k, l, 1, is_matched(k, l, patterns.vv_patterns[i][j].pattern, 1));
            if (rough && matched_flags[i][j].get_count() > 0) {
              break;
            }
          }
          if (rough && matched_flags[i][j].get_count() > 0) {
            break;
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
        for (int k = 0; k < i; ++k) {
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
    for (int i = MAX_PATTERN_LENGTH + 1 - 1; i >= 0; --i) {
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
    for (int i = 0; i < 6; ++i) {
      prefixMaps.push_back(vector<vector<P>>(1 << (i * 3)));
      if (i == 0) {
        continue;
      }
      for (int j = 0; j < vv.size(); ++j) {
        for (int k = 0; k < vv[j].size(); ++k) {
          int num = 0;
          for (int l = 0; l < i; ++l) {
            num *= CHARACTER_SIZE;
            num += vv[j][(k + l) % vv[j].size()];
          }
          if (num >= prefixMaps[i].size()) {
            cout << "NGa" << endl;
            continue;
          }
          prefixMaps[i][num].emplace_back(j, k);
        }
      }
    }

    int best_g[N][N];
    int best_score = -1;

    int g[N][N];
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        g[i][j] = CHARACTER_SIZE;
        best_g[i][j] = CHARACTER_SIZE;
        grid[i][j] = 0;
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

    int ii_num = vv.size();

    for (int ii = 0; ii < ii_num; ++ii) {
      int no_start_count = 0;

      for (int i = 0; i < vv.size(); ++i) {
        for (int j = 0; j < vv.size(); ++j) {
          if (i == j || ii == i || ii == j) {
            continue;
          }
          if (get_elapsed_time() > time_limit) {
            break;
          }
          for (int k = 0; k < N; ++k) {
            if (get_elapsed_time() > time_limit) {
              break;
            }
            for (int l = 0; l < N; ++l) {
              // 2行目と3行目をセット
              for (int m = 0; m < N; ++m) {
                g[0][m] = CHARACTER_SIZE;
                g[1][m] = CHARACTER_SIZE;
                g[2][m] = CHARACTER_SIZE;
              }
              for (int j = 0; j < vv[ii].size(); ++j) {
                g[0][j] = vv[ii][j];
              }
              for (int m = 0; m < vv[i].size(); ++m) {
                g[1][(k + m) % N] = vv[i][m];
              }
              for (int m = 0; m < vv[j].size(); ++m) {
                g[2][(l + m) % N] = vv[j][m];
              }

              // 各列を決定していく
              for (int m = 0; m < N; ++m) {
                decided_col[m] = -1;
              }
              decided_col_count = 0;
              for (int m = 0; m < N * 2; ++m) {
                used[m] = -1;
              }
              used_version++;
              used[ii] = used_version;
              used[i] = used_version;
              used[j] = used_version;
              int ng_count = 0;
              const int NG_SAFE = 0;
              while (true) {
                bool has_change = false;
                for (int m = 0; m < N; ++m) {
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

                    if (idx1 == -1) {
                      idx1 = pre.first;
                      idx2 = pre.second;
                    }
                    else {
                      idx1 = -2;
                      break;
                    }
                  }

                  if (idx1 == -1) {
                    ng_count++;
                    if (ng_count > NG_SAFE) {
                      break;
                    }
                  }
                  else if (idx1 >= 0) {
                    decided_col[m] = idx1;
                    decided_col_2[m] = idx2;
                    decided_col_count++;
                    used[idx1] = used_version;
                    has_change = true;
                    break;
                  }
                }

                if (ng_count > NG_SAFE) {
                  break;
                }
                if (!has_change) {
                  break;
                }
              }

              // 失敗したら元に戻す
              if (ng_count > NG_SAFE) {
                for (int m = 0; m < vv[i].size(); ++m) {
                  g[1][(k + m) % N] = CHARACTER_SIZE;
                }
                for (int m = 0; m < vv[j].size(); ++m) {
                  g[2][(l + m) % N] = CHARACTER_SIZE;
                }
                used_version++;
                continue;
              }

              // 行決定に使う列を決定する
              const int COL_LENGTH = 1;
              int start_col = -1;
              for (int m = 0; m < N; ++m) {
                bool ok = true;
                for (int n = 0; n < COL_LENGTH; ++n) {
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
                for (int mm = 0; mm < COL_LENGTH; ++mm) {
                  int m = (start_col + mm) % N;
                  for (int n = 0; n < N; ++n) {
                    int idx = (decided_col_2[m] + n) % N;
                    if (idx < vv[decided_col[m]].size()) {
                      g[n][m] = vv[decided_col[m]][idx];
                    }
                    else {
                      g[n][m] = CHARACTER_SIZE;
                    }

                  }
                }
              }

              if (start_col == -1) {
                no_start_count++;
              }

              // 各行を決定していく
              for (int m = 0; m < N; ++m) {
                decided_row[m] = -1;
              }
              decided_row_count = 0;
              ng_count = start_col == -1;
              while (ng_count == 0) {
                bool has_change = false;
                for (int m = 3; m < N; ++m) {
                  if (decided_row[m] != -1) {
                    continue;
                  }

                  bool is_blank = false;
                  int num = 0;
                  for (int n = 0; n < COL_LENGTH; ++n) {
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

                    bool is_ok = true;

                    for (int n = 0; n < N; ++n) {
                      if (decided_col[n] != -1) {
                        int idx_col = (decided_col_2[n] + m) % N;
                        int idx_row = (pre.second + n + N - start_col) % N;
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

                  if (idx1 == -1) {
                    ng_count++;
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

              int score = 0;
              if (ng_count == 0) {
                for (int n = 0; n < N; ++n) {
                  if (n < 3) {
                    for (int m = 0; m < N; ++m) {
                      grid[n][m] = g[n][m];
                    }
                  }
                  else {
                    for (int m = 0; m < N; ++m) {
                      if (decided_col[m] != -1) {
                        int idx = (decided_col_2[m] + n) % N;
                        grid[n][m] = idx < vv[decided_col[m]].size() ? vv[decided_col[m]][idx] : CHARACTER_SIZE;
                      }
                      else {
                        grid[n][m] = CHARACTER_SIZE;
                      }
                    }
                  }
                }

                for (int n = 3; n < N; ++n) {
                  if (decided_row[n] != -1) {
                    for (int m = 0; m < N; ++m) {
                      int idx = (decided_row_2[n] + m) % N;
                      if (idx < vv[decided_row[n]].size()) {
                        grid[n][(start_col + m) % N] = vv[decided_row[n]][idx];
                      }
                    }
                  }
                }

                recalc_all(patterns, true);
                score = get_score(patterns);
              }

              if (score > best_score) {
                best_score = score;
                for (int n = 0; n < N; ++n) {
                  if (n < 3) {
                    for (int m = 0; m < N; ++m) {
                      best_g[n][m] = g[n][m];
                    }
                  }
                  else {
                    for (int m = 0; m < N; ++m) {
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

                for (int n = 3; n < N; ++n) {
                  if (decided_row[n] != -1) {
                    for (int m = 0; m < N; ++m) {
                      int idx = (decided_row_2[n] + m) % N;
                      if (idx < vv[decided_row[n]].size()) {
                        best_g[n][(start_col + m) % N] = vv[decided_row[n]][idx];
                      }
                    }
                  }
                }

                if (best_score == PERFECT_SCORE) {
                  break;
                }
              }

              // 後片付け
              for (int m = 0; m < vv[i].size(); ++m) {
                g[1][(k + m) % N] = CHARACTER_SIZE;
              }
              for (int m = 0; m < vv[j].size(); ++m) {
                g[2][(l + m) % N] = CHARACTER_SIZE;
              }
              if (start_col != -1) {
                for (int m = 0; m < COL_LENGTH; ++m) {
                  int m2 = (start_col + m) % N;
                  for (int n = 1; n < N; ++n) {
                    g[n][m2] = CHARACTER_SIZE;
                  }
                }
              }
            }
          }
        }
      }

      if (best_score == PERFECT_SCORE) {
        break;
      }
      if (get_elapsed_time() > time_limit) {
        break;
      }
    }

    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        grid[i][j] = best_g[i][j];
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

void input_data(int case_num, PatternsManager& patterns_manager) {
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  int _n, _m;
  vector<string> vs;
  if (exec_mode == 4) {
    int _l = 10;
    int _m = 800;
    generate_true_genome(_l, _m);
    for (int i = 0; i < _m; i++) {
      vs.push_back(true_genome_strs[i]);
    }
  }
  else if (!ifs.is_open()) {
    // 標準入力
    cin >> _n >> _m;
    for (int i = 0; i < _m; ++i) {
      string s;
      cin >> s;
      vs.push_back(s);
    }
  }
  else {
    // ファイル入力
    ifs >> _n >> _m;
    for (int i = 0; i < _m; ++i) {
      string s;
      ifs >> s;
      vs.push_back(s);
    }
  }

  patterns_manager.initialize(vs);

  int ma_len = 0;
  int mi_len = 999;
  for (auto str : vs) {
    int len = str.size();
    if (len > ma_len) {
      ma_len = len;
    }
    if (len < mi_len) {
      mi_len = len;
    }
  }
  L = (ma_len + mi_len) / 2;
}

void open_ofs(int case_num, ofstream& ofs) {
  if (exec_mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofs.open(oss.str());
  }
}

void output_data(ofstream& ofs, State& state) {
  if (exec_mode == 0) {
    // 標準出力
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
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
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
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

struct AnnealingParams {
  double start_temperature[10];
  double end_temperature;
  double score_scale;
  int operation_thresholds[10];
};

void run_simulated_annealing(AnnealingParams annealingParams, State& state, const Patterns& patterns) {
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

  if (exec_mode >= 3 && exec_mode != 4) {
    cerr << "iteration_count = " << iteration_count << endl;
  }
  if (exec_mode >= 2 && exec_mode != 4) {
    cerr << "焼きなまし処理終了 : " << get_elapsed_time() << " sec" << endl;
  }
}

ll execute_solver(AnnealingParams annealingParams, PatternsManager& patterns_manager, State& state) {
  patterns_manager.build_merge_patterns(TIME_LIMIT * 0.4, 0);
  cout << "build_merge_patterns : " << get_elapsed_time() << " sec" << endl;
  state.assemble(TIME_LIMIT * 0.5, patterns_manager.merged_patterns);
  cout << "assemble : " << get_elapsed_time() << " sec" << endl;
  state.recalc_all(patterns_manager.merged_patterns);
  cout << "recalc_all : " << get_elapsed_time() << " sec" << endl;
  ll score = state.get_score(patterns_manager.merged_patterns);
  if (score == PERFECT_SCORE) {
    return PERFECT_SCORE;
  }

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
  if (exec_mode >= 2 && exec_mode != 4) {
    cerr << "パターンマージ処理終了 : " << get_elapsed_time() << " sec" << endl;
  }

  if (exec_mode >= 3 && exec_mode != 4) {
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
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
  if (exec_mode >= 2 && exec_mode != 4) {
    cerr << "完全復元処理終了 : " << get_elapsed_time() << " sec" << endl;
  }

  state.greedy_after_assemble(patterns_manager.merged_patterns);

  state.greedy_after_assemble(patterns_manager.initial_patterns);

  state.generate_random_empty();
  state.recalc_all(patterns_manager.initial_patterns);
  if (state.get_score(patterns_manager.merged_patterns) == PERFECT_SCORE) {
    return PERFECT_SCORE;
  }

  // 焼きなまし実行
  run_simulated_annealing(annealingParams, state, patterns_manager.initial_patterns);

  state.recalc_all(patterns_manager.initial_patterns);
  return state.get_score(patterns_manager.initial_patterns);
}

ll solve_case(int case_num, AnnealingParams annealingParams) {
  start_timer();

  PatternsManager patterns_manager;

  input_data(case_num, patterns_manager);

  State state(patterns_manager);

  ll score = execute_solver(annealingParams, patterns_manager, state);

  // 解答を出力
  ofstream ofs;
  open_ofs(case_num, ofs);
  output_data(ofs, state);

  if (ofs.is_open()) {
    ofs.close();
  }

  return score;
}

int main() {
  exec_mode = 3;

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
    for (int i = 0; i < 10; ++i) {
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
