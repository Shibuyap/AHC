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
#define srep(i, s, t) for (int i = (s); i < (t); ++i)
#define drep(i, n) for (int i = (n) - 1; i >= 0; --i)

using namespace std;

typedef long long int ll;
typedef pair<int, int> PAIR;

// -------------------- constants --------------------
const int INF = 1001001001;
const int DIR4_X[4] = { -1, 0, 1, 0 };
const int DIR4_Y[4] = { 0, -1, 0, 1 };
const int F_ARRAY_SIZE = 10000;
const double TIME_LIMIT_SEC = 3.8;
const int SET_SIZE = 20;

// -------------------- random helpers --------------------
namespace /* random helpers */
{
  static uint32_t rand_xorshift()
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

  static double rand_unit() { return (rand_xorshift() + 0.5) * (1.0 / UINT_MAX); }
} // namespace

// -------------------- global state --------------------
int io_mode;                           // 0: interact, 1: file
std::mt19937 rng_engine;               // RNG engine
std::normal_distribution<> normal_dist; // N(0,S)

constexpr int MAX_N = 100;
constexpr int MAX_L = 50;

// input
int board_size;
int wormhole_count;
int square_area;
int effective_size;
int sqrt_s_minus1;
int wormhole_y[MAX_N];
int wormhole_x[MAX_N];
int wormhole_perm[MAX_N];
int noise_table[F_ARRAY_SIZE];

// grid & answer tracking
int grid[MAX_L][MAX_L];
int guessed_perm[MAX_N];
int answer_perm[MAX_N];

// annealing helpers
int SEARCH_WINDOW_SIZE = 9;
int sa_measurement_noise_cache[MAX_N][SET_SIZE][11][11];

ll max_score;
ll max_layout_cost;
ll max_measurement_cost;
int candidate_perm_set[SET_SIZE][MAX_N];
int diff_cache[SET_SIZE][MAX_N][MAX_N];

// rollback / best snapshot
ll best_score;
ll best_layout_cost;
ll best_measurement_cost;
int best_grid[MAX_L][MAX_L];
int best_candidate_perm_set[SET_SIZE][MAX_N];
int best_diff_cache[SET_SIZE][MAX_N][MAX_N];

// runtime stats
int method_stat[20][10];
int measurement_count;
ll measure_cost;

// helpers for window tests
int is_within_window_iy;
int is_within_window_ix;

// -------------------- utility functions --------------------
void init_normal_dist()
{
  std::normal_distribution<>::param_type param(0.0, square_area);
  normal_dist.param(param);
}

inline int sample_normal_int() {
  return round(normal_dist(rng_engine));
}

void init_counters()
{
  measurement_count = 0;
  rep(i, 20) {
    rep(j, 10) {
      method_stat[i][j] = 0;
    }
  }
}

void load_problem(int problem_num)
{
  string file_name = "./in/";
  string num = to_string(10000 + problem_num).substr(1);
  file_name += num + ".txt";

  ifstream ifs(file_name);

  if (!ifs.is_open()) {
    cin >> board_size >> wormhole_count >> square_area;
    rep(i, wormhole_count) {
      cin >> wormhole_y[i] >> wormhole_x[i];
    }
  }
  else {
    ifs >> board_size >> wormhole_count >> square_area;
    rep(i, wormhole_count) {
      ifs >> wormhole_y[i] >> wormhole_x[i];
    }
    rep(i, wormhole_count) {
      ifs >> wormhole_perm[i];
    }
    rep(i, F_ARRAY_SIZE) {
      ifs >> noise_table[i];
    }
  }

  effective_size = board_size - 10;
  sqrt_s_minus1 = 0;
  srep(i, 1, 31) {
    if (i * i == square_area) {
      sqrt_s_minus1 = i - 1;
      break;
    }
  }

  init_normal_dist();
}

void open_output_stream(int prob_num, ofstream& ofs)
{
  if (io_mode == 0) {
    return;
  }
  string file_name = "./out/" + to_string(10000 + prob_num).substr(1) + ".txt";
  ofs.open(file_name);
}

ll layout_cost()
{
  ll cost = 0;
  rep(i, board_size) {
    rep(j, board_size) {
      cost += (grid[i][j] - grid[(i + 1) % board_size][j]) * (grid[i][j] - grid[(i + 1) % board_size][j]);
      cost += (grid[i][j] - grid[i][(j + 1) % board_size]) * (grid[i][j] - grid[i][(j + 1) % board_size]);
    }
  }
  return cost;
}

void init_layout_random_binary() {
  rep(i, board_size) {
    rep(j, board_size) {
      grid[i][j] = rand_xorshift() % 2 * 1000;
    }
  }
}

void output_layout(ofstream& ofs)
{
  auto& out = (io_mode == 0 ? cout : ofs);
  rep(i, board_size) {
    rep(j, board_size) {
      out << grid[i][j] << ' ';
    }
    out << '\n';
  }
  if (io_mode == 0) {
    fflush(stdout);
  }
}

// -------- measurement wrappers --------
int measure_cell(int i, int y, int x, ofstream& ofs)
{
  if (io_mode == 0) {
    cout << i << ' ' << y << ' ' << x << '\n';
    fflush(stdout);
  }
  else {
    ofs << i << ' ' << y << ' ' << x << '\n';
  }

  int m = 0;
  if (io_mode == 0) {
    cin >> m;
  }
  else {
    int yy = (wormhole_y[wormhole_perm[i]] + y + board_size * 10) % board_size;
    int xx = (wormhole_x[wormhole_perm[i]] + x + board_size * 10) % board_size;
    m = std::max(0, std::min(1000, (int)round(grid[yy][xx] + noise_table[measurement_count])));
    measurement_count++;
  }

  measure_cost += 100LL * (10LL + abs(y) + abs(x));
  return m;
}

int sa_measure_cell(int i, int y, int x, int se)
{
  int slide = (SEARCH_WINDOW_SIZE - 1) / 2;
  int yy = (wormhole_y[answer_perm[i]] + y + board_size * 10) % board_size;
  int xx = (wormhole_x[answer_perm[i]] + x + board_size * 10) % board_size;
  return std::max(0, std::min(1000, (int)round(grid[yy][xx] + sa_measurement_noise_cache[i][se][y + slide][x + slide])));
}

inline int sa_measure_value(int i, int y, int x, int se, int value)
{
  int slide = (SEARCH_WINDOW_SIZE - 1) / 2;
  return std::max(0, std::min(1000, (int)round(value + sa_measurement_noise_cache[i][se][y + slide][x + slide])));
}

void output_measurements(ofstream& ofs)
{
  int slide = (SEARCH_WINDOW_SIZE - 1) / 2;
  rep(i, wormhole_count) {
    int scores[31][31];
    srep(j, -slide, -slide + SEARCH_WINDOW_SIZE) {
      srep(k, -slide, -slide + SEARCH_WINDOW_SIZE) {
        scores[j + slide][k + slide] = measure_cell(i, j, k, ofs);
      }
    }

    int minDiff = INF;
    rep(j, wormhole_count) {
      int sumDiff = 0;
      srep(k, -slide, -slide + SEARCH_WINDOW_SIZE) {
        srep(l, -slide, -slide + SEARCH_WINDOW_SIZE) {
          int y = (wormhole_y[j] + k) % board_size;
          int x = (wormhole_x[j] + l) % board_size;
          sumDiff += abs(grid[y][x] - scores[k + slide][l + slide]);
        }
      }
      if (sumDiff < minDiff) {
        minDiff = sumDiff;
        guessed_perm[i] = j;
      }
    }
  }
}

void output_answer(ofstream& ofs)
{
  auto& out = (io_mode == 0 ? cout : ofs);
  out << "-1 -1 -1\n";
  rep(i, wormhole_count) {
    out << guessed_perm[i] << '\n';
  }
  if (io_mode == 0) {
    fflush(stdout);
  }
}

ll calc_offline_score()
{
  double score = 1e14;
  rep(i, wormhole_count) {
    if (guessed_perm[i] != wormhole_perm[i]) {
      score *= 0.8;
    }
  }
  score = score / (1e5 + layout_cost() + measure_cost);
  return ceil(score);
}

void snapshot_best_state()
{
  best_score = max_score;
  best_layout_cost = max_layout_cost;
  best_measurement_cost = max_measurement_cost;
  rep(i, board_size) {
    rep(j, board_size) {
      best_grid[i][j] = grid[i][j];
    }
  }
  rep(s, SET_SIZE) {
    rep(i, wormhole_count) {
      best_candidate_perm_set[s][i] = candidate_perm_set[s][i];
    }
    rep(i, wormhole_count) {
      rep(j, wormhole_count) {
        best_diff_cache[s][i][j] = diff_cache[s][i][j];
      }
    }
  }
}

void restore_best_state()
{
  max_score = best_score;
  max_layout_cost = best_layout_cost;
  max_measurement_cost = best_measurement_cost;
  rep(i, board_size) {
    rep(j, board_size) {
      grid[i][j] = best_grid[i][j];
    }
  }
  rep(s, SET_SIZE) {
    rep(i, wormhole_count) {
      candidate_perm_set[s][i] = best_candidate_perm_set[s][i];
    }
    rep(i, wormhole_count) {
      rep(j, wormhole_count) {
        diff_cache[s][i][j] = best_diff_cache[s][i][j];
      }
    }
  }
}

int scores_sa_measurements[SET_SIZE][31][31];
void init_sa_measurements()
{
  max_measurement_cost = 0;
  int slide = (SEARCH_WINDOW_SIZE - 1) / 2;
  rep(i, wormhole_count) {
    srep(k, -slide, -slide + SEARCH_WINDOW_SIZE) {
      srep(l, -slide, -slide + SEARCH_WINDOW_SIZE) {
        rep(m, SET_SIZE) {
          scores_sa_measurements[m][k + slide][l + slide] = sa_measure_cell(i, k, l, m);
        }
        max_measurement_cost += 100LL * (10LL + abs(k) + abs(l));
      }
    }
    rep(m, SET_SIZE) {
      int minDiff = INF;
      rep(j, wormhole_count) {
        diff_cache[m][i][j] = 0;
        srep(k, -slide, -slide + SEARCH_WINDOW_SIZE) {
          srep(l, -slide, -slide + SEARCH_WINDOW_SIZE) {
            int y = (wormhole_y[j] + k + board_size) % board_size;
            int x = (wormhole_x[j] + l + board_size) % board_size;
            diff_cache[m][i][j] += abs(grid[y][x] - scores_sa_measurements[m][k + slide][l + slide]);
          }
        }
        if (diff_cache[m][i][j] < minDiff) {
          minDiff = diff_cache[m][i][j];
          candidate_perm_set[m][i] = j;
        }
      }
    }
  }
}

void init_simulated_annealing()
{
  rep(i, wormhole_count) {
    answer_perm[i] = i;
  }
  init_normal_dist();
  rep(i, wormhole_count) {
    rep(j, SET_SIZE) {
      rep(k, SEARCH_WINDOW_SIZE) {
        rep(l, SEARCH_WINDOW_SIZE) {
          sa_measurement_noise_cache[i][j][k][l] = sample_normal_int();
        }
      }
    }
  }

  max_layout_cost = layout_cost();
  init_sa_measurements();

  max_score = 0;
  rep(m, SET_SIZE) {
    double score = 1e14;
    rep(i, wormhole_count) {
      if (candidate_perm_set[m][i] != answer_perm[i]) {
        score *= 0.8;
      }
    }
    score /= (1e5 + max_layout_cost + max_measurement_cost);
    max_score += ceil(score);
  }
  snapshot_best_state();
}

bool is_within_window(int i, int y, int x)
{
  int slide = (SEARCH_WINDOW_SIZE - 1) / 2;
  int iy = wormhole_y[i];
  int ix = wormhole_x[i];
  int iU = iy - slide;
  int iD = iU + SEARCH_WINDOW_SIZE;
  int iL = ix - slide;
  int iR = iL + SEARCH_WINDOW_SIZE;
  srep(j, -1, 2) {
    srep(k, -1, 2) {
      int yy = y + board_size * j;
      int xx = x + board_size * k;
      if (iU <= yy && yy < iD && iL <= xx && xx < iR) {
        is_within_window_iy = yy - iU - slide;
        is_within_window_ix = xx - iL - slide;
        return true;
      }
    }
  }
  return false;
}

ll calc_delta_layout_cost(int y, int x, int beforeP, int afterP)
{
  ll ret = 0;
  rep(i, 4) {
    int ny = (y + DIR4_Y[i] + board_size) % board_size;
    int nx = (x + DIR4_X[i] + board_size) % board_size;
    ret -= (beforeP - grid[ny][nx]) * (beforeP - grid[ny][nx]);
    ret += (afterP - grid[ny][nx]) * (afterP - grid[ny][nx]);
  }
  return ret;
}

// ----- tweak method 1: single cell change -----
int tweak_single_cell_vector[MAX_N];
void tweak_single_cell(double temperature)
{
  method_stat[1][0]++;

  int y = rand_xorshift() % board_size;
  int x = rand_xorshift() % board_size;
  int diff = rand() % 201 - 100;
  int beforeP = grid[y][x];
  int afterP = max(0, min(1000, grid[y][x] + diff));

  int cnt = 0;

  rep(i, wormhole_count) {
    if (!is_within_window(i, y, x)) {
      continue;
    }
    int iy = is_within_window_iy;
    int ix = is_within_window_ix;

    rep(m, SET_SIZE) {
      rep(j, wormhole_count) {
        diff_cache[m][i][j] -= abs(grid[(wormhole_y[j] + iy + board_size) % board_size][(wormhole_x[j] + ix + board_size) % board_size] - sa_measure_value(i, iy, ix, m, beforeP));
      }
    }

    grid[y][x] = afterP;

    rep(m, SET_SIZE) {
      rep(j, wormhole_count) {
        diff_cache[m][i][j] += abs(grid[(wormhole_y[j] + iy + board_size) % board_size][(wormhole_x[j] + ix + board_size) % board_size] - sa_measure_value(i, iy, ix, m, afterP));
      }
    }
    grid[y][x] = beforeP;

    rep(m, SET_SIZE) {
      int minDiff = INF;
      rep(j, wormhole_count) {
        if (diff_cache[m][i][j] < minDiff) {
          minDiff = diff_cache[m][i][j]; candidate_perm_set[m][i] = j;
        }
      }
    }
    tweak_single_cell_vector[cnt++] = i;
  }

  if (cnt == 0) {
    return;
  }
  method_stat[1][1]++;

  grid[y][x] = afterP;

  ll delta_layout = calc_delta_layout_cost(y, x, beforeP, afterP);
  ll delta_measure = 0; // simplified (kept as 0)

  ll tmp_total = 0;
  rep(m, SET_SIZE) {
    double score = 1e14;
    rep(i, wormhole_count) {
      if (candidate_perm_set[m][i] != answer_perm[i]) {
        score *= 0.8;
      }
    }
    score /= (1e5 + max_layout_cost + delta_layout + max_measurement_cost + delta_measure);
    tmp_total += ceil(score);
  }

  ll diffScore = tmp_total - max_score;
  if (exp((double)diffScore / temperature) > rand_unit()) {
    method_stat[1][2]++;
    max_score += diffScore;
    max_layout_cost += delta_layout;
    max_measurement_cost += delta_measure;
    if (max_score > best_score) {
      method_stat[1][3]++;
      if (method_stat[1][3] % 100 == 0) {
        cout << max_score << ' ' << best_score << '\n';
      }
      snapshot_best_state();
    }
  }
  else {
    grid[y][x] = beforeP; // rollback grid
    // rollback diff_cache and candidate_perm_set
    rep(ii, cnt) {
      int i_idx = tweak_single_cell_vector[ii];
      // recompute minimal EE for i_idx
      rep(m, SET_SIZE) {
        int minDiff = INF;
        rep(j, wormhole_count) {
          diff_cache[m][i_idx][j] = 0;
          int slide = (SEARCH_WINDOW_SIZE - 1) / 2;
          srep(k, -slide, -slide + SEARCH_WINDOW_SIZE) {
            srep(l, -slide, -slide + SEARCH_WINDOW_SIZE) {
              int y_ = (wormhole_y[j] + k + board_size) % board_size;
              int x_ = (wormhole_x[j] + l + board_size) % board_size;
              diff_cache[m][i_idx][j] += abs(grid[y_][x_] - sa_measure_cell(i_idx, k, l, m));
            }
          }
          if (diff_cache[m][i_idx][j] < minDiff) {
            minDiff = diff_cache[m][i_idx][j];
            candidate_perm_set[m][i_idx] = j;
          }
        }
      }
    }
  }
}

// -------------- main solver per instance --------------
ll solve_problem(int prob_num)
{
  init_counters();
  ofstream ofs;
  load_problem(prob_num);
  open_output_stream(prob_num, ofs);

  init_layout_random_binary(); // choose one of the initial layouts
  init_simulated_annealing();

  clock_t start_clock = clock();
  int loopCount = 0;
  while (true) {
    double elapsed = (double)(clock() - start_clock) / CLOCKS_PER_SEC;
    if (elapsed > TIME_LIMIT_SEC) {
      break;
    }
    double progress = elapsed / TIME_LIMIT_SEC;
    double temperature = 2000 * (1.0 - progress);
    tweak_single_cell(temperature);
    ++loopCount;
  }

  restore_best_state();

  output_layout(ofs);
  output_measurements(ofs);
  output_answer(ofs);

  if (ofs.is_open()) {
    ofs.close();
  }
  if (io_mode == 0) {
    return 0;
  }
  return calc_offline_score();
}

// ----------------------------- main -----------------------------
int main_old()
{
  io_mode = 1; // 0: interactive, 1: file batch

  if (io_mode == 0) {
    solve_problem(0);
  }
  else if (io_mode == 1) {
    srep(i, 0, 10) {
      ll score = solve_problem(i);
      cout << "num = " << i << ", score = " << score << '\n';
    }
  }
  return 0;
}
