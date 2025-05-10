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

// タイマー
namespace
{
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
namespace
{
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

const double TIME_LIMIT = 2.9;
int exec_mode;

class Board
{
public:
  static constexpr int MIN_LINE = 5;
  static constexpr int MAX_LINE = 50;

  int n;
  int a[11];
  int a_sum;
  vector<vector<int>> x, y;

  Board() : n(0), x(20000), y(20000) {
  }

  void init() {
    for (int i = 0; i < 20000; i++) {
      sort(x[i].begin(), x[i].end());
      sort(y[i].begin(), y[i].end());
    }

    a_sum = 0;
    for (int i = 1; i <= 10; i++) {
      a_sum += a[i];
    }
  }
};

class Answer
{
public:
  static constexpr int MIN = 0;
  static constexpr int MAX = 20000;

  vector<int> xs, ys; // 半開区間 : [xs[i], xs[i+1])

  vector<vector<int>> counts;
  int b[11];

  void initialize(int v_num, int h_num) {
    xs.clear();
    xs.push_back(MIN);
    for (int i = 1; i < v_num + 1; i++) {
      xs.push_back(MIN + (MAX - MIN) * i / (v_num + 1));
    }
    xs.push_back(MAX);

    ys.clear();
    ys.push_back(MIN);
    for (int i = 1; i < h_num + 1; i++) {
      ys.push_back(MIN + (MAX - MIN) * i / (h_num + 1));
    }
    ys.push_back(MAX);
  }

  void calc_counts(const Board& board) {
    for (int i = 0; i < 11; i++) {
      b[i] = 0;
    }

    counts.resize(xs.size() - 1, vector<int>(ys.size() - 1));
    for (int i = 0; i < xs.size() - 1; i++) {
      int sx = xs[i];
      int gx = xs[i + 1];
      for (int j = 0; j < ys.size() - 1; j++) {
        int sy = ys[j];
        int gy = ys[j + 1];

        counts[i][j] = 0;
        for (int k = sx; k < gx; k++) {
          for (auto l : board.x[k]) {
            if (sy <= l && l < gy) {
              counts[i][j]++;
            }
          }
        }

        if (1 <= counts[i][j] && counts[i][j] <= 10) {
          b[counts[i][j]]++;
        }
      }
    }
  }

  void update_xs(const Board& board, int num, int diff) {
    if (diff < 0) {
      for (int i = xs[num] + diff; i < xs[num]; ++i) {
        for (auto y : board.x[i]) {
          for (int j = 0; j < ys.size() - 1; j++) {
            if (ys[j] <= y && y < ys[j + 1]) {
              // 上側を減らす
              if (1 <= counts[num - 1][j] && counts[num - 1][j] <= 10) {
                b[counts[num - 1][j]]--;
              }
              counts[num - 1][j]--;
              if (1 <= counts[num - 1][j] && counts[num - 1][j] <= 10) {
                b[counts[num - 1][j]]++;
              }

              // 下側を増やす
              if (1 <= counts[num][j] && counts[num][j] <= 10) {
                b[counts[num][j]]--;
              }
              counts[num][j]++;
              if (1 <= counts[num][j] && counts[num][j] <= 10) {
                b[counts[num][j]]++;
              }
            }
          }
        }
      }
    }
    else {
      for (int i = xs[num]; i <= xs[num] + diff - 1; ++i) {
        for (auto y : board.x[i]) {
          for (int j = 0; j < ys.size() - 1; j++) {
            if (ys[j] <= y && y < ys[j + 1]) {
              // 上側を増やす
              if (1 <= counts[num - 1][j] && counts[num - 1][j] <= 10) {
                b[counts[num - 1][j]]--;
              }
              counts[num - 1][j]++;
              if (1 <= counts[num - 1][j] && counts[num - 1][j] <= 10) {
                b[counts[num - 1][j]]++;
              }

              // 下側を減らす
              if (1 <= counts[num][j] && counts[num][j] <= 10) {
                b[counts[num][j]]--;
              }
              counts[num][j]--;
              if (1 <= counts[num][j] && counts[num][j] <= 10) {
                b[counts[num][j]]++;
              }
            }
          }
        }
      }
    }

    xs[num] += diff;
  }

  void update_ys(const Board& board, int num, int diff) {
    if (diff < 0) {
      for (int i = ys[num] + diff; i < ys[num]; ++i) {
        for (auto x : board.y[i]) {
          for (int j = 0; j < xs.size() - 1; j++) {
            if (xs[j] <= x && x < xs[j + 1]) {
              // 左側を減らす
              if (1 <= counts[j][num - 1] && counts[j][num - 1] <= 10) {
                b[counts[j][num - 1]]--;
              }
              counts[j][num - 1]--;
              if (1 <= counts[j][num - 1] && counts[j][num - 1] <= 10) {
                b[counts[j][num - 1]]++;
              }
              // 右側を増やす
              if (1 <= counts[j][num] && counts[j][num] <= 10) {
                b[counts[j][num]]--;
              }
              counts[j][num]++;
              if (1 <= counts[j][num] && counts[j][num] <= 10) {
                b[counts[j][num]]++;
              }
            }
          }
        }
      }
    }
    else {
      for (int i = ys[num]; i <= ys[num] + diff - 1; ++i) {
        for (auto x : board.y[i]) {
          for (int j = 0; j < xs.size() - 1; j++) {
            if (xs[j] <= x && x < xs[j + 1]) {
              // 左側を増やす
              if (1 <= counts[j][num - 1] && counts[j][num - 1] <= 10) {
                b[counts[j][num - 1]]--;
              }
              counts[j][num - 1]++;
              if (1 <= counts[j][num - 1] && counts[j][num - 1] <= 10) {
                b[counts[j][num - 1]]++;
              }
              // 右側を減らす
              if (1 <= counts[j][num] && counts[j][num] <= 10) {
                b[counts[j][num]]--;
              }
              counts[j][num]--;
              if (1 <= counts[j][num] && counts[j][num] <= 10) {
                b[counts[j][num]]++;
              }
            }
          }
        }
      }
    }

    ys[num] += diff;
  }
};

Board input_data(int case_num) {
  Board board;

  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    int  _k;
    cin >> board.n >> _k;
    for (int i = 1; i <= 10; i++) {
      cin >> board.a[i];
    }
    for (int i = 0; i < board.n; i++) {
      int x, y;
      cin >> x >> y;
      x += 10000;
      y += 10000;
      board.x[x].push_back(y);
      board.y[y].push_back(x);
    }
  }
  else {
    // ファイル入力
    int _k;
    ifs >> board.n >> _k;
    for (int i = 1; i <= 10; i++) {
      ifs >> board.a[i];
    }
    for (int i = 0; i < board.n; i++) {
      int x, y;
      ifs >> x >> y;
      x += 10000;
      y += 10000;
      board.x[x].push_back(y);
      board.y[y].push_back(x);
    }
  }

  board.init();

  return board;
}

int calculate_score(const Board& board, const Answer& answer) {
  int ok_cnt = 0;
  for (int i = 1; i <= 10; i++) {
    ok_cnt += min(answer.b[i], board.a[i]);
  }

  int res = round(1e6 * ok_cnt / board.a_sum);
  return res;
}

void output_data(int case_num, const Answer& answer) {
  if (exec_mode == 0) {
    // 標準出力
    cout << answer.xs.size() + answer.ys.size() - 4 << endl;
    for (int i = 1; i < answer.xs.size() - 1; i++) {
      cout << answer.xs[i] - 10000 - 1 << " " << -10000000 << " " << answer.xs[i] - 10000 << " " << 10000000 << endl;
    }
    for (int i = 1; i < answer.ys.size() - 1; i++) {
      cout << -10000000 << " " << answer.ys[i] - 10000 - 1 << " " << 10000000 << " " << answer.ys[i] - 10000 << endl;
    }
  }
  else {
    // ファイル出力
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofstream ofs(oss.str());

    ofs << answer.xs.size() + answer.ys.size() - 4 << endl;
    for (int i = 1; i < answer.xs.size() - 1; i++) {
      ofs << answer.xs[i] - 10000 - 1 << " " << -10000000 << " " << answer.xs[i] - 10000 << " " << 10000000 << endl;
    }
    for (int i = 1; i < answer.ys.size() - 1; i++) {
      ofs << -10000000 << " " << answer.ys[i] - 10000 - 1 << " " << 10000000 << " " << answer.ys[i] - 10000 << endl;
    }

    if (ofs.is_open()) {
      ofs.close();
    }
  }
}

Answer build_initial_solution(const Board& board) {
  Answer best_answer;
  int best_score = -1;

  int iter = 0;
  while (get_elapsed_time() < TIME_LIMIT / 3) {
    iter++;
    int sum = rand_xorshift() % 71 + 30;
    int v_num = rand_xorshift() % (sum - 1) + 1;
    int h_num = sum - v_num;
    Answer answer;
    answer.initialize(v_num, h_num);
    answer.calc_counts(board);
    int score = calculate_score(board, answer);
    if (score > best_score) {
      best_score = score;
      best_answer = answer;
    }
  }
  cerr << "iter = " << iter << ", score = " << best_score << endl;

  return best_answer;
}

struct AnnealingParams
{
  double start_temperature;
  double end_temperature;
  double score_scale;
};

void annealing(const Board& board, Answer& answer, const AnnealingParams& params) {
  Answer best_answer = answer;
  answer.calc_counts(board);
  int best_score = calculate_score(board, answer);

  double start_time = get_elapsed_time();
  double now_time = start_time;
  int iter = 0;
  while (true) {
    iter++;

    if (iter % 100 == 0) {
      now_time = get_elapsed_time();
      if (now_time > TIME_LIMIT) {
        break;
      }
    }

    int choice = rand_xorshift() % 2;
    int num = 1;
    if (choice == 0) {
      num = rand_xorshift() % (answer.xs.size() - 2) + 1;
    }
    else {
      num = rand_xorshift() % (answer.ys.size() - 2) + 1;
    }
    int diff = rand_xorshift() % 61 - 30;
    while (diff == 0) {
      diff = rand_xorshift() % 61 - 30;
    }

    int current_score = calculate_score(board, answer);
    if (choice == 0) {
      if (answer.xs[num] + diff <= answer.xs[num - 1]) {
        continue;
      }
      if (answer.xs[num] + diff >= answer.xs[num + 1]) {
        continue;
      }
      answer.update_xs(board, num, diff);
    }
    else {
      if (answer.ys[num] + diff <= answer.ys[num - 1]) {
        continue;
      }
      if (answer.ys[num] + diff >= answer.ys[num + 1]) {
        continue;
      }
      answer.update_ys(board, num, diff);
    }

    double progress_ratio = (now_time - start_time) / (TIME_LIMIT - start_time);
    double temp = params.start_temperature + (params.end_temperature - params.start_temperature) * progress_ratio;

    int new_score = calculate_score(board, answer);
    double diff_score = (new_score - current_score) * params.score_scale;
    double prob = exp(diff_score / temp);
    if (prob > rand_01()) {
      // 採用

      if (new_score > best_score) {
        best_score = new_score;
        best_answer = answer;
      }
    }
    else {
      if (choice == 0) {
        answer.update_xs(board, num, -diff);
      }
      else {
        answer.update_ys(board, num, -diff);
      }
    }
  }

  cerr << "iter = " << iter << ", score = " << best_score << endl;

  answer = best_answer;
}

ll solve_case(int case_num) {
  start_timer();

  Board board = input_data(case_num);

  Answer answer = build_initial_solution(board);

  AnnealingParams params;
  params.start_temperature = 10000048.0;
  params.end_temperature = 0.0;
  params.score_scale = 12345.0;
  annealing(board, answer, params);

  output_data(case_num, answer);

  int score = 0;
  if (exec_mode != 0) {
    score = calculate_score(board, answer);
    cerr << answer.xs.size() + answer.ys.size() - 4 << endl;
  }
  return score;
}

int main() {
  exec_mode = 2;

  if (exec_mode == 0) {
    solve_case(0);
  }
  else if (exec_mode <= 2) {
    int sum_score = 0;
    for (int i = 0; i < 15; i++) {
      int score = solve_case(i);
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
