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

  vector<int> xs, ys;

  void initialize(int v_num, int h_num) {
    xs.clear();
    xs.push_back(MIN);
    for (int i = 1; i <= v_num + 1; i++) {
      xs.push_back(MIN + (MAX - MIN) * i / (v_num + 1));
    }
    xs.push_back(MAX);

    ys.clear();
    ys.push_back(MIN);
    for (int i = 1; i <= h_num + 1; i++) {
      ys.push_back(MIN + (MAX - MIN) * i / (h_num + 1));
    }
    ys.push_back(MAX);
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
  int b[11] = {};
  for (int i = 0; i < answer.xs.size() - 1; i++) {
    int sx = answer.xs[i];
    int gx = answer.xs[i + 1];
    for (int j = 0; j < answer.ys.size() - 1; j++) {
      int sy = answer.ys[j];
      int gy = answer.ys[j + 1];
      int cnt = 0;
      for (int k = sx; k < gx; k++) {
        for (auto l : board.x[k]) {
          if (sy <= l && l < gy) {
            cnt++;
          }
        }
        if (cnt > 10) {
          break;
        }
      }
      if (1 <= cnt && cnt <= 10) {
        b[cnt]++;
      }
    }
  }

  int ok_cnt = 0;
  for (int i = 1; i <= 10; i++) {
    ok_cnt += min(b[i], board.a[i]);
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
    int sum = rand_xorshift() % 70 + 30;
    int v_num = rand_xorshift() % (sum - 1) + 1;
    int h_num = sum - v_num;
    Answer answer;
    answer.initialize(v_num, h_num);
    int score = calculate_score(board, answer);
    if (score > best_score) {
      best_score = score;
      best_answer = answer;
    }
  }
  cerr << "iter = " << iter << ", score = " << best_score << endl;

  return best_answer;
}

void annealing(const Board& board, Answer& answer) {
  Answer best_answer = answer;
  int best_score = calculate_score(board, answer);

  int iter = 0;
  while (get_elapsed_time() < TIME_LIMIT) {
    iter++;

    int choice = rand_xorshift() % 2;
    int num = 1;
    if (choice == 0) {
      num = rand_xorshift() % (answer.xs.size() - 2) + 1;
    }
    else {
      num = rand_xorshift() % (answer.ys.size() - 2) + 1;
    }
    int diff = rand_xorshift() % 61 - 30;

    if (choice == 0) {
      if (answer.xs[num] + diff <= answer.xs[num - 1]) {
        continue;
      }
      if (answer.xs[num] + diff >= answer.xs[num + 1]) {
        continue;
      }
      answer.xs[num] += diff;
    }
    else {
      if (answer.ys[num] + diff <= answer.ys[num - 1]) {
        continue;
      }
      if (answer.ys[num] + diff >= answer.ys[num + 1]) {
        continue;
      }
      answer.ys[num] += diff;
    }

    int new_score = calculate_score(board, answer);
    if (new_score >= best_score) {
      best_score = new_score;
      best_answer = answer;
    }
    else {
      if (choice == 0) {
        answer.xs[num] -= diff;
      }
      else {
        answer.ys[num] -= diff;
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

  annealing(board, answer);

  output_data(case_num, answer);

  int score = 0;
  if (exec_mode != 0) {
    score = calculate_score(board, answer);
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
