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

const double TIME_LIMIT = 2.8;
int exec_mode;

const int N = 20;
const int H = 1000;

class Input {
public:
  int n, k, h, t, d;
  vector<vector<double>> owns;
  vector<vector<double>> targets;
};

class Board {
public:
  int n;
  vector<vector<int>> v;
  vector<vector<int>> h;

  Board(int _n) : n(_n) {
    v.resize(n);
    for (int i = 0; i < n; i++) {
      v[i].resize(n - 1, 0);
    }
    h.resize(n - 1);
    for (int i = 0; i < n - 1; i++) {
      h[i].resize(n, 0);
    }
  }
};

class Answer {
public:
  Board initial_board;
  Board board;

  int max_t;
  int t;
  vector<vector<int>> turns;

  Answer(int _n, int _max_t) : max_t(_max_t), initial_board(_n), board(_n) {
    t = 0;
    for (int i = 0; i < max_t; i++) {
      turns.push_back(vector<int>(5, 0));
    }
  }

  void add_turn_1(int x, int y, int k) {
    if (t >= max_t) {
      cerr << "Error: add_turn_1 called after max_t reached." << endl;
      return;
    }

    turns[t][0] = 1;
    turns[t][1] = x;
    turns[t][2] = y;
    turns[t][3] = k;
    turns[t][4] = 0; // dummy value
    t++;
  }

  void add_turn_2(int x, int y) {
    if (t >= max_t) {
      cerr << "Error: add_turn_2 called after max_t reached." << endl;
      return;
    }
    turns[t][0] = 2;
    turns[t][1] = x;
    turns[t][2] = y;
    turns[t][3] = 0; // dummy value
    turns[t][4] = 0; // dummy value
    t++;
  }

  void add_turn_3(int x, int y) {
    if (t >= max_t) {
      cerr << "Error: add_turn_3 called after max_t reached." << endl;
      return;
    }
    turns[t][0] = 3;
    turns[t][1] = x;
    turns[t][2] = y;
    turns[t][3] = 0; // dummy value
    turns[t][4] = 0; // dummy value
    t++;
  }

  void add_turn_4(int x1, int y1, int x2, int y2) {
    if (t >= max_t) {
      cerr << "Error: add_turn_4 called after max_t reached." << endl;
      return;
    }
    turns[t][0] = 4;
    turns[t][1] = x1;
    turns[t][2] = y1;
    turns[t][3] = x2;
    turns[t][4] = y2;
    t++;
  }
};

Input input_data(int case_num) {
  Input input;

  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    cin >> input.n >> input.k >> input.h >> input.t >> input.d;
    for (int i = 0; i < input.k; i++) {
      vector<double> own(3);
      for (int j = 0; j < 3; j++) {
        cin >> own[j];
      }
      input.owns.push_back(own);
    }
    for (int i = 0; i < input.h; i++) {
      vector<double> target(3);
      for (int j = 0; j < 3; j++) {
        cin >> target[j];
      }
      input.targets.push_back(target);
    }
  }
  else {
    // ファイル入力
    ifs >> input.n >> input.k >> input.h >> input.t >> input.d;
    for (int i = 0; i < input.k; i++) {
      vector<double> own(3);
      for (int j = 0; j < 3; j++) {
        ifs >> own[j];
      }
      input.owns.push_back(own);
    }
    for (int i = 0; i < input.h; i++) {
      vector<double> target(3);
      for (int j = 0; j < 3; j++) {
        ifs >> target[j];
      }
      input.targets.push_back(target);
    }
  }

  return input;
}

void output_data(int case_num, const Answer& answer) {
  if (exec_mode == 0) {
    // 標準出力
    for (int i = 0; i < answer.initial_board.n; i++) {
      for (int j = 0; j < answer.initial_board.n - 1; j++) {
        cout << answer.initial_board.v[i][j] << " ";
      }
      cout << endl;
    }
    for (int i = 0; i < answer.initial_board.n - 1; i++) {
      for (int j = 0; j < answer.initial_board.n; j++) {
        cout << answer.initial_board.h[i][j] << " ";
      }
      cout << endl;
    }
    for (int i = 0; i < answer.t; i++) {
      if (answer.turns[i][0] == 1) {
        cout << "1 " << answer.turns[i][1] << " " << answer.turns[i][2] << " " << answer.turns[i][3] << endl;
      }
      else if (answer.turns[i][0] == 2) {
        cout << "2 " << answer.turns[i][1] << " " << answer.turns[i][2] << endl;
      }
      else if (answer.turns[i][0] == 3) {
        cout << "3 " << answer.turns[i][1] << " " << answer.turns[i][2] << endl;
      }
      else if (answer.turns[i][0] == 4) {
        cout << "4 " << answer.turns[i][1] << " " << answer.turns[i][2] << " " << answer.turns[i][3] << " " << answer.turns[i][4] << endl;
      }
    }
  }
  else {
    // ファイル出力
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofstream ofs(oss.str());

    for (int i = 0; i < answer.initial_board.n; i++) {
      for (int j = 0; j < answer.initial_board.n - 1; j++) {
        ofs << answer.initial_board.v[i][j] << " ";
      }
      ofs << endl;
    }
    for (int i = 0; i < answer.initial_board.n - 1; i++) {
      for (int j = 0; j < answer.initial_board.n; j++) {
        ofs << answer.initial_board.h[i][j] << " ";
      }
      ofs << endl;
    }
    for (int i = 0; i < answer.t; i++) {
      if (answer.turns[i][0] == 1) {
        ofs << "1 " << answer.turns[i][1] << " " << answer.turns[i][2] << " " << answer.turns[i][3] << endl;
      }
      else if (answer.turns[i][0] == 2) {
        ofs << "2 " << answer.turns[i][1] << " " << answer.turns[i][2] << endl;
      }
      else if (answer.turns[i][0] == 3) {
        ofs << "3 " << answer.turns[i][1] << " " << answer.turns[i][2] << endl;
      }
      else if (answer.turns[i][0] == 4) {
        ofs << "4 " << answer.turns[i][1] << " " << answer.turns[i][2] << " " << answer.turns[i][3] << " " << answer.turns[i][4] << endl;
      }
    }

    if (ofs.is_open()) {
      ofs.close();
    }
  }
}

ll calculate_score() {
  ll res = 0;
  return res;
}

ll solve_case(int case_num) {
  start_timer();

  Input input = input_data(case_num);

  Answer answer(input.n, input.t);

  for (int i = 0; i < input.h; i++) {
    answer.add_turn_1(0, 0, 0);
    answer.add_turn_2(0, 0);
  }

  output_data(case_num, answer);

  ll score = 0;
  if (exec_mode != 0) {
    score = calculate_score();
  }
  return score;
}

int main() {
  exec_mode = 2;

  if (exec_mode == 0) {
    solve_case(0);
  }
  else if (exec_mode <= 2) {
    ll sum_score = 0;
    for (int i = 0; i < 15; i++) {
      ll score = solve_case(i);
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
