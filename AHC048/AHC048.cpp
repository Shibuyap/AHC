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

const int DX[4] = { -1, 0, 1, 0 };
const int DY[4] = { 0, -1, 0, 1 };

// 2次元キューのクラス
class Queue2D {
private:
  static const int MAX_SIZE = 10000;
  int arr[MAX_SIZE][2];
  int head;
  int tail;

public:
  // コンストラクタ
  Queue2D() : head(0), tail(0) {}

  void clear_queue() {
    head = 0;
    tail = 0;
  }

  int front_x() const {
    return arr[head][0];
  }

  int front_y() const {
    return arr[head][1];
  }

  void push(int x, int y) {
    arr[tail][0] = x;
    arr[tail][1] = y;
    tail++;
  }

  void pop() {
    head++;
  }

  int size() const {
    return tail - head;
  }
};
Queue2D queue2d;

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

  vector<vector<int>> counts;
  vector<vector<double>> volumes;
  vector<vector<vector<double>>> colors;

  bool is_ng(int x, int y, int dir) {
    int nx = x + DX[dir];
    int ny = y + DY[dir];
    if (nx < 0 || nx >= n || ny < 0 || ny >= n) {
      return true;
    }
    if (dir == 0) { // 上
      return h[nx][ny] == 1;
    }
    else if (dir == 1) { // 左
      return v[nx][ny] == 1;
    }
    else if (dir == 2) { // 下
      return h[x][y] == 1;
    }
    else if (dir == 3) { // 右
      return v[x][y] == 1;
    }
  }

  void calc_counts() {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        counts[i][j] = 0;
      }
    }

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n - 1; j++) {
        if (counts[i][j] != 0) {
          continue;
        }
        vector<pair<int, int>> vec;
        queue2d.clear_queue();
        queue2d.push(i, j);
        vec.push_back(make_pair(i, j));
        counts[i][j] = 1;
        while (queue2d.size()) {
          int x = queue2d.front_x();
          int y = queue2d.front_y();
          queue2d.pop();
          for (int d = 0; d < 4; d++) {
            if (is_ng(x, y, d)) {
              continue;
            }
            int nx = x + DX[d];
            int ny = y + DY[d];
            if (counts[nx][ny] != 0) {
              continue;
            }
            counts[nx][ny] = 1;
            queue2d.push(nx, ny);
            vec.push_back(make_pair(nx, ny));
          }
        }
        for (const auto& p : vec) {
          counts[p.first][p.second] = vec.size();
        }
      }
    }
  }

  Board(int _n) : n(_n) {
    v.resize(n);
    for (int i = 0; i < n; i++) {
      v[i].resize(n - 1, 0);
    }
    h.resize(n - 1);
    for (int i = 0; i < n - 1; i++) {
      h[i].resize(n, 0);
    }

    counts.resize(n);
    for (int i = 0; i < n; i++) {
      counts[i].resize(n, 0);
    }

    volumes.resize(n);
    for (int i = 0; i < n; i++) {
      volumes[i].resize(n, 0.0);
    }

    colors.resize(n);
    for (int i = 0; i < n; i++) {
      colors[i].resize(n);
      for (int j = 0; j < n; j++) {
        colors[i][j].resize(3, 0.0); // RGB
      }
    }
  }

  int get_wall(int x, int y, int d) {
    if (d == 0) { // 上
      return h[x - 1][y];
    }
    else if (d == 1) { // 左
      return v[x][y - 1];
    }
    else if (d == 2) { // 下
      return h[x][y];
    }
    else if (d == 3) { // 右
      return v[x][y];
    }
  }

  void toggle_wall(int x, int y, int d) {
    if (d == 0) { // 上
      h[x - 1][y] = 1 - h[x - 1][y];
    }
    else if (d == 1) { // 左
      v[x][y - 1] = 1 - v[x][y - 1];
    }
    else if (d == 2) { // 下
      h[x][y] = 1 - h[x][y];
    }
    else if (d == 3) { // 右
      v[x][y] = 1 - v[x][y];
    }
  }

  vector<pair<int, int>> get_well_cells(int x, int y) {
    vector<pair<int, int>> vec;
    queue2d.clear_queue();
    queue2d.push(x, y);
    vec.push_back(make_pair(x, y));
    counts[x][y] = 0;
    while (queue2d.size()) {
      int cx = queue2d.front_x();
      int cy = queue2d.front_y();
      queue2d.pop();
      for (int d = 0; d < 4; d++) {
        if (is_ng(cx, cy, d)) {
          continue;
        }
        int nx = cx + DX[d];
        int ny = cy + DY[d];
        if (counts[nx][ny] == 0) {
          continue;
        }
        counts[nx][ny] = 0;
        queue2d.push(nx, ny);
        vec.push_back(make_pair(nx, ny));
      }
    }
    for (const auto& p : vec) {
      counts[p.first][p.second] = vec.size();
    }
    return vec;
  }

  void calc_one_well_count(int x, int y) {
    get_well_cells(x, y);
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

  // 絵具をウェルに追加する
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

    // 絵具をウェルに追加
    auto cells = board.get_well_cells(x, y);
    for (const auto& cell : cells) {
      int cx = cell.first;
      int cy = cell.second;
      board.volumes[cx][cy] += 1.0;
    }
  }

  bool can_turn_2(int x, int y) {
    return board.volumes[x][y] >= 1.0 - 1e-6;
  }

  // 絵具を画伯に渡す
  void add_turn_2(int x, int y) {
    if (!can_turn_2(x, y)) {
      cerr << "Error: add_turn_2 called when not enough paint is available." << endl;
      return;
    }

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

    // 画伯に絵具を渡す
    auto cells = board.get_well_cells(x, y);
    for (const auto& cell : cells) {
      int cx = cell.first;
      int cy = cell.second;
      board.volumes[cx][cy] = max(0.0, board.volumes[cx][cy] - 1.0);
    }
  }

  // 絵具を破棄する
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

    auto cells = board.get_well_cells(x, y);
    for (const auto& cell : cells) {
      int cx = cell.first;
      int cy = cell.second;
      board.volumes[cx][cy] = max(0.0, board.volumes[cx][cy] - 1.0);
    }
  }

  // 仕切りを出し入れする
  void add_turn_4(int x, int y, int d) {
    if (t >= max_t) {
      cerr << "Error: add_turn_4 called after max_t reached." << endl;
      return;
    }

    int nx = x + DX[d];
    int ny = y + DY[d];

    turns[t][0] = 4;
    turns[t][1] = x;
    turns[t][2] = y;
    turns[t][3] = nx;
    turns[t][4] = ny;
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

void initialize_board_1x1(Answer& answer) {
  for (int i = 0; i < answer.initial_board.n; i++) {
    for (int j = 0; j < answer.initial_board.n - 1; j++) {
      answer.initial_board.v[i][j] = 1;
    }
  }
  for (int i = 0; i < answer.initial_board.n - 1; i++) {
    for (int j = 0; j < answer.initial_board.n; j++) {
      answer.initial_board.h[i][j] = 1;
    }
  }

  answer.board = answer.initial_board;
}

void initialize_board_4x1(Answer& answer) {
  for (int i = 0; i < answer.initial_board.n; i++) {
    for (int j = 0; j < answer.initial_board.n - 1; j++) {
      if (j % 4 == 3) {
        answer.initial_board.v[i][j] = 1;
      }
      else {
        answer.initial_board.v[i][j] = 0;
      }
    }
  }
  for (int i = 0; i < answer.initial_board.n - 1; i++) {
    for (int j = 0; j < answer.initial_board.n; j++) {
      answer.initial_board.h[i][j] = 1;
    }
  }

  answer.board = answer.initial_board;
}

void initialize_board_20x20(Answer& answer) {
  for (int i = 0; i < answer.initial_board.n; i++) {
    for (int j = 0; j < answer.initial_board.n - 1; j++) {
      answer.initial_board.v[i][j] = 0;
    }
  }
  for (int i = 0; i < answer.initial_board.n - 1; i++) {
    for (int j = 0; j < answer.initial_board.n; j++) {
      answer.initial_board.h[i][j] = 0;
    }
  }

  answer.board = answer.initial_board;
}

void method_1(Answer& answer, const Input& input) {
  for (int i = 0; i < input.h; i++) {
    answer.add_turn_1(0, 0, 0);
    answer.add_turn_2(0, 0);
  }
}

ll solve_case(int case_num) {
  start_timer();

  Input input = input_data(case_num);

  Answer answer(input.n, input.t);

  initialize_board_4x1(answer);

  method_1(answer, input);

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
