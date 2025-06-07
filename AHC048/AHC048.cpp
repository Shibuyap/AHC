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

constexpr int INF = 1001001001;
constexpr double EPS = 1e-6;

const int DX[4] = { -1, 0, 1, 0 };
const int DY[4] = { 0, -1, 0, 1 };
constexpr int UP = 0;    // 上
constexpr int LEFT = 1;  // 左
constexpr int DOWN = 2;  // 下
constexpr int RIGHT = 3; // 右

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

int exec_mode;

// 誤差を計算
inline double calc_error(const vector<double>& col1, const vector<double>& col2) {
  return sqrt(
    (col1[0] - col2[0]) * (col1[0] - col2[0])
    + (col1[1] - col2[1]) * (col1[1] - col2[1])
    + (col1[2] - col2[2]) * (col1[2] - col2[2])
  );
}

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
  vector<vector<int>> v; // 垂直の壁
  vector<vector<int>> h; // 水平の壁

  vector<vector<int>> counts; // セルが含まれるウェルのサイズ
  vector<vector<double>> volumes; // セルが含まれるウェルの絵具の量
  vector<vector<vector<double>>> colors; // セルが含まれるウェルの色(CMY)

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

    calc_mixed_color_vec.resize(3, 0.0); // RGB
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

  vector<double> calc_mixed_color_vec;
  inline vector<double> calc_mixed_color(const vector<double>& col1, const vector<double>& col2, double vol1, double vol2) {
    double total_volume = vol1 + vol2;
    for (int i = 0; i < 3; i++) {
      if (total_volume < EPS) {
        calc_mixed_color_vec[i] = 0.0;
      }
      else {
        calc_mixed_color_vec[i] = (vol1 * col1[i] + vol2 * col2[i]) / total_volume;
      }
    }
    return calc_mixed_color_vec;
  }

  void add_turn_1(int x, int y, const vector<double>& col) {
    // 実際に加えることのできる量
    double w = min(1.0, counts[x][y] - volumes[x][y]);
    double after_vol = volumes[x][y] + w;

    // 混ぜた後の絵の具の色
    vector<double> mixed_color = calc_mixed_color(colors[x][y], col, volumes[x][y], w);

    // セルの絵具量を更新
    auto cells = get_well_cells(x, y);
    for (const auto& cell : cells) {
      int cx = cell.first;
      int cy = cell.second;
      volumes[cx][cy] = after_vol;
      colors[cx][cy] = mixed_color;
    }
  }

  void add_turn_2(int x, int y) {
    auto cells = get_well_cells(x, y);
    for (const auto& cell : cells) {
      int cx = cell.first;
      int cy = cell.second;
      volumes[cx][cy] = max(0.0, volumes[cx][cy] - 1.0);
    }
  }

  void add_turn_3(int x, int y) {
    add_turn_2(x, y);
  }

  void add_turn_4(int x, int y, int d) {
    int nx = x + DX[d];
    int ny = y + DY[d];

    int before_wall = get_wall(x, y, d);
    if (before_wall == 0) {
      // 仕切りを上げる
      int before_count = counts[x][y];
      double before_volume = volumes[x][y];
      toggle_wall(x, y, d);
      auto after_cells = get_well_cells(x, y);
      if (after_cells.size() == before_count) {
        // 変化なし
      }
      else {
        auto other_cells = get_well_cells(nx, ny);
        for (const auto& cell : after_cells) {
          int cx = cell.first;
          int cy = cell.second;
          volumes[cx][cy] = before_volume * after_cells.size() / before_count;
        }
        for (const auto& cell : other_cells) {
          int cx = cell.first;
          int cy = cell.second;
          volumes[cx][cy] = before_volume * other_cells.size() / before_count;
        }
      }
    }
    else {
      // 仕切りを下す
      int before_count = counts[x][y];
      int before_other_count = counts[nx][ny];
      toggle_wall(x, y, d);
      auto after_cells = get_well_cells(x, y);
      if (after_cells.size() == before_count) {
        // 変化なし
      }
      else {
        // 混ぜた後の絵の具の色
        auto mixed_color = calc_mixed_color(colors[x][y], colors[nx][ny], volumes[x][y], volumes[nx][ny]);
        double after_volume = volumes[x][y] + volumes[nx][ny];
        for (const auto& cell : after_cells) {
          int cx = cell.first;
          int cy = cell.second;
          counts[cx][cy] = after_cells.size();
          volumes[cx][cy] = after_volume;
          colors[cx][cy] = mixed_color;
        }
      }
    }
  }
};

class Answer {
public:
  Board initial_board;
  Board board;

  int max_t;
  int t;
  bool is_over;
  vector<vector<int>> turns;

  Answer(int _n, int _max_t) : max_t(_max_t), initial_board(_n), board(_n) {
    t = 0;
    is_over = false;
    for (int i = 0; i < max_t; i++) {
      turns.push_back(vector<int>(5, 0));
    }
  }

  void clear() {
    t = 0;
    board = initial_board;
  }

  void sim_turn_1(int x, int y, int k, const Input& input) {
    board.add_turn_1(x, y, input.owns[k]);
  }

  // 絵具をウェルに追加する
  void add_turn_1(int x, int y, int k, const Input& input) {
    if (t >= max_t) {
      if (exec_mode != 778) {
        cerr << "Error: add_turn_1 called after max_t reached." << endl;
      }
      is_over = true;
      return;
    }

    turns[t][0] = 1;
    turns[t][1] = x;
    turns[t][2] = y;
    turns[t][3] = k;
    turns[t][4] = 0; // dummy value
    t++;

    sim_turn_1(x, y, k, input);
  }

  bool can_turn_2(int x, int y) {
    return board.volumes[x][y] >= 1.0 - 1e-6;
  }

  void sim_turn_2(int x, int y) {
    board.add_turn_2(x, y);
  }

  // 絵具を画伯に渡す
  void add_turn_2(int x, int y) {
    if (!can_turn_2(x, y)) {
      if (exec_mode != 778) {
        cerr << "Error: add_turn_2 called when not enough paint is available." << endl;
      }
      is_over = true;
      return;
    }

    if (t >= max_t) {
      if (exec_mode != 778) {
        cerr << "Error: add_turn_2 called after max_t reached." << endl;
      }
      is_over = true;
      return;
    }
    turns[t][0] = 2;
    turns[t][1] = x;
    turns[t][2] = y;
    turns[t][3] = 0; // dummy value
    turns[t][4] = 0; // dummy value
    t++;

    sim_turn_2(x, y);
  }

  void sim_turn_3(int x, int y) {
    board.add_turn_3(x, y);
  }

  // 絵具を破棄する
  void add_turn_3(int x, int y) {
    if (t >= max_t) {
      cerr << "Error: add_turn_3 called after max_t reached." << endl;
      is_over = true;
      return;
    }
    turns[t][0] = 3;
    turns[t][1] = x;
    turns[t][2] = y;
    turns[t][3] = 0; // dummy value
    turns[t][4] = 0; // dummy value
    t++;

    sim_turn_3(x, y);
  }

  void sim_turn_4(int x, int y, int d) {
    board.add_turn_4(x, y, d);
  }

  // 仕切りを出し入れする
  void add_turn_4(int x, int y, int d) {
    if (t >= max_t) {
      if (exec_mode != 778) {
        cerr << "Error: add_turn_4 called after max_t reached." << endl;
      }
      is_over = true;
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

    sim_turn_4(x, y, d);
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

struct Score {
  int score;
  int d_score;
  int e_scrore;
  double max_e_score;

  vector<double> e_score_list; // 各ターンの誤差を格納するリスト

  Score() : score(0), d_score(0), e_scrore(0), max_e_score(0.0) {}
};

Score calculate_score(Answer& ans, const Input& input) {
  Score score;
  score.e_score_list.resize(input.h, 0.0);

  ans.board = ans.initial_board;
  int count_add_1 = 0;
  int count_add_2 = 0;
  double sum_e_score = 0.0;
  for (int i = 0; i < ans.t; i++) {
    if (ans.turns[i][0] == 1) { // add_turn_1
      int x = ans.turns[i][1];
      int y = ans.turns[i][2];
      int k = ans.turns[i][3];
      ans.sim_turn_1(x, y, k, input);
      count_add_1++;
    }
    else if (ans.turns[i][0] == 2) { // add_turn_2
      int x = ans.turns[i][1];
      int y = ans.turns[i][2];
      if (!ans.can_turn_2(x, y)) {
        cerr << i << " Error: add_turn_2 called when not enough paint is available." << endl;
      }
      ans.sim_turn_2(x, y);
      count_add_2++;

      // 誤差を計算
      double e_score = calc_error(ans.board.colors[x][y], input.targets[count_add_2 - 1]);
      sum_e_score += e_score;
      score.e_score_list[count_add_2 - 1] = e_score;
      if (e_score > score.max_e_score) {
        score.max_e_score = e_score;
      }
    }
    else if (ans.turns[i][0] == 3) { // add_turn_3
      int x = ans.turns[i][1];
      int y = ans.turns[i][2];
      ans.sim_turn_3(x, y);
    }
    else if (ans.turns[i][0] == 4) { // add_turn_4
      int x = ans.turns[i][1];
      int y = ans.turns[i][2];
      int nx = ans.turns[i][3];
      int ny = ans.turns[i][4];
      for (int d = 0; d < 4; d++) {
        if (x + DX[d] == nx && y + DY[d] == ny) {
          ans.sim_turn_4(x, y, d);
          break;
        }
      }
    }
  }

  if (count_add_2 != input.h) {
  }

  score.d_score = input.d * (count_add_1 - input.h);
  score.e_scrore = round(sum_e_score * 1e4);
  score.score = 1 + score.d_score + score.e_scrore;

  sort(score.e_score_list.begin(), score.e_score_list.end());

  return score;
}

void initialize_board_axb(Answer& answer, int a, int b) {
  for (int i = 0; i < answer.initial_board.n; i++) {
    for (int j = 0; j < answer.initial_board.n - 1; j++) {
      if (j % a == a - 1) {
        answer.initial_board.v[i][j] = 1;
      }
      else {
        answer.initial_board.v[i][j] = 0;
      }
    }
  }
  for (int i = 0; i < answer.initial_board.n - 1; i++) {
    for (int j = 0; j < answer.initial_board.n; j++) {
      if (i % b == b - 1) {
        answer.initial_board.h[i][j] = 1;
      }
      else {
        answer.initial_board.h[i][j] = 0;
      }
    }
  }
  answer.initial_board.calc_counts();

  answer.board = answer.initial_board;
}

void initialize_board_solver_56(Answer& answer) {
  for (int i = 0; i < answer.initial_board.n; i++) {
    for (int j = 0; j < answer.initial_board.n - 1; j++) {
      if (j == 0) {
        answer.initial_board.v[i][j] = 1;
      }
      else {
        answer.initial_board.v[i][j] = 0;
      }
    }
  }
  for (int i = 0; i < answer.initial_board.n - 1; i++) {
    for (int j = 0; j < answer.initial_board.n; j++) {
      if (j == 0) {
        answer.initial_board.h[i][j] = 0;
      }
      else {
        answer.initial_board.h[i][j] = 1;
      }
    }
  }
  answer.initial_board.calc_counts();

  answer.board = answer.initial_board;
}

class Solver_1 {
public:
  Answer answer;
  const Input& input;
  Score score;

  Solver_1(Answer ans, const Input& in) : answer(ans), input(in) {
    score.score = INF;
  }

  void solve() {
    initialize_board_axb(answer, 1, 1);

    answer.clear();

    for (int i = 0; i < input.h; i++) {
      answer.add_turn_1(0, 0, rand_xorshift() % input.k, input);
      answer.add_turn_1(0, 0, rand_xorshift() % input.k, input);
      answer.add_turn_2(0, 0);
      answer.add_turn_3(0, 0);
    }

    score = calculate_score(answer, input);
  }
};

ll nCr(ll n, ll r) {
  if (r < 0 || r > n) return 0;
  r = min(r, n - r);
  ll res = 1;
  for (ll i = 1; i <= r; ++i) {
    res = res * (n - r + i) / i;
  }
  return res;
}

vector<vector<int>> create_candidates(int max_num, int max_count) {
  vector<vector<int>> res;
  vector<int> cur;

  // 深さ優先で非減少列を生成
  function<void(int, int)> dfs = [&](int start, int remain)
    {
      if (remain == 0) {            // ちょうど max_count 個そろった
        res.push_back(cur);
        return;
      }
      for (int i = start; i < max_num; ++i) {
        cur.push_back(i);         // i を選ぶ
        dfs(i, remain - 1);       // 同じ i を重複して選べるので start は i
        cur.pop_back();           // 戻す
      }
    };

  if (max_count >= 0 && max_num > 0) dfs(0, max_count);
  return res;
}

vector<pair<vector<double>, vector<int>>> create_candidate_pairs(int max_num, int max_count, const Input& input) {
  vector<pair<vector<double>, vector<int>>> res;
  vector<int> cur;
  // 深さ優先で非減少列を生成
  function<void(int, int)> dfs = [&](int start, int remain)
    {
      if (remain == 0) {            // ちょうど max_count 個そろった
        vector<double> col(3, 0.0);
        for (int i = 0; i < cur.size(); i++) {
          for (int j = 0; j < 3; j++) {
            col[j] += input.owns[cur[i]][j];
          }
        }
        for (int j = 0; j < 3; j++) {
          col[j] /= cur.size();
        }
        res.push_back(make_pair(col, cur));
        return;
      }
      for (int i = start; i < max_num; ++i) {
        cur.push_back(i);         // i を選ぶ
        dfs(i, remain - 1);       // 同じ i を重複して選べるので start は i
        cur.pop_back();           // 戻す
      }
    };
  if (max_count >= 0 && max_num > 0) dfs(0, max_count);

  sort(res.begin(), res.end());
  return res;
}

inline double calc_one_cost(int vol, const vector<double>& col1, const vector<double>& col2, int d) {
  return vol * d + calc_error(col1, col2) * 1e4;
}

class Solver_2 {
public:
  double time_limit;
  const Input& input;
  Answer answer;
  Score score;

  Answer best_answer;
  Score best_score;

  Solver_2(Answer ans, const Input& in, double tl) : input(in), answer(ans), best_answer(ans) {
    score.score = INF;
    best_score.score = INF;
    time_limit = tl;
  }

  void solve() {
    initialize_board_axb(answer, 20, 20);

    vector<pair<vector<double>, vector<int>>> candidates;

    for (int paint_count = 1; paint_count < 100; paint_count++) {
      //cerr << paint_count << endl;

      // 重複組み合わせが大きい時break
      if (nCr(input.k + paint_count - 1, paint_count) > 2e6) {
        break;
      }

      auto this_candidates = create_candidate_pairs(input.k, paint_count, input);
      candidates.insert(candidates.end(), this_candidates.begin(), this_candidates.end());

      sort(candidates.begin(), candidates.end());

      answer.clear();
      bool ok = true;
      int max_count = 0;

      for (int i = 0; i < input.h; i++) {
        if (get_elapsed_time() > time_limit) {
          ok = false;
          break;
        }

        vector<double> tar = input.targets[i];
        double min_cost = INF;
        int min_iter = -1;

        int start_iter = lower_bound(candidates.begin(), candidates.end(), make_pair(tar, vector<int>())) - candidates.begin();
        if (start_iter == candidates.size()) {
          start_iter--;
        }

        vector<double> best_cand_col = tar;
        for (int j = start_iter; j >= 0; j--) {
          double cost = calc_one_cost(candidates[j].second.size() - 1, tar, candidates[j].first, input.d);
          if (cost < min_cost) {
            min_cost = cost;
            min_iter = j;
          }

          best_cand_col[0] = candidates[j].first[0];
          if (calc_one_cost(0, tar, best_cand_col, input.d) > min_cost) {
            break;
          }
        }

        for (int j = start_iter; j < candidates.size(); j++) {
          double cost = calc_one_cost(candidates[j].second.size() - 1, tar, candidates[j].first, input.d);
          if (cost < min_cost) {
            min_cost = cost;
            min_iter = j;
          }

          best_cand_col[0] = candidates[j].first[0];
          if (calc_one_cost(0, tar, best_cand_col, input.d) > min_cost) {
            break;
          }
        }

        int need_turn_count = candidates[min_iter].second.size() * 2;
        if (answer.t + need_turn_count > answer.max_t) {
          ok = false;
          break;
        }

        max_count = max(max_count, (int)candidates[min_iter].second.size());

        for (auto num : candidates[min_iter].second) {
          answer.add_turn_1(0, 0, num, input);
        }
        for (int j = 0; j < candidates[min_iter].second.size(); j++) {
          if (j == 0) {
            answer.add_turn_2(0, 0);
          }
          else {
            answer.add_turn_3(0, 0);
          }
        }
      }

      if (ok) {
        score = calculate_score(answer, input);
        if (score.score < best_score.score) {
          best_answer = answer;
          best_score = score;
        }
      }

      if (!ok) {
        break;
      }

      if (max_count < paint_count) {
        break;
      }
    }

    answer = best_answer;
    score = best_score;
  }
};

// 仕切り固定(T小向け)
class Solver_3 {
public:
  double time_limit;
  const Input& input;
  Answer answer;
  Score score;

  Answer best_answer;
  Score best_score;
  int best_num = 0;

  Solver_3(Answer ans, const Input& in, double tl) : input(in), answer(ans), best_answer(ans) {
    score.score = INF;
    best_score.score = INF;
    time_limit = tl;
  }

  void solve() {
    vector<pair<vector<double>, vector<int>>> candidates[6];
    for (int num = 1; num <= 5; num++) {
      candidates[num] = create_candidate_pairs(input.k, num, input);
    }

    for (int num = 2; num <= 5; num++) {
      initialize_board_axb(answer, num, 1);

      answer.clear();
      bool ok = true;
      int max_count = 0;

      for (int i = 0; i < input.h; i++) {
        if (get_elapsed_time() > time_limit) {
          ok = false;
          break;
        }

        double min_cost = INF;
        int query = 0;
        int min_x = -1;
        int min_y = -1;
        int col_indices[5] = { -1, -1, -1, -1, -1 };
        bool visited_vacant = false;

        for (int x = 0; x < input.n; x++) {
          for (int y = 0; y + num - 1 < input.n; y += num) {
            if (answer.board.volumes[x][y] < EPS) {
              if (visited_vacant) {
                continue;
              }
              visited_vacant = true;

              for (int col_count = 1; col_count <= num; col_count++) {
                for (int j = 0; j < candidates[col_count].size(); j++) {
                  double cost = calc_one_cost(0, input.targets[i], candidates[col_count][j].first, input.d);
                  if (cost < min_cost) {
                    min_cost = cost;
                    query = 0;
                    min_x = x;
                    min_y = y;
                    for (int k = 0; k < 5; k++) {
                      if (k < col_count) {
                        col_indices[k] = candidates[col_count][j].second[k];
                      }
                      else {
                        col_indices[k] = -1;
                      }
                    }
                  }
                }
              }
            }
            else {
              // 足さない
              double cost = calc_one_cost(0, input.targets[i], answer.board.colors[x][y], input.d);
              if (cost < min_cost) {
                min_cost = cost;
                query = 2;
                min_x = x;
                min_y = y;
                for (int k = 0; k < 5; k++) {
                  col_indices[k] = -1;
                }
              }

              int vol_int = (int)(answer.board.volumes[x][y] + 0.5);
              for (int col_count = 1; col_count <= num - vol_int; col_count++) {
                for (int j = 0; j < candidates[col_count].size(); j++) {
                  auto mixed_color = answer.board.calc_mixed_color(answer.board.colors[x][y], candidates[col_count][j].first, answer.board.volumes[x][y], col_count);
                  double cost = calc_one_cost(0, input.targets[i], mixed_color, input.d);
                  if (cost < min_cost) {
                    min_cost = cost;
                    query = 1;
                    min_x = x;
                    min_y = y;
                    for (int k = 0; k < 5; k++) {
                      if (k < col_count) {
                        col_indices[k] = candidates[col_count][j].second[k];
                      }
                      else {
                        col_indices[k] = -1;
                      }
                    }
                  }
                }
              }
            }
          }
        }

        if (query == 0) {
          // 空のウェル
          for (int j = 0; j < num; j++) {
            if (col_indices[j] == -1) {
              continue;
            }
            answer.add_turn_1(min_x, min_y, col_indices[j], input);
          }
        }
        else if (query == 1) {
          // 空でないウェルに追加
          for (int j = 0; j < num; j++) {
            if (col_indices[j] == -1) {
              continue;
            }
            answer.add_turn_1(min_x, min_y, col_indices[j], input);
          }
        }
        else if (query == 2) {
          // 空でないウェルから取得
          // 何もしない
        }
        answer.add_turn_2(min_x, min_y);
      }

      if (ok) {
        score = calculate_score(answer, input);
        if (score.score < best_score.score) {
          best_answer = answer;
          best_score = score;
          best_num = num;
        }
      }

      if (!ok) {
        break;
      }
    }

    answer = best_answer;
    score = best_score;
  }
};

class Solver_4 {
public:
  double time_limit;
  const Input& input;
  Answer answer;
  Score score;

  Answer best_answer;
  Score best_score;
  int best_num = 0;

  static const int MAX_NUM = 5;

  Solver_4(Answer ans, const Input& in, double tl) : input(in), answer(ans), best_answer(ans) {
    score.score = INF;
    best_score.score = INF;
    time_limit = tl;
  }

  void solve() {
    vector<pair<vector<double>, vector<int>>> candidates[MAX_NUM + 1];
    for (int num = 1; num <= MAX_NUM; num++) {
      candidates[num] = create_candidate_pairs(input.k, num, input);
    }

    for (int num = 2; num <= MAX_NUM; num++) {
      initialize_board_axb(answer, num, 1);

      answer.clear();
      bool ok = true;
      int max_count = 0;

      for (int i = 0; i < input.h; i++) {
        if (get_elapsed_time() > time_limit) {
          ok = false;
          break;
        }

        double min_cost = INF;
        int query = 0;
        int min_x = -1;
        int min_y = -1;
        int col_indices[5] = { -1, -1, -1, -1, -1 };
        int eliminate_count = 0;
        bool visited_vacant = false;

        for (int x = 0; x < input.n; x++) {
          for (int y = 0; y + num - 1 < input.n; y += num) {
            if (answer.board.volumes[x][y] < EPS) {
              if (visited_vacant) {
                continue;
              }
              visited_vacant = true;

              for (int col_count = 1; col_count <= num; col_count++) {
                for (int j = 0; j < candidates[col_count].size(); j++) {
                  double cost = calc_one_cost(0, input.targets[i], candidates[col_count][j].first, input.d);
                  if (cost < min_cost) {
                    min_cost = cost;
                    query = 0;
                    min_x = x;
                    min_y = y;
                    for (int k = 0; k < 5; k++) {
                      if (k < col_count) {
                        col_indices[k] = candidates[col_count][j].second[k];
                      }
                      else {
                        col_indices[k] = -1;
                      }
                    }
                  }
                }
              }
            }
            else {
              // 足さない
              double cost = calc_one_cost(0, input.targets[i], answer.board.colors[x][y], input.d);
              if (cost < min_cost) {
                min_cost = cost;
                query = 2;
                min_x = x;
                min_y = y;
                for (int k = 0; k < 5; k++) {
                  col_indices[k] = -1;
                }
              }

              int vol_int = (int)(answer.board.volumes[x][y] + 0.5);
              for (int elim_count = 0; elim_count <= vol_int - 1; elim_count++) {
                for (int col_count = 1; col_count <= num - vol_int; col_count++) {
                  for (int j = 0; j < candidates[col_count].size(); j++) {
                    auto mixed_color = answer.board.calc_mixed_color(answer.board.colors[x][y], candidates[col_count][j].first, answer.board.volumes[x][y] - elim_count, col_count);
                    double cost = calc_one_cost(elim_count, input.targets[i], mixed_color, input.d);
                    if (cost < min_cost) {
                      min_cost = cost;
                      query = 1;
                      min_x = x;
                      min_y = y;
                      for (int k = 0; k < 5; k++) {
                        if (k < col_count) {
                          col_indices[k] = candidates[col_count][j].second[k];
                        }
                        else {
                          col_indices[k] = -1;
                        }
                      }
                      eliminate_count = elim_count;
                    }
                  }
                }
              }
            }
          }
        }

        if (!visited_vacant) {
          int vol_int = (int)(answer.board.volumes[0][0] + 0.5);
          for (int col_count = 1; col_count <= num; col_count++) {
            for (int j = 0; j < candidates[col_count].size(); j++) {
              double cost = calc_one_cost(vol_int, input.targets[i], candidates[col_count][j].first, input.d);
              if (cost < min_cost) {
                min_cost = cost;
                query = 3;
                min_x = 0;
                min_y = 0;
                for (int k = 0; k < 5; k++) {
                  if (k < col_count) {
                    col_indices[k] = candidates[col_count][j].second[k];
                  }
                  else {
                    col_indices[k] = -1;
                  }
                }
              }
            }
          }
        }

        if (query == 0) {
          // 空のウェル
          for (int j = 0; j < num; j++) {
            if (col_indices[j] == -1) {
              continue;
            }
            answer.add_turn_1(min_x, min_y, col_indices[j], input);
          }
        }
        else if (query == 1) {
          // 空でないウェルに追加
          for (int j = 0; j < eliminate_count; j++) {
            answer.add_turn_3(min_x, min_y);
          }
          for (int j = 0; j < num; j++) {
            if (col_indices[j] == -1) {
              continue;
            }
            answer.add_turn_1(min_x, min_y, col_indices[j], input);
          }
        }
        else if (query == 2) {
          // 空でないウェルから取得
          // 何もしない
        }
        else if (query == 3) {
          // 空にしてから追加
          while (answer.board.volumes[min_x][min_y] > EPS) {
            answer.add_turn_3(min_x, min_y);
            if (answer.is_over) {
              ok = false;
              break;
            }
          }
          for (int j = 0; j < num; j++) {
            if (col_indices[j] == -1) {
              continue;
            }
            answer.add_turn_1(min_x, min_y, col_indices[j], input);
          }
        }
        answer.add_turn_2(min_x, min_y);

        if (answer.is_over) {
          ok = false;
          break;
        }
      }

      if (ok) {
        score = calculate_score(answer, input);
        if (score.score < best_score.score) {
          best_answer = answer;
          best_score = score;
          best_num = num;
        }
      }

      if (!ok) {
        break;
      }
    }

    answer = best_answer;
    score = best_score;
  }
};

class Solver_5 {
public:
  class Solver5State {
  public:
    vector<int> vertical_lines;
    int change_vertical_lines_count = 0;
    double mixed_volume;
    vector<double> mixed_colors;
    double score;

    Solver5State(int n) : vertical_lines(n, 0), change_vertical_lines_count(0), mixed_volume(0.0), mixed_colors(3, 0.0), score(1e12) {}

    inline double calc_eval(const vector<double>& target_color, int change_vertical_lines_count_limit, int input_d) {
      if (mixed_volume < 1.0 || 3.0 < mixed_volume || change_vertical_lines_count > change_vertical_lines_count_limit) {
        score = 1e9 + abs(1.1 - mixed_volume) + change_vertical_lines_count_limit * 100;
      }
      else {
        score = calc_error(mixed_colors, target_color) * 1e4 + (mixed_volume - 1.0) * min(input_d, 400) * 0.0001;
      }
      return score;
    }

    inline double calc_eval_2(const vector<double>& target_color, int change_vertical_lines_count_limit, int input_d) {
      if (mixed_volume < 1.0 || 1.2 < mixed_volume) {
        score = 1e9 + calc_error(mixed_colors, target_color) * 1e4 + abs(mixed_volume - 1.0) * min(input_d, 400) * 100000;
      }
      else {
        score = calc_error(mixed_colors, target_color) * 1e4;
      }
      return score;
    }
  };

  class Solver5Params {
  public:
    double minimum_volume = 20 + EPS;
    int reset_vertical_line_num = 5;
    double discard_value = 0.4;
    int start_change_vertical_lines_count_limit = 999;

    int initial_set_count = 40;
    int initial_loop_count = 20;
    int initial_operation_thresholds[10] = { 0, 100, 300, 400, 500, 600, 700, 800, 900, 1000 };
    int initial_max_threshold = 100;

    int main_loop_count = 500000; // 1000
    int main_operation_thresholds[10] = { 90, 100, 150, 400, 500, 600, 700, 800, 900, 1000 };
    int main_max_threshold = 150;
    int main_break_last_update_iter = 999999; // 50
  };

  double time_limit;
  const Input& input;
  Answer answer;
  Score score;

  Answer best_answer;
  Score best_score;
  int best_num = 0;

  Solver_5(Answer ans, const Input& in, double tl) : input(in), answer(ans), best_answer(ans) {
    score.score = INF;
    best_score.score = INF;
    time_limit = tl;
  }

  int change_vertical_lines_count_limit = 999;

  inline void move_diff(const Solver5State& state, Solver5State& next_state, int idx1, int diff1, double diff_vol) {
    next_state.mixed_colors = answer.board.calc_mixed_color(next_state.mixed_colors, answer.board.colors[idx1][input.n - 1], next_state.mixed_volume, diff_vol);
    next_state.mixed_volume += diff_vol;
    if (next_state.vertical_lines[idx1] == state.vertical_lines[idx1]) {
      next_state.change_vertical_lines_count++;
    }
    next_state.vertical_lines[idx1] += diff1;
    if (next_state.vertical_lines[idx1] == state.vertical_lines[idx1]) {
      next_state.change_vertical_lines_count--;
    }
  }

  void solve_inner(int i, const Solver5Params& params, const Solver5State& state, Solver5State& next_state, int& count1, int& count2, int mode, vector<int>& used_indices, const vector<vector<int>>& each_colors, const vector<int>& each_row_colors) {
    int thresholds[10];
    int max_threshold = 0;
    int loop_count = 0;
    if (mode == 0) {
      for (int j = 0; j < 10; j++) {
        thresholds[j] = params.initial_operation_thresholds[j];
      }
      max_threshold = params.initial_max_threshold;
      loop_count = params.initial_loop_count;
    }
    else {
      for (int j = 0; j < 10; j++) {
        thresholds[j] = params.main_operation_thresholds[j];
      }
      max_threshold = params.main_max_threshold;
      loop_count = params.main_loop_count;
    }

    used_indices = calc_used_indices(state, next_state);
    next_state.calc_eval(input.targets[i], change_vertical_lines_count_limit, input.d);

    int last_update_iter = 0;
    for (int iter = 0; iter < loop_count; iter++) {
      double current_score = next_state.score;

      int idx1 = 0;
      int idx2 = 0;
      int diff1 = 0;
      int col1 = 0;
      double diff_vol = 0.0;
      double diff_vol2 = 0.0;

      int ra = rand_xorshift() % max_threshold;
      if (ra < thresholds[0]) {
        // 仕切りを+-1する
        idx1 = rand_xorshift() % input.n;
        diff1 = rand_xorshift() % 2 == 0 ? 1 : -1;
        if (next_state.vertical_lines[idx1] + diff1 < state.vertical_lines[idx1] || next_state.vertical_lines[idx1] + diff1 > input.n - 2) {
          diff1 *= -1;
        }

        diff_vol = answer.board.volumes[idx1][input.n - 1] / answer.board.counts[idx1][input.n - 1] * diff1;
        move_diff(state, next_state, idx1, diff1, diff_vol);
      }
      else if (ra < thresholds[1]) {
        // 仕切りをランダムに変える
        idx1 = rand_xorshift() % input.n;
        // [vertical_lines[idx1], input.n - 2]
        int new1 = rand_xorshift() % (input.n - 2 - state.vertical_lines[idx1] + 1) + state.vertical_lines[idx1];
        if (rand_xorshift() % 3 == 0) {
          new1 = state.vertical_lines[idx1];
        }
        while (new1 == next_state.vertical_lines[idx1]) {
          new1 = rand_xorshift() % (input.n - 2 - state.vertical_lines[idx1] + 1) + state.vertical_lines[idx1];
        }
        diff1 = new1 - next_state.vertical_lines[idx1];

        diff_vol = answer.board.volumes[idx1][input.n - 1] / answer.board.counts[idx1][input.n - 1] * diff1;
        move_diff(state, next_state, idx1, diff1, diff_vol);
      }
      else if (ra < thresholds[2]) {
        // 同じ色の仕切りを一つ-1して、別の仕切りを+1する
        diff1 = -1;
        if (used_indices.size() == 0) {
          //iter--;
          continue;
        }
        for (int j = 0; j < 10; j++) {
          idx1 = used_indices[rand_xorshift() % used_indices.size()];
          col1 = each_row_colors[idx1];
          if (each_colors[col1].size() > 1) {
            break;
          }
        }
        if (each_colors[col1].size() <= 1) {
          //iter--;
          continue;
        }
        idx2 = each_colors[col1][rand_xorshift() % each_colors[col1].size()];
        while (idx1 == idx2) {
          idx2 = each_colors[col1][rand_xorshift() % each_colors[col1].size()];
        }
        if (next_state.vertical_lines[idx2] >= input.n - 2) {
          //iter--;
          continue;
        }

        if (next_state.vertical_lines[idx1] <= state.vertical_lines[idx1]) {
          cerr << "Error: next_state.vertical_lines[idx1] <= state.vertical_lines[idx1]: " << next_state.vertical_lines[idx1] << " <= " << state.vertical_lines[idx1] << endl;
        }

        diff_vol = answer.board.volumes[idx1][input.n - 1] / answer.board.counts[idx1][input.n - 1] * diff1;
        move_diff(state, next_state, idx1, diff1, diff_vol);
        diff_vol2 = answer.board.volumes[idx2][input.n - 1] / answer.board.counts[idx2][input.n - 1] * -diff1;
        move_diff(state, next_state, idx2, -diff1, diff_vol2);
      }

      next_state.calc_eval(input.targets[i], change_vertical_lines_count_limit, input.d);
      count1++;
      double diff_score = (current_score - next_state.score) * 12345.6;
      const double START_TEMP = 1000000; // 初期温度
      const double END_TEMP = 0.01; // 終了温度
      double temp = START_TEMP - (START_TEMP - END_TEMP) * (double)iter / loop_count;
      bool accept = false;
      if (mode == 0) {
        accept = next_state.score < current_score;
      }
      else {
        if (temp <= 0.0) {
          accept = next_state.score < current_score;
        }
        else {
          double prob = exp(diff_score / temp);
          accept = prob > rand_01();
        }
      }
      if (accept) {
        last_update_iter = iter;
        // スコアが改善した
        count2++;
        used_indices = calc_used_indices(state, next_state);
        //cerr << "Iter: " << iter << ", Score: " << next_state.score << ", Change: " << diff_score << ", Temp: " << temp;
      }
      else {
        // 戻す
        next_state.score = current_score;
        if (ra < thresholds[0]) {
          move_diff(state, next_state, idx1, -diff1, -diff_vol);
        }
        else if (ra < thresholds[1]) {
          move_diff(state, next_state, idx1, -diff1, -diff_vol);
        }
        else if (ra < thresholds[2]) {
          move_diff(state, next_state, idx2, diff1, -diff_vol2);
          move_diff(state, next_state, idx1, -diff1, -diff_vol);
        }
      }
      if (iter > last_update_iter + params.main_break_last_update_iter) {
        break;
      }
    }
  }

  vector<int> calc_used_indices(const Solver5State& state, const Solver5State& next_state) {
    vector<int> used_indices;
    for (int j = 0; j < input.n; j++) {
      if (next_state.vertical_lines[j] != state.vertical_lines[j]) {
        used_indices.push_back(j);
      }
    }
    return used_indices;
  }

  void solve(const Solver5Params& params) {

    initialize_board_solver_56(answer);

    answer.clear();
    bool ok = true;
    int max_count = 0;

    vector<int> each_row_colors(input.n, 0);
    for (int i = 0; i < input.n; i++) {
      answer.add_turn_1(i, input.n - 1, i % input.k, input);
      each_row_colors[i] = i % input.k;
    }

    Solver5State state(input.n);

    int count1 = 0, count2 = 0;

    for (int i = 0; i < input.h; i++) {
      //if (i % 100 == 0) {
      //  cerr << "Processing target " << i << " / " << input.h << endl;
      //}
      //cerr << i << endl;
      if (get_elapsed_time() > time_limit) {
        ok = false;
        break;
      }

      for (int j = 0; j < input.n; j++) {
        if (state.vertical_lines[j] >= params.reset_vertical_line_num) {
          double diff_vol = answer.board.volumes[0][0] / answer.board.counts[0][0] * state.vertical_lines[j];
          state.mixed_volume -= diff_vol;
          answer.add_turn_4(j, 0, RIGHT);
          answer.add_turn_4(j, state.vertical_lines[j], RIGHT);
          state.vertical_lines[j] = 0;
        }
      }

      vector<double> each_color_volumes(input.k, 0.0);
      while (true) {
        int idx = rand_xorshift() % input.n;
        double min_vol = 1e9;
        for (int j = 0; j < input.n; j++) {
          double vol = answer.board.volumes[j][input.n - 1];
          if (vol < min_vol) {
            min_vol = vol;
            idx = j;
          }
        }
        if (min_vol >= params.discard_value) {
          break;
        }
        for (int j = 0; j < input.k; j++) {
          each_color_volumes[j] = 0.0;
        }
        for (int j = 0; j < input.n; j++) {
          each_color_volumes[each_row_colors[j]] += answer.board.volumes[j][input.n - 1];
        }
        int col = 0;
        double min_color_vol = 1e9;
        for (int j = 0; j < input.k; j++) {
          if (each_color_volumes[j] < min_color_vol) {
            min_color_vol = each_color_volumes[j];
            col = j;
          }
        }
        answer.add_turn_3(idx, input.n - 1);
        answer.add_turn_1(idx, input.n - 1, col, input);
        each_row_colors[idx] = col;
        if (answer.is_over) {
          ok = false;
          break;
        }
      }
      if (answer.is_over) {
        ok = false;
        break;
      }

      vector<vector<int>> each_colors(input.k);
      for (int j = 0; j < input.n; j++) {
        each_colors[each_row_colors[j]].push_back(j);
      }

      Solver5State next_state = state;

      // 山登り
      change_vertical_lines_count_limit = params.start_change_vertical_lines_count_limit;
      int attempt_count = 0;
      int attempt_count_2 = 0;

      Solver5State best_best_state = state;
      Solver5State best_state = state;
      auto used_indices = calc_used_indices(state, next_state);
      while (true) {
        attempt_count++;

        best_state = state;

        for (int set_iter = 0; set_iter < params.initial_set_count; set_iter++) {
          next_state = state;
          used_indices = calc_used_indices(state, next_state);

          if (false && attempt_count == 1 && set_iter == 0) {
            while (next_state.mixed_volume < 1.0) {
              int best_idx = -1;
              double best_vol = 0.0;
              double best_score = 1e12;
              for (int j = 0; j < input.n; j++) {
                if (next_state.vertical_lines[j] >= params.reset_vertical_line_num) {
                  continue;
                }
                if (next_state.vertical_lines[j] + 1 > input.n - 2) {
                  continue;
                }
                double diff_vol = answer.board.volumes[j][input.n - 1] / answer.board.counts[j][input.n - 1];
                auto mixed_color = answer.board.calc_mixed_color(next_state.mixed_colors, answer.board.colors[j][input.n - 1], next_state.mixed_volume, diff_vol);
                double new_score = calc_error(mixed_color, input.targets[i]) * 1e4;
                if (new_score < best_score) {
                  best_score = new_score;
                  best_idx = j;
                  best_vol = diff_vol;
                }
              }
              move_diff(state, next_state, best_idx, 1, best_vol);
            }
            used_indices = calc_used_indices(state, next_state);
          }
          else {
            solve_inner(i, params, state, next_state, count1, count2, 0, used_indices, each_colors, each_row_colors);
          }

          next_state.calc_eval(input.targets[i], change_vertical_lines_count_limit, input.d);
          if (next_state.score < best_state.score) {
            // bestに保存
            best_state = next_state;
          }
        }

        // bestから戻す
        next_state = best_state;
        used_indices = calc_used_indices(state, next_state);

        solve_inner(i, params, state, next_state, count1, count2, 1, used_indices, each_colors, each_row_colors);
        cerr << "i = " << i << ", attempt_count = " << attempt_count << ", error = " << calc_error(next_state.mixed_colors, input.targets[i]) * 1e4 << ", mixed_volume = " << next_state.mixed_volume << endl;

        if (next_state.mixed_volume < 1.0 - EPS) {
          change_vertical_lines_count_limit++;
          cerr << change_vertical_lines_count_limit << " " << next_state.change_vertical_lines_count << endl;
          cerr << "next_state.mixed_volume = " << next_state.mixed_volume << endl;
          for (int i = 0; i < input.n; i++) {
            cerr << next_state.vertical_lines[i] << " ";
          }
          cerr << endl;
          continue;
        }

        next_state.calc_eval(input.targets[i], change_vertical_lines_count_limit, input.d);
        //cerr << "i = " << i << ", attempt_count = " << attempt_count << ", next_mixed_volume = " << next_state.mixed_volume << ", new_score = " << next_state.score << ", best_best_score = " << best_best_state.score << endl;
        if (next_state.score < best_best_state.score) {
          best_best_state = next_state;
        }

        if (best_best_state.score > (input.d / 100.0 + 2.0) * (attempt_count + 25) / 26) {
          if (attempt_count == 3 && attempt_count_2 == 0) {
            while (state.mixed_volume > EPS) {
              answer.add_turn_3(0, 0);
              state.mixed_volume = max(0.0, state.mixed_volume - 1.0);
            }

            attempt_count = 0;
            attempt_count_2++;

            best_best_state = state;
          }

          next_state = state;
          used_indices = calc_used_indices(state, next_state);

          if (attempt_count_2 > 0 && attempt_count > 20) {
            change_vertical_lines_count_limit++;
          }

          if (attempt_count_2 == 0 || attempt_count <= 30) {
            continue;
          }
        }

        // best_bestから戻す
        next_state = best_best_state;
        used_indices = calc_used_indices(state, next_state);

        break;
      }

      //cerr << i << " : attempt_count = " << attempt_count << endl;
      //cerr << attempt_count << ' ';
      //if(i == input.h - 1) {
      //  cerr << endl;
      //}

      for (int j = 0; j < input.n; j++) {
        if (next_state.vertical_lines[j] != state.vertical_lines[j]) {
          answer.add_turn_4(j, next_state.vertical_lines[j], RIGHT);
          answer.add_turn_4(j, state.vertical_lines[j], RIGHT);
        }
      }

      answer.add_turn_2(0, 0);
      next_state.mixed_volume = max(0.0, next_state.mixed_volume - 1.0);

      state = next_state;
      state.score = 1e12;
      state.change_vertical_lines_count = 0;

      //cout << "i = " << i << ", mixed_volume = " << mixed_volume << ", " << answer.board.volumes[0][0] << endl;

      if (answer.is_over) {
        ok = false;
        break;
      }
    }

    if (ok) {
      //cerr << "count1 = " << count1 << ", count2 = " << count2 << endl;
      score = calculate_score(answer, input);
    }
    else {
      score.score = INF;
    }
  }
};

class Solver_6 {
public:
  class Solver6Params {
  public:
    double minimum_volume = 20 + EPS;
    int reset_vertical_line_num = 5;
    double discard_value = 0.4;
    int start_change_vertical_lines_count_limit = 999;
    bool is_mixed_volume_reset = true;
    bool is_each_row_volume_reset = true;

    int initial_set_count = 40;
    int initial_loop_count = 20;
    int initial_operation_thresholds[10] = { 0, 100, 300, 400, 500, 600, 700, 800, 900, 1000 };
    int initial_max_threshold = 100;

    int main_loop_count = 500000; // 1000
    int main_operation_thresholds[10] = { 90, 100, 150, 400, 500, 600, 700, 800, 900, 1000 };
    int main_max_threshold = 150;
    int main_break_last_update_iter = 999999; // 50
  };

  class Solver6State {
  public:
    vector<int> vertical_lines;
    int change_vertical_lines_count = 0;
    double mixed_volume;
    vector<double> mixed_colors;
    double score;

    Solver6State(int n) : vertical_lines(n, 0), change_vertical_lines_count(0), mixed_volume(0.0), mixed_colors(3, 0.0), score(1e12) {}

    inline double calc_eval(const vector<double>& target_color, int change_vertical_lines_count_limit, int input_d) {
      if (mixed_volume < 1.0 || 3.0 < mixed_volume || change_vertical_lines_count > change_vertical_lines_count_limit) {
        score = 1e9 + abs(1.1 - mixed_volume) + change_vertical_lines_count_limit * 100;
      }
      else {
        score = calc_error(mixed_colors, target_color) * 1e4 + (mixed_volume - 1.0) * min(input_d, 400) * 0.0001;
      }
      return score;
    }

    inline double calc_eval_2(const vector<double>& target_color, int change_vertical_lines_count_limit, int input_d) {
      if (mixed_volume < 1.0 || 1.2 < mixed_volume) {
        score = 1e9 + calc_error(mixed_colors, target_color) * 1e4 + abs(mixed_volume - 1.0) * min(input_d, 400) * 100000;
      }
      else {
        score = calc_error(mixed_colors, target_color) * 1e4;
      }
      return score;
    }
  };

  class Solver6SAState {
  public:
    int n;
    vector<int> used_indices;
    vector<vector<int>> each_colors;
    vector<int> each_row_colors;

    vector<vector<pair<double, pair<int, int>>>> each_row_each_operation_volumes;
    vector<int> each_row_each_operation_indices;

    void calc_used_indices() {
      used_indices.clear();
      for (int i = 0; i < n; i++) {
        if (each_row_each_operation_indices[i] != 0) {
          used_indices.push_back(i);
        }
      }
    }
  };

  double time_limit;
  const Input& input;
  Answer answer;
  Score score;

  Answer best_answer;
  Score best_score;
  int best_num = 0;

  Solver_6(Answer ans, const Input& in, double tl) : input(in), answer(ans), best_answer(ans) {
    score.score = INF;
    best_score.score = INF;
    time_limit = tl;
  }

  int change_vertical_lines_count_limit = 999;

  inline void move_diff(const Solver6State& state, Solver6State& next_state, int idx1, int diff1, double diff_vol) {
    next_state.mixed_colors = answer.board.calc_mixed_color(next_state.mixed_colors, answer.board.colors[idx1][input.n - 1], next_state.mixed_volume, diff_vol);
    next_state.mixed_volume += diff_vol;
    if (next_state.vertical_lines[idx1] == state.vertical_lines[idx1]) {
      next_state.change_vertical_lines_count++;
    }
    next_state.vertical_lines[idx1] += diff1;
    if (next_state.vertical_lines[idx1] == state.vertical_lines[idx1]) {
      next_state.change_vertical_lines_count--;
    }
  }

  void solve_inner(int i, const Solver6Params& params, const Solver6State& state, Solver6State& next_state, int mode, Solver6SAState& saState) {
    int thresholds[10];
    int max_threshold = 0;
    int loop_count = 0;
    if (mode == 0) {
      for (int j = 0; j < 10; j++) {
        thresholds[j] = params.initial_operation_thresholds[j];
      }
      max_threshold = params.initial_max_threshold;
      loop_count = params.initial_loop_count;
    }
    else {
      for (int j = 0; j < 10; j++) {
        thresholds[j] = params.main_operation_thresholds[j];
      }
      max_threshold = params.main_max_threshold;
      loop_count = params.main_loop_count;
    }

    saState.calc_used_indices();
    next_state.calc_eval(input.targets[i], change_vertical_lines_count_limit, input.d);

    int last_update_iter = 0;
    for (int iter = 0; iter < loop_count; iter++) {
      double current_score = next_state.score;

      int idx1 = 0;
      int idx2 = 0;
      int diff1 = 0;
      int col1 = 0;
      double diff_vol = 0.0;
      double diff_vol2 = 0.0;

      int ra = rand_xorshift() % max_threshold;
      if (ra < thresholds[0]) {
        // 仕切りを+-1する
        idx1 = rand_xorshift() % input.n;
        diff1 = rand_xorshift() % 2 == 0 ? 1 : -1;
        if (next_state.vertical_lines[idx1] + diff1 < state.vertical_lines[idx1] || next_state.vertical_lines[idx1] + diff1 > input.n - 2) {
          diff1 *= -1;
        }

        diff_vol = answer.board.volumes[idx1][input.n - 1] / answer.board.counts[idx1][input.n - 1] * diff1;
        move_diff(state, next_state, idx1, diff1, diff_vol);
      }
      else if (ra < thresholds[1]) {
        // 仕切りをランダムに変える
        idx1 = rand_xorshift() % input.n;
        // [vertical_lines[idx1], input.n - 2]
        int new1 = rand_xorshift() % (input.n - 2 - state.vertical_lines[idx1] + 1) + state.vertical_lines[idx1];
        if (rand_xorshift() % 3 == 0) {
          new1 = state.vertical_lines[idx1];
        }
        while (new1 == next_state.vertical_lines[idx1]) {
          new1 = rand_xorshift() % (input.n - 2 - state.vertical_lines[idx1] + 1) + state.vertical_lines[idx1];
        }
        diff1 = new1 - next_state.vertical_lines[idx1];

        diff_vol = answer.board.volumes[idx1][input.n - 1] / answer.board.counts[idx1][input.n - 1] * diff1;
        move_diff(state, next_state, idx1, diff1, diff_vol);
      }
      else if (ra < thresholds[2]) {
        // 同じ色の仕切りを一つ-1して、別の仕切りを+1する
        diff1 = -1;
        if (saState.used_indices.size() == 0) {
          continue;
        }
        for (int j = 0; j < 10; j++) {
          idx1 = saState.used_indices[rand_xorshift() % saState.used_indices.size()];
          col1 = saState.each_row_colors[idx1];
          if (saState.each_colors[col1].size() > 1) {
            break;
          }
        }
        if (saState.each_colors[col1].size() <= 1) {
          continue;
        }
        idx2 = saState.each_colors[col1][rand_xorshift() % saState.each_colors[col1].size()];
        while (idx1 == idx2) {
          idx2 = saState.each_colors[col1][rand_xorshift() % saState.each_colors[col1].size()];
        }
        if (next_state.vertical_lines[idx2] >= input.n - 2) {
          continue;
        }

        if (next_state.vertical_lines[idx1] <= state.vertical_lines[idx1]) {
          cerr << "Error: next_state.vertical_lines[idx1] <= state.vertical_lines[idx1]: " << next_state.vertical_lines[idx1] << " <= " << state.vertical_lines[idx1] << endl;
        }

        diff_vol = answer.board.volumes[idx1][input.n - 1] / answer.board.counts[idx1][input.n - 1] * diff1;
        move_diff(state, next_state, idx1, diff1, diff_vol);
        diff_vol2 = answer.board.volumes[idx2][input.n - 1] / answer.board.counts[idx2][input.n - 1] * -diff1;
        move_diff(state, next_state, idx2, -diff1, diff_vol2);
      }

      next_state.calc_eval(input.targets[i], change_vertical_lines_count_limit, input.d);
      double diff_score = (current_score - next_state.score) * 12345.6;
      const double START_TEMP = 1000000; // 初期温度
      const double END_TEMP = 0.01; // 終了温度
      double temp = START_TEMP - (START_TEMP - END_TEMP) * (double)iter / loop_count;
      bool accept = false;
      if (mode == 0) {
        accept = next_state.score < current_score;
      }
      else {
        if (temp <= 0.0) {
          accept = next_state.score < current_score;
        }
        else {
          double prob = exp(diff_score / temp);
          accept = prob > rand_01();
        }
      }
      if (accept) {
        // スコアが改善した
        last_update_iter = iter;
        saState.calc_used_indices();
      }
      else {
        // 戻す
        next_state.score = current_score;
        if (ra < thresholds[0]) {
          move_diff(state, next_state, idx1, -diff1, -diff_vol);
        }
        else if (ra < thresholds[1]) {
          move_diff(state, next_state, idx1, -diff1, -diff_vol);
        }
        else if (ra < thresholds[2]) {
          move_diff(state, next_state, idx2, diff1, -diff_vol2);
          move_diff(state, next_state, idx1, -diff1, -diff_vol);
        }
      }
      if (iter > last_update_iter + params.main_break_last_update_iter) {
        break;
      }
    }
  }

  void solve(const Solver6Params& params) {
    initialize_board_solver_56(answer);

    answer.clear();
    bool ok = true;

    Solver6SAState saState;
    saState.n = input.n;
    saState.each_row_colors.resize(input.n);
    for (int i = 0; i < input.n; i++) {
      answer.add_turn_1(i, input.n - 1, i % input.k, input);
      saState.each_row_colors[i] = i % input.k;
    }
    saState.each_row_each_operation_volumes.resize(input.n);
    saState.each_row_each_operation_indices.resize(input.n);

    Solver6State state(input.n);

    for (int i = 0; i < input.h; i++) {
      for (int j = 0; j < input.n; j++) {
        if (state.vertical_lines[j] >= params.reset_vertical_line_num) {
          double diff_vol = answer.board.volumes[0][0] / answer.board.counts[0][0] * state.vertical_lines[j];
          state.mixed_volume -= diff_vol;
          answer.add_turn_4(j, 0, RIGHT);
          answer.add_turn_4(j, state.vertical_lines[j], RIGHT);
          state.vertical_lines[j] = 0;
        }
      }

      // each_row_each_operation_volumesを計算
      // answer.board.volumes[0][0]は0.0の前提
      for (int j = 0; j < input.n; j++) {
        saState.each_row_each_operation_indices[j] = 0;
        saState.each_row_each_operation_volumes[j].clear();
        saState.each_row_each_operation_volumes[j].emplace_back(0.0, make_pair(state.vertical_lines[j], state.vertical_lines[j]));

        double current_row_volume = answer.board.volumes[j][input.n - 1] / answer.board.counts[j][input.n - 1];

        // 広げる → 縮める
        for (int k = state.vertical_lines[j]; k >= 0; k--) {
          double one_block_volume = current_row_volume / (input.n - 1 - k);
          for (int l = k + 1; l <= input.n - 3; l++) {
            double diff_vol = one_block_volume * (l - k);
            saState.each_row_each_operation_volumes[j].emplace_back(diff_vol, make_pair(k, l));
          }
        }

        sort(saState.each_row_each_operation_volumes[j].begin(), saState.each_row_each_operation_volumes[j].end());
      }

      vector<double> each_color_volumes(input.k, 0.0);
      while (true) {
        int idx = rand_xorshift() % input.n;
        double min_vol = 1e9;
        for (int j = 0; j < input.n; j++) {
          double vol = answer.board.volumes[j][input.n - 1];
          if (vol < min_vol) {
            min_vol = vol;
            idx = j;
          }
        }
        if (min_vol >= params.discard_value) {
          break;
        }
        for (int j = 0; j < input.k; j++) {
          each_color_volumes[j] = 0.0;
        }
        for (int j = 0; j < input.n; j++) {
          each_color_volumes[saState.each_row_colors[j]] += answer.board.volumes[j][input.n - 1];
        }
        int col = 0;
        double min_color_vol = 1e9;
        for (int j = 0; j < input.k; j++) {
          if (each_color_volumes[j] < min_color_vol) {
            min_color_vol = each_color_volumes[j];
            col = j;
          }
        }
        if (params.is_each_row_volume_reset) {
          answer.add_turn_3(idx, input.n - 1);
        }
        answer.add_turn_1(idx, input.n - 1, col, input);
        saState.each_row_colors[idx] = col;
        if (answer.is_over) {
          ok = false;
          break;
        }
      }
      if (answer.is_over) {
        ok = false;
        break;
      }

      saState.each_colors.clear();
      saState.each_colors.resize(input.k);
      for (int j = 0; j < input.n; j++) {
        saState.each_colors[saState.each_row_colors[j]].push_back(j);
      }

      Solver6State next_state = state;

      // 山登り
      change_vertical_lines_count_limit = params.start_change_vertical_lines_count_limit;
      int attempt_count = 0;
      int attempt_count_2 = 0;

      Solver6State best_best_state = state;
      Solver6State best_state = state;
      saState.calc_used_indices();
      while (true) {
        attempt_count++;

        best_state = state;

        for (int set_iter = 0; set_iter < params.initial_set_count; set_iter++) {
          next_state = state;
          saState.calc_used_indices();

          if (false && attempt_count == 1 && set_iter == 0) {
            while (next_state.mixed_volume < 1.0) {
              int best_idx = -1;
              double best_vol = 0.0;
              double best_score = 1e12;
              for (int j = 0; j < input.n; j++) {
                if (next_state.vertical_lines[j] >= params.reset_vertical_line_num) {
                  continue;
                }
                if (next_state.vertical_lines[j] + 1 > input.n - 2) {
                  continue;
                }
                double diff_vol = answer.board.volumes[j][input.n - 1] / answer.board.counts[j][input.n - 1];
                auto mixed_color = answer.board.calc_mixed_color(next_state.mixed_colors, answer.board.colors[j][input.n - 1], next_state.mixed_volume, diff_vol);
                double new_score = calc_error(mixed_color, input.targets[i]) * 1e4;
                if (new_score < best_score) {
                  best_score = new_score;
                  best_idx = j;
                  best_vol = diff_vol;
                }
              }
              move_diff(state, next_state, best_idx, 1, best_vol);
            }
            saState.calc_used_indices();
          }
          else {
            solve_inner(i, params, state, next_state, 0, saState);
          }

          next_state.calc_eval(input.targets[i], change_vertical_lines_count_limit, input.d);
          if (next_state.score < best_state.score) {
            // bestに保存
            best_state = next_state;
          }
        }

        // bestから戻す
        next_state = best_state;
        saState.calc_used_indices();

        solve_inner(i, params, state, next_state, 1, saState);

        if (next_state.mixed_volume < 1.0 - EPS) {
          cerr << "next_state.mixed_volume = " << next_state.mixed_volume << endl;
          continue;
        }

        next_state.calc_eval(input.targets[i], change_vertical_lines_count_limit, input.d);
        if (next_state.score < best_best_state.score) {
          best_best_state = next_state;
        }

        if (best_best_state.score > (input.d / 100.0 + 2.0) * (attempt_count + 25) / 26) {
          if (attempt_count == 3 && attempt_count_2 == 0) {
            if (params.is_mixed_volume_reset) {
              while (state.mixed_volume > EPS) {
                answer.add_turn_3(0, 0);
                state.mixed_volume = max(0.0, state.mixed_volume - 1.0);
              }
            }

            attempt_count = 0;
            attempt_count_2++;

            best_best_state = state;
          }

          next_state = state;
          saState.calc_used_indices();

          if (attempt_count_2 > 0 && attempt_count > 20) {
            change_vertical_lines_count_limit++;
          }

          if (attempt_count_2 == 0 || attempt_count <= 30) {
            continue;
          }
        }

        // best_bestから戻す
        next_state = best_best_state;
        saState.calc_used_indices();

        break;
      }

      cerr << "i = " << i << ", attempt_count = " << attempt_count << ", error = " << calc_error(next_state.mixed_colors, input.targets[i]) * 1e4 << ", mixed_volume = " << next_state.mixed_volume << endl;

      for (int j = 0; j < input.n; j++) {
        if (next_state.vertical_lines[j] != state.vertical_lines[j]) {
          answer.add_turn_4(j, next_state.vertical_lines[j], RIGHT);
          answer.add_turn_4(j, state.vertical_lines[j], RIGHT);
        }
      }

      answer.add_turn_2(0, 0);
      next_state.mixed_volume = max(0.0, next_state.mixed_volume - 1.0);

      state = next_state;
      state.score = 1e12;
      state.change_vertical_lines_count = 0;

      if (answer.is_over) {
        ok = false;
        break;
      }
    }

    if (ok) {
      score = calculate_score(answer, input);
    }
    else {
      score.score = INF;
    }
  }
};

ll solve_case(int case_num) {
  const double TIME_LIMIT = 2.8;
  start_timer();

  Input input = input_data(case_num);

  //if (input.t >= 8000) {
  //  return 0;
  //}

  Answer answer(input.n, input.t);
  Score best_score;
  best_score.score = INF * 2;
  int best_solver = -1;

  //{
  //  auto solver = Solver_1(answer, input);
  //  solver.solve();
  //  if (solver.score.score < best_score.score) {
  //    answer = solver.answer;
  //    best_score = solver.score;
  //    best_solver = 1;
  //  }
  //}

  //{
  //  auto solver = Solver_3(answer, input, 99.9);
  //  solver.solve();
  //  if (solver.score.score < best_score.score) {
  //    answer = solver.answer;
  //    best_score = solver.score;
  //    best_solver = 3 * 10 + solver.best_num;
  //  }
  //}

  //{
  //  auto solver = Solver_4(answer, input, 99.9);
  //  solver.solve();
  //  if (solver.score.score < best_score.score) {
  //    answer = solver.answer;
  //    best_score = solver.score;
  //    best_solver = 4 * 10 + solver.best_num;
  //  }
  //}

  input.d = 1;

  //if (true || (input.t > 20000 && input.d > 2000)) {
  //  int arr[2] = { 20,4 };
  //  for (int i = 0; i < 1; i++) {
  //    auto solver = Solver_5(answer, input, 99999.9);

  //    Solver_5::Solver5Params params;
  //    params.start_change_vertical_lines_count_limit = arr[i];

  //    solver.solve(params);
  //    if (solver.score.score < best_score.score) {
  //      answer = solver.answer;
  //      best_score = solver.score;
  //      best_solver = 5 * 100 + arr[i];
  //    }
  //    if (!solver.answer.is_over) {
  //      break;
  //    }
  //  }
  //}

  if (true || (input.t > 20000 && input.d > 2000)) {
    int arr[2] = { 20,4 };
    for (int i = 0; i < 1; i++) {
      auto solver = Solver_6(answer, input, 99999.9);

      Solver_6::Solver6Params params;
      params.start_change_vertical_lines_count_limit = arr[i];

      solver.solve(params);
      if (solver.score.score < best_score.score) {
        answer = solver.answer;
        best_score = solver.score;
        best_solver = 6 * 100 + arr[i];
      }
      if (!solver.answer.is_over) {
        break;
      }
    }
  }

  //{
  //  auto solver = Solver_2(answer, input, 19.9);
  //  solver.solve();
  //  if (solver.score.score < best_score.score) {
  //    answer = solver.answer;
  //    best_score = solver.score;
  //    best_solver = 2;
  //  }
  //}

  //if (best_solver / 100 != 5) {
  //  return 0;
  //}

  output_data(case_num, answer);

  int score = 0;
  if (exec_mode != 0) {
    score = calculate_score(answer, input).score;
  }

  if (exec_mode == 1) {
    cerr << score << endl;
  }
  else if (exec_mode == 777 || exec_mode == 778) {
    if (best_score.score < 1e7) {
      //cerr << "Turn = " << answer.t << endl;
      double e_ruisekiwa = 0;
      for (int i = 0; i < input.h; i++) {
        e_ruisekiwa += best_score.e_score_list[input.h - 1 - i] * 1e4;
        if (i % 100 == 99) {
          cerr << setw(6) << (int)e_ruisekiwa << ", ";
        }
      }
      //cerr << endl;
      cerr
        << setw(3) << case_num << ", "
        << setw(2) << input.k << ", "
        << setw(5) << input.t << ", "
        << setw(4) << input.d << ", "
        << setw(7) << best_score.score << ", "
        << setw(7) << best_score.d_score << ", "
        << setw(7) << best_score.e_scrore << ", "
        << setw(7) << fixed << setprecision(7) << best_score.max_e_score << ", "
        << setw(3) << best_solver << ", "
        << setw(5) << get_elapsed_time() << ", "
        << setw(5) << answer.t << ", "
        << endl;
    }
  }
  else {
    cerr
      << "case = " << setw(3) << case_num << ", "
      << "k = " << setw(2) << input.k << ", "
      << "t = " << setw(5) << input.t << ", "
      << "d = " << setw(4) << input.d << ", "
      << "score = " << setw(7) << best_score.score << ", "
      << "d_score = " << setw(7) << best_score.d_score << ", "
      << "e_score = " << setw(7) << best_score.e_scrore << ", "
      << "max_e_score = " << setw(7) << fixed << setprecision(7) << best_score.max_e_score << ", "
      << "solver = " << setw(3) << best_solver << ", "
      << "time = " << setw(5) << get_elapsed_time() << " ms"
      << endl;
  }

  return score;
}

int main() {
  exec_mode = 777;

  if (exec_mode == 0) {
    solve_case(0);
  }
  else if (exec_mode <= 999) {
    ll sum_score = 0;
    for (int i = 1; i < 2; i++) {
      ll score = solve_case(i * 5);
      sum_score += score;
    }
  }

  return 0;
}
