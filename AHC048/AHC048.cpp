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

  inline vector<double> calc_mixed_color(const vector<double>& col1, const vector<double>& col2, double vol1, double vol2) {
    vector<double> mixed_color(3, 0.0);
    double total_volume = vol1 + vol2;
    for (int i = 0; i < 3; i++) {
      mixed_color[i] = (vol1 * col1[i] + vol2 * col2[i]) / total_volume;
    }
    return mixed_color;
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
  vector<vector<int>> turns;

  Answer(int _n, int _max_t) : max_t(_max_t), initial_board(_n), board(_n) {
    t = 0;
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
      cerr << "Error: add_turn_1 called after max_t reached." << endl;
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

    sim_turn_2(x, y);
  }

  void sim_turn_3(int x, int y) {
    board.add_turn_3(x, y);
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

    sim_turn_3(x, y);
  }

  void sim_turn_4(int x, int y, int d) {
    board.add_turn_4(x, y, d);
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

  Score() : score(0), d_score(0), e_scrore(0), max_e_score(0.0) {}
};

Score calculate_score(Answer& ans, const Input& input) {
  Score score;

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

  return score;
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
  answer.initial_board.calc_counts();

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
  answer.initial_board.calc_counts();

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
    initialize_board_1x1(answer);

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
    initialize_board_20x20(answer);

    vector<pair<vector<double>, vector<int>>> candidates;

    for (int paint_count = 1; paint_count < 100; paint_count++) {
      //cerr << paint_count << endl;

      // 重複組み合わせが大きい時break
      if (nCr(input.k + paint_count - 1, paint_count) > 2e5) {
        break;
      }

      auto this_candidates = create_candidates(input.k, paint_count);
      for (const auto& c : this_candidates) {
        vector<double> col(3, 0.0);
        for (int i = 0; i < paint_count; i++) {
          for (int j = 0; j < 3; j++) {
            col[j] += input.owns[c[i]][j];
          }
        }
        for (int j = 0; j < 3; j++) {
          col[j] /= paint_count;
        }
        candidates.push_back(make_pair(col, c));
      }

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

ll solve_case(int case_num) {
  const double TIME_LIMIT = 2.8;
  start_timer();

  Input input = input_data(case_num);

  Answer answer(input.n, input.t);
  Score best_score;
  best_score.score = INF;

  auto solver_1 = Solver_1(answer, input);
  solver_1.solve();
  if (solver_1.score.score < best_score.score) {
    answer = solver_1.answer;
    best_score = solver_1.score;
  }

  auto solver_2 = Solver_2(answer, input, 2.0);
  solver_2.solve();
  if (solver_2.score.score < best_score.score) {
    answer = solver_2.answer;
    best_score = solver_2.score;
  }

  output_data(case_num, answer);

  int score = 0;
  if (exec_mode != 0) {
    score = calculate_score(answer, input).score;
  }

  if (exec_mode == 1) {
    cerr << score << endl;
  }
  else if (exec_mode == 777) {
    cerr
      << setw(3) << case_num << ", "
      << setw(2) << input.k << ", "
      << setw(5) << input.t << ", "
      << setw(4) << input.d << ", "
      << setw(7) << best_score.score << ", "
      << setw(7) << best_score.d_score << ", "
      << setw(7) << best_score.e_scrore << ", "
      << setw(7) << fixed << setprecision(7) << best_score.max_e_score << ", "
      << setw(5) << get_elapsed_time() << ", "
      << endl;
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
      << "max_e_score = " << setw(7) << best_score.max_e_score << ", "
      << "time = " << setw(5) << get_elapsed_time() << ", "
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
    for (int i = 0; i < 1000; i++) {
      ll score = solve_case(i);
      sum_score += score;
    }
  }

  return 0;
}
