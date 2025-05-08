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
typedef pair<int, int> P;
typedef pair<P, P> PP;

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

// 2次元キューのクラス
class Queue2D
{
private:
  static const int MAX_SIZE = 10000;
  int arr[MAX_SIZE][2];
  int head;
  int tail;

public:
  // コンストラクタ
  Queue2D() : head(0), tail(0) {
  }

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

const ll INF = 1001001001001001001LL;
const int INT_INF = 1001001001;

const int DX[4] = { -1, 0, 1, 0 };
const int DY[4] = { 0, -1, 0, 1 };

const double TIME_LIMIT = 1.9;
int exec_mode;

const int n = 20;
const int m = 40;

const char C[3] = { 'M','S','A' };
const char MOVE[4] = { 'U','L','D','R' };

class Answer
{
public:
  int current_score;
  int ans[2000][2];
  int ans_count;
  int add_flags[m][4];
  int add_flags_2[n + 2][n + 2][m][4];

  void initialize_state() {
    current_score = 0;
    ans_count = 0;
    for (int i = 0; i < (m); ++i) {
      for (int j = 0; j < (4); ++j) {
        add_flags[i][j] = 0;
      }
    }
    for (int i = 0; i < (n + 2); ++i) {
      for (int j = 0; j < (n + 2); ++j) {
        for (int k = 0; k < (m); ++k) {
          for (int l = 0; l < (4); ++l) {
            add_flags_2[i][j][k][l] = 0;
          }
        }
      }
    }
  }

  void Copy(const Answer& src) {
    current_score = src.current_score;
    ans_count = src.ans_count;
    for (int i = 0; i < (ans_count); ++i) {
      for (int j = 0; j < (2); ++j) {
        ans[i][j] = src.ans[i][j];
      }
    }
    for (int i = 0; i < (m); ++i) {
      for (int j = 0; j < (4); ++j) {
        add_flags[i][j] = src.add_flags[i][j];
      }
    }
    for (int i = 0; i < (n + 2); ++i) {
      for (int j = 0; j < (n + 2); ++j) {
        for (int k = 0; k < (m); ++k) {
          for (int l = 0; l < (4); ++l) {
            add_flags_2[i][j][k][l] = src.add_flags_2[i][j][k][l];
          }
        }
      }
    }
  }
};

class Board
{
public:
  int X[m], Y[m];
  int board[n + 2][n + 2];
  int rock_row[n + 2][n + 2];
  int rock_row_count[n + 2];
  int rock_col[n + 2][n + 2];
  int rock_col_count[n + 2];

  bool is_out_of_range(int x, int y) const {
    if (board[x][y] < 0) return true;
    return false;
  }

  void init_board() {
    for (int i = 0; i < (n + 2); ++i) {
      for (int j = 0; j < (n + 2); ++j) {
        if (i == 0 || i == n + 1 || j == 0 || j == n + 1) {
          board[i][j] = -1;
        }
        else {
          board[i][j] = 0;
        }
      }
    }
    for (int i = 0; i < (m); ++i) {
      board[X[i]][Y[i]] = i;
    }
    for (int i = 1; i < n + 1; ++i) {
      rock_row_count[i] = 0;
      rock_col_count[i] = 0;
    }
  }

  bool attempt_skate(int dir, int x, int y, int& nx, int& ny) const {
    nx = x;
    ny = y;
    if (dir == 0) {
      int ma = 0;
      for (int i = 0; i < (rock_col_count[y]); ++i) {
        if (ma < rock_col[y][i] && rock_col[y][i] < x) {
          ma = rock_col[y][i];
        }
      }
      nx = ma + 1;
    }
    else if (dir == 1) {
      int ma = 0;
      for (int i = 0; i < (rock_row_count[x]); ++i) {
        if (ma < rock_row[x][i] && rock_row[x][i] < y) {
          ma = rock_row[x][i];
        }
      }
      ny = ma + 1;
    }
    else if (dir == 2) {
      int mi = n + 1;
      for (int i = 0; i < (rock_col_count[y]); ++i) {
        if (x < rock_col[y][i] && rock_col[y][i] < mi) {
          mi = rock_col[y][i];
        }
      }
      nx = mi - 1;
    }
    else if (dir == 3) {
      int mi = n + 1;
      for (int i = 0; i < (rock_row_count[x]); ++i) {
        if (y < rock_row[x][i] && rock_row[x][i] < mi) {
          mi = rock_row[x][i];
        }
      }
      ny = mi - 1;
    }
    if (nx == x && ny == y) {
      return false;
    }
    return true;
  }

  bool attempt_move(int dir, int x, int y, int& nx, int& ny) const {
    nx = x;
    ny = y;
    if (!is_out_of_range(nx + DX[dir], ny + DY[dir])) {
      nx += DX[dir];
      ny += DY[dir];
    }
    if (nx == x && ny == y)return false;
    return true;
  }
};

static void move_one(Board& board, Answer& answer, int dir, int& x, int& y) {
  int nx = x;
  int ny = y;
  if (!board.attempt_move(dir, x, y, nx, ny)) {
    return;
  }
  answer.ans[answer.ans_count][0] = 0;
  answer.ans[answer.ans_count][1] = dir;
  x += DX[dir];
  y += DY[dir];
  answer.ans_count++;
}

static void skate_one(Board& board, Answer& answer, int dir, int& x, int& y) {
  int nx = x;
  int ny = y;
  if (!board.attempt_skate(dir, x, y, nx, ny)) {
    return;
  }
  answer.ans[answer.ans_count][0] = 1;
  answer.ans[answer.ans_count][1] = dir;
  x = nx;
  y = ny;
  answer.ans_count++;
}

static void add_one(Board& board, Answer& answer, int dir, int x, int y, bool isSim = false) {
  int nx = x + DX[dir];
  int ny = y + DY[dir];
  if (board.board[nx][ny] != 0) {
    return;
  }
  board.board[nx][ny] = -2;
  board.rock_row[nx][board.rock_row_count[nx]] = ny;
  board.rock_row_count[nx]++;
  board.rock_col[ny][board.rock_col_count[ny]] = nx;
  board.rock_col_count[ny]++;
  if (!isSim) {
    answer.ans[answer.ans_count][0] = 2;
    answer.ans[answer.ans_count][1] = dir;
    answer.ans_count++;
  }
}

static Board input_data(int case_num) {
  Board board;

  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  int _n, _m;
  if (!ifs.is_open()) {
    // 標準入力
    cin >> _n >> _m;
    for (int i = 0; i < (m); ++i) {
      cin >> board.X[i] >> board.Y[i];
    }
  }
  else {
    // ファイル入力
    ifs >> _n >> _m;
    for (int i = 0; i < (m); ++i) {
      ifs >> board.X[i] >> board.Y[i];
    }
  }
  for (int i = 0; i < (m); ++i) {
    board.X[i]++;
    board.Y[i]++;
  }
  board.init_board();
  return board;
}

static void open_ofs(int case_num, ofstream& ofs) {
  if (exec_mode != 0) {
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofs.open(oss.str());
  }
}

static int calculate_score(const Answer& answer) {
  int res = m + 2 * n * m - answer.ans_count;
  return res;
}

static void output_data(const Answer& answer, ofstream& ofs) {
  if (exec_mode == 0) {
    // 標準出力
    for (int i = 0; i < (answer.ans_count); ++i) {
      cout << C[answer.ans[i][0]] << ' ' << MOVE[answer.ans[i][1]] << endl;
    }
  }
  else {
    // ファイル出力
    for (int i = 0; i < (answer.ans_count); ++i) {
      ofs << C[answer.ans[i][0]] << ' ' << MOVE[answer.ans[i][1]] << endl;
    }
  }
}

Queue2D que2d;
int dp[n + 2][n + 2];
int dp2[n + 2][n + 2][4];
static vector<P> bfs(Board& board, int sx, int sy, int gx, int gy) {
  for (int i = 0; i < (n + 2); ++i) {
    for (int j = 0; j < (n + 2); ++j) {
      dp[i][j] = INT_INF;
    }
  }
  que2d.clear_queue();
  que2d.push(sx, sy);
  dp[sx][sy] = 0;
  while (que2d.size()) {
    int x = que2d.front_x();
    int y = que2d.front_y();
    que2d.pop();
    if (x == gx && y == gy) {
      break;
    }
    int nx, ny;
    for (int i = 3; i >= 0; --i) {
      if (board.attempt_skate(i, x, y, nx, ny)) {
        if (dp[x][y] + 1 < dp[nx][ny]) {
          dp[nx][ny] = dp[x][y] + 1;
          dp2[nx][ny][0] = x;
          dp2[nx][ny][1] = y;
          dp2[nx][ny][2] = 1;
          dp2[nx][ny][3] = i;
          que2d.push(nx, ny);
        }
      }
      if (board.attempt_move(i, x, y, nx, ny)) {
        if (dp[x][y] + 1 < dp[nx][ny]) {
          dp[nx][ny] = dp[x][y] + 1;
          dp2[nx][ny][0] = x;
          dp2[nx][ny][1] = y;
          dp2[nx][ny][2] = 0;
          dp2[nx][ny][3] = i;
          que2d.push(nx, ny);
        }
      }
    }
  }

  vector<P> vp;
  if (dp[gx][gy] == INT_INF) {
    return vp;
  }
  int x = gx;
  int y = gy;
  while (x != sx || y != sy) {
    int nx = dp2[x][y][0];
    int ny = dp2[x][y][1];
    vp.push_back(P(dp2[x][y][2], dp2[x][y][3]));
    x = nx;
    y = ny;
  }
  reverse(vp.begin(), vp.end());
  return vp;
}

static int create_ans(Answer& answer, Board& board) {
  answer.ans_count = 0;
  int x = board.X[0];
  int y = board.Y[0];
  for (int i = 0; i < (m); ++i) {
    auto vp = bfs(board, x, y, board.X[i], board.Y[i]);
    for (const auto& p : vp) {
      if (p.first == 0) {
        move_one(board, answer, p.second, x, y);
      }
      else {
        skate_one(board, answer, p.second, x, y);
      }
      for (int j = 0; j < (4); ++j) {
        if (answer.add_flags_2[x][y][i][j] == 1) {
          add_one(board, answer, j, x, y);
        }
      }
    }
    if (x != board.X[i] || y != board.Y[i]) {
      return -1;
    }
    if (i == m - 1) {
      break;
    }
    for (int j = 0; j < (4); ++j) {
      if (answer.add_flags[i][j] == 1) {
        add_one(board, answer, j, x, y);
      }
    }
  }
  return calculate_score(answer);
}

static void simulate_best(Board& board, Answer& answer, const Answer& best_answer, int& x, int& y, int& mm, int turn, int& lastRockX, int& lastRockY, int& lastRockDir, int& lastRockMM) {
  x = board.X[0];
  y = board.Y[0];
  mm = 0;
  lastRockX = -1;
  lastRockY = -1;
  board.init_board();
  for (int i = 0; i < (turn); ++i) {
    if (mm < m && x == board.X[mm] && y == board.Y[mm]) {
      mm++;
    }
    if (best_answer.ans[i][0] == 0) {
      x += DX[best_answer.ans[i][1]];
      y += DY[best_answer.ans[i][1]];
    }
    if (best_answer.ans[i][0] == 1) {
      int nx, ny;
      board.attempt_skate(best_answer.ans[i][1], x, y, nx, ny);
      x = nx;
      y = ny;
    }
    if (best_answer.ans[i][0] == 2) {
      add_one(board, answer, best_answer.ans[i][1], x, y, true);
      int nx = x + DX[best_answer.ans[i][1]];
      int ny = y + DY[best_answer.ans[i][1]];
      if (board.board[nx][ny] == 2) {
        lastRockX = x;
        lastRockX = y;
        lastRockDir = best_answer.ans[i][1];
        lastRockMM = mm;
      }
    }
    if (mm < m && x == board.X[mm] && y == board.Y[mm]) {
      mm++;
    }
  }
}

static Answer build_initial_solution_2(Board& board) {
  Answer answer;
  answer.initialize_state();
  create_ans(answer, board);
  return answer;
}

struct AnnealingParams
{
  double start_temperature[10];
  double end_temperature;
  double score_scale;
  int operation_thresholds[10];
};

static void run_simulated_annealing(Board& board, Answer& answer, AnnealingParams annealingParams) {
  Answer best_answer;
  best_answer.Copy(answer);

  double now_time = get_elapsed_time();
  const double START_TEMP = annealingParams.start_temperature[0];
  const double END_TEMP = annealingParams.end_temperature;

  int loop = 0;
  while (true) {
    loop++;

    if (loop % 1 == 0) {
      now_time = get_elapsed_time();
      if (now_time > TIME_LIMIT) break;
    }

    double progress_ratio = now_time / TIME_LIMIT;
    double temp = START_TEMP + (END_TEMP - START_TEMP) * progress_ratio;

    // 近傍解作成
    int ra_exec_mode = rand_xorshift() % annealingParams.operation_thresholds[2];
    int ra1, ra2, ra3, ra4, ra5;

    if (ra_exec_mode < annealingParams.operation_thresholds[0]) {
      // 近傍操作1
      ra1 = rand_xorshift() % (m - 1);
      ra2 = rand_xorshift() % 4;
      answer.add_flags[ra1][ra2] = 1 - answer.add_flags[ra1][ra2];
    }
    else if (ra_exec_mode < annealingParams.operation_thresholds[1]) {
      // 近傍操作2
      ra1 = rand_xorshift() % (best_answer.ans_count - 10);
      ra2 = rand_xorshift() % 4;
      int _1, _2, _3, _4;
      simulate_best(board, answer, best_answer, ra3, ra4, ra5, ra1, _1, _2, _3, _4);
      answer.add_flags_2[ra3][ra4][ra5][ra2] = 1 - answer.add_flags_2[ra3][ra4][ra5][ra2];
    }
    else if (ra_exec_mode < annealingParams.operation_thresholds[2]) {
      // 近傍操作3
      ra1 = rand_xorshift() % (best_answer.ans_count - 10);
      int _1, _2, _3;
      simulate_best(board, answer, best_answer, _1, _2, _3, ra1, ra2, ra3, ra4, ra5);
      if (ra2 == -1)continue;
      answer.add_flags_2[ra2][ra3][ra5][ra4] = 1 - answer.add_flags_2[ra2][ra3][ra5][ra4];
    }

    // スコア計算
    board.init_board();
    double tmp_score = create_ans(answer, board);

    // 焼きなましで採用判定
    double diff_score = (tmp_score - answer.current_score) * annealingParams.score_scale;
    double prob = exp(diff_score / temp);
    if (prob > rand_01()) {
      // 採用
      answer.current_score = tmp_score;

      // ベスト更新
      if (answer.current_score > best_answer.current_score) {
        best_answer.Copy(answer);
      }
    }
    else {
      // 元に戻す
      if (ra_exec_mode < annealingParams.operation_thresholds[0]) {
        // 近傍操作1 の巻き戻し
        answer.add_flags[ra1][ra2] = 1 - answer.add_flags[ra1][ra2];
      }
      else if (ra_exec_mode < annealingParams.operation_thresholds[1]) {
        // 近傍操作2 の巻き戻し
        answer.add_flags_2[ra3][ra4][ra5][ra2] = 1 - answer.add_flags_2[ra3][ra4][ra5][ra2];
      }
      else if (ra_exec_mode < annealingParams.operation_thresholds[2]) {
        // 近傍操作3 の巻き戻し
        answer.add_flags_2[ra2][ra3][ra5][ra4] = 1 - answer.add_flags_2[ra2][ra3][ra5][ra4];
      }
    }
  }

  if (exec_mode != 0 && exec_mode != 3) {
    cerr << loop << endl;
  }

  answer.Copy(best_answer);
}

static ll solve_case(int case_num, AnnealingParams annealingParams) {
  start_timer();

  Board board = input_data(case_num);

  ofstream ofs;
  open_ofs(case_num, ofs);

  Answer answer = build_initial_solution_2(board);
  answer.current_score = calculate_score(answer);

  // 焼きなまし実行
  run_simulated_annealing(board, answer, annealingParams);

  // 解答を出力
  output_data(answer, ofs);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (exec_mode != 0) {
    score = calculate_score(answer);
  }
  return score;
}

int main() {
  exec_mode = 2;

  AnnealingParams annealingParams;
  annealingParams.start_temperature[0] = 20.0;
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
  else if (exec_mode <= 2) {
    ll sum_score = 0;
    for (int i = 0; i < 15; ++i) {
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
  else if (exec_mode == 3) {
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
