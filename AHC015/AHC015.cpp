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
#define MAX_N 200005

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

namespace /* 乱数ライブラリ */
{
  static uint32_t xor_shift32()
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


  static double rand_unit()
  {
    return (xor_shift32() + 0.5) * (1.0 / UINT_MAX);
  }
}  // namespace

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
  Queue2D() : head(0), tail(0) {}

  void clear_queue()
  {
    head = 0;
    tail = 0;
  }

  int front_x() const
  {
    return arr[head][0];
  }

  int front_y() const
  {
    return arr[head][1];
  }

  void push(int x, int y)
  {
    arr[tail][0] = x;
    arr[tail][1] = y;
    tail++;
  }

  void pop()
  {
    head++;
  }

  int size() const
  {
    return tail - head;
  }
};
Queue2D que2d;

const int DX[4] = { -1, 0, 1, 0 };
const int DY[4] = { 0, -1, 0, 1 };
const char DC[4] = { 'F', 'L', 'B', 'R' };

const int N = 100;
const int H = 10;
const int W = 10;
int flavor_at_turn[N];
int input_index_list[N];
char output_dir_list[N];
int flavor_total[4];

using board_t = std::array<std::array<int, W>, H>;

board_t board{};
board_t work_board{};
board_t backup_board_buf{};

const double TL = 1.8;

int exec_mode = 0;

inline void copy_board_to_work()
{
  work_board = board;
}

inline void copy_work_to_board()
{
  board = work_board;
}

inline void backup_board()
{
  backup_board_buf = board;
}

inline void restore_board()
{
  board = backup_board_buf;
}

inline int random_empty_index(int turn)
{
  return xor_shift32() % (N - turn) + 1;
}

// 入力受け取り（実行中一度しか呼ばれないことを想定）
void load_testcase(int problemNum)
{
  string fileNameIfs = "./in/";
  string strNum;
  for (int i = 0; i < 4; ++i) {
    strNum += (char)(problemNum % 10 + '0');
    problemNum /= 10;
  }
  reverse(strNum.begin(), strNum.end());
  fileNameIfs += strNum + ".txt";

  ifstream ifs(fileNameIfs);

  for (int i = 0; i < N; ++i) {
    ifs >> flavor_at_turn[i];
  }
  for (int i = 0; i < N; ++i) {
    ifs >> input_index_list[i];
  }
}

// 解答出力
void save_answer(int problemNum)
{
  string fileNameOfs = "./out/";
  string strNum;
  for (int i = 0; i < 4; ++i) {
    strNum += (char)(problemNum % 10 + '0');
    problemNum /= 10;
  }
  reverse(strNum.begin(), strNum.end());
  fileNameOfs += strNum + ".txt";

  ofstream ofs(fileNameOfs);

  for (int i = 0; i < N; ++i) {
    ofs << output_dir_list[i] << endl;
  }

  ofs.close();
}

inline bool is_out_of_board(int x, int y)
{
  if (x < 0 || x >= 10 || y < 0 || y >= 10) {
    return true;
  }
  return false;
}

P find_empty_cell(int ite)
{
  int cnt = 1;
  for (int i = 0; i < H; ++i) {
    for (int j = 0; j < W; ++j) {
      if (board[i][j] == 0) {
        if (cnt == ite) {
          return P(i, j);
        }
        else {
          cnt++;
        }
      }
    }
  }

  cerr << "NG" << endl;
  return P(-1, -1);
}

int visited[H][W];
int visited_version;

int score_board(const board_t& board)
{
  //――― 事前初期化 ―――//
  que2d.clear_queue();
  double score_sum = 0.0;
  visited_version++;

  //――― 全セル走査 ―――//
  for (int start_row = 0; start_row < H; ++start_row) {
    for (int start_col = 0; start_col < W; ++start_col) {
      if (board[start_row][start_col] == 0) {
        continue;
      }
      if (visited[start_row][start_col] == visited_version) {
        continue;
      }

      que2d.push(start_row, start_col);
      visited[start_row][start_col] = visited_version;

      int component_size = 1;

      //――― BFS で同フレーバ連結成分を数える ―――//
      while (que2d.size()) {
        int cur_row = que2d.front_x();
        int cur_col = que2d.front_y();
        que2d.pop();

        for (int dir_idx = 0; dir_idx < 4; ++dir_idx) {
          int next_row = cur_row + DX[dir_idx];
          int next_col = cur_col + DY[dir_idx];

          if (is_out_of_board(next_row, next_col)) {
            continue;
          }
          if (board[next_row][next_col] == board[start_row][start_col] && visited[next_row][next_col] != visited_version) {
            visited[next_row][next_col] = visited_version;
            que2d.push(next_row, next_col);
            ++component_size;
          }
        }
      }

      score_sum += component_size * component_size;
    }
  }

  //――― 正規化 ―――//
  score_sum *= 1'000'000.0;

  int denom_sum = 0;
  for (int flavor = 1; flavor < 4; ++flavor) {
    denom_sum += flavor_total[flavor] * flavor_total[flavor];
  }

  score_sum /= denom_sum;
  return static_cast<int>(std::round(score_sum));
}

void tilt_board(int dir, board_t& board)
{
  if (dir == 0) {                // 上方向（Front）
    for (int col = 0; col < W; ++col) {
      int write_row = 0;
      for (int row = 0; row < H; ++row) {
        if (board[row][col] != 0) {
          int cell_value = board[row][col];
          board[row][col] = 0;
          board[write_row][col] = cell_value;
          ++write_row;
        }
      }
    }
  }
  else if (dir == 1) {           // 左方向（Left）
    for (int row = 0; row < H; ++row) {
      int write_col = 0;
      for (int col = 0; col < W; ++col) {
        if (board[row][col] != 0) {
          int cell_value = board[row][col];
          board[row][col] = 0;
          board[row][write_col] = cell_value;
          ++write_col;
        }
      }
    }
  }
  else if (dir == 2) {           // 下方向（Back）
    for (int col = 0; col < W; ++col) {
      int write_row = H - 1;
      for (int row = (H)-1; row >= 0; --row) {
        if (board[row][col] != 0) {
          int cell_value = board[row][col];
          board[row][col] = 0;
          board[write_row][col] = cell_value;
          --write_row;
        }
      }
    }
  }
  else if (dir == 3) {           // 右方向（Right）
    for (int row = 0; row < H; ++row) {
      int write_col = W - 1;
      for (int col = (W)-1; col >= 0; --col) {
        if (board[row][col] != 0) {
          int cell_value = board[row][col];
          board[row][col] = 0;
          board[row][write_col] = cell_value;
          --write_col;
        }
      }
    }
  }
}

int tilt_and_score_work(int dir)
{
  tilt_board(dir, work_board);
  return score_board(work_board);
}

int simulate_rollout(int start_turn, int sim_bias)
{
  int sim_turn = 10 + sim_bias;
  for (int turn = start_turn; turn < start_turn + sim_turn; ++turn) {
    if (turn == N) {
      break;
    }

    //――― キャンディー配置 ―――//
    int empty_index = random_empty_index(turn);
    P   cell = find_empty_cell(empty_index);
    int row = cell.first;
    int col = cell.second;

    if (board[row][col] != 0) {
      std::cout << "NG" << std::endl;
    }
    board[row][col] = flavor_at_turn[turn];

    //――― 4 方向シミュレートして最良を選択 ―――//
    int best_score = -1;
    int best_dir = -1;
    for (int dir = 0; dir < 4; ++dir) {
      copy_board_to_work();
      int score = tilt_and_score_work(dir);
      if (score > best_score) {
        best_score = score;
        best_dir = dir;
      }
    }

    tilt_board(best_dir, board);   // 本盤面を更新
  }

  return score_board(board);
}

int run_solver(int case_num)
{
  std::ofstream debug_ofs("debug.txt");

  start_timer();

  for (int i = 0; i < H; ++i) {
    for (int j = 0; j < W; ++j) {
      board[i][j] = 0;
    }
  }
  for (int i = 0; i < 4; ++i) {
    flavor_total[i] = 0;
  }

  //――― 入力 ―――//
  if (exec_mode == 0) {
    for (int i = 0; i < N; ++i) {
      std::cin >> flavor_at_turn[i];
    }
  }
  else {
    load_testcase(case_num);
  }
  for (int i = 0; i < N; ++i) {
    flavor_total[flavor_at_turn[i]]++;
  }

  int iter = 0;

  //――― 100 ターン回す ―――//
  for (int turn = 0; turn < N; ++turn) {
    /*---------- 空マスインデックス取得 ----------*/
    int empty_index;
    if (exec_mode == 0) {
      std::cin >> empty_index;
    }
    else {
      empty_index = input_index_list[turn];
    }

    P   cell = find_empty_cell(empty_index);   // (row,col)
    int row = cell.first;
    int col = cell.second;

    if (board[row][col] != 0) {
      cerr << "NG" << endl;
    }
    board[row][col] = flavor_at_turn[turn];

    backup_board();            // 現盤面を保存

    /*---------- モンテカルロ探索 ----------*/
    start_timer();
    const double turn_time_limit = TL / N;

    int    trial_cnt = 0;
    double dir_score_avg[4] = {};
    int    dir_trial_cnt[4] = {};
    int    sim_bias = 0;

    while (true) {
      if (trial_cnt % 8 == 0) {
        if (get_elapsed_time() > turn_time_limit) {
          break;
        }
      }
      iter++;
      if (trial_cnt % 4 == 0) {
        sim_bias = static_cast<int>(xor_shift32() % 3) - 1;
      }

      restore_board();

      int dir = trial_cnt % 4;
      tilt_board(dir, board);

      double score = simulate_rollout(turn + 1, sim_bias);

      dir_score_avg[dir] = (dir_score_avg[dir] * dir_trial_cnt[dir] + score) / (dir_trial_cnt[dir] + 1);
      dir_trial_cnt[dir]++;

      ++trial_cnt;
    }

    /*---------- ベスト方向の決定 ----------*/
    int best_dir = -1;
    int best_score = -1;
    for (int dir = 0; dir < 4; ++dir) {
      if (dir_score_avg[dir] > best_score) {
        best_score = static_cast<int>(dir_score_avg[dir]);
        best_dir = dir;
      }
    }

    /*---------- 出力 ----------*/
    if (exec_mode == 0) {
      std::cout << DC[best_dir] << '\n';
      std::fflush(stdout);
    }
    else {
      output_dir_list[turn] = DC[best_dir];
    }

    /*---------- 盤面を進める ----------*/
    restore_board();
    tilt_board(best_dir, board);

    if (exec_mode != 0) {
      debug_ofs << std::right << std::setw(7) << trial_cnt << ' ' << score_board(board) << '\n';
    }
  }

  /*---------- ファイル保存＆最終スコア表示 ----------*/
  if (exec_mode != 0) {
    save_answer(case_num);
    cout << "iter = " << iter << endl;
    std::cout << "Score = " << score_board(board) << '\n';
  }
  return 0;
}

int main()
{
  exec_mode = 2;

  if (exec_mode == 0) {
    run_solver(0);
  }
  else if (exec_mode <= 2) {
    ll sum_score = 0;
    for (int i = 0; i < 15; i++) {
      ll score = run_solver(i);
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
