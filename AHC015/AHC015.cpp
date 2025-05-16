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
  static const int MAX_SIZE = 500;
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

char output_dir_list[N];


int visited[H][W];
int visited_version;

class Input
{
private:
  int mode;
  int turn;
  int input_index_list[N];
public:
  int flavor_total[4];
  int flavor_at_turn[N];

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

  Input(int _mode, int _case_num) : mode(_mode), turn(0)
  {
    if (mode == 0) {
      for (int i = 0; i < N; ++i) {
        std::cin >> flavor_at_turn[i];
      }
    }
    else {
      load_testcase(_case_num);
    }

    for (int i = 0; i < 4; ++i) {
      flavor_total[i] = 0;
    }
    for (int i = 0; i < N; ++i) {
      flavor_total[flavor_at_turn[i]]++;
    }
  }

  int get_empty_index()
  {
    int empty_index = 0;
    if (mode == 0) {
      std::cin >> empty_index;
    }
    else {
      empty_index = input_index_list[turn];
    }
    turn++;
    return empty_index;
  }
};

class Board
{
public:
  std::array<std::array<int, W>, H> board;

  void init()
  {
    for (int i = 0; i < H; ++i) {
      for (int j = 0; j < W; ++j) {
        board[i][j] = 0;
      }
    }
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

  inline bool is_out_of_board(int x, int y)
  {
    if (x < 0 || x >= H || y < 0 || y >= W) {
      return true;
    }
    return false;
  }

  int calculate_score(const Input& input)
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
      denom_sum += input.flavor_total[flavor] * input.flavor_total[flavor];
    }

    score_sum /= denom_sum;
    return static_cast<int>(std::round(score_sum));
  }

  void tilt_board(int dir)
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

  void set_candy(int empty_index, int flavor)
  {
    P   cell = find_empty_cell(empty_index);
    int row = cell.first;
    int col = cell.second;
    board[row][col] = flavor;
  }
};

const double TL = 1.8;

int exec_mode = 0;

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

// ルールベースに従う操作列を記録する配列
vector<int> rule_actions(N, -1);

// イチゴは後ろに、スイカは右前に、パンプキンは左前に集める
vector<vector<int>> rule_action_samples = {
    {0, 2, 2},
    {0, 1, 3},
    {0, 1, 3},
};

void init_rule_actions(const Input& input)
{
  for (int t = 0; t < N - 1; t++) {
    int now = input.flavor_at_turn[t];
    int next = input.flavor_at_turn[t + 1];
    rule_actions[t] = rule_action_samples[now - 1][next - 1];
  }
  rule_actions[N - 1] = 0; // 最後のターンは何でもよい
}

inline int generate_random_empty_index(int turn)
{
  return xor_shift32() % (N - turn) + 1;
}

int simulate_rollout_with_rules(int start_turn, Board& board, const Input& input)
{
  for (int turn = start_turn; turn < N; ++turn) {
    int empty_index = generate_random_empty_index(turn);
    board.set_candy(empty_index, input.flavor_at_turn[turn]);
    board.tilt_board(rule_actions[turn]);
  }
  return board.calculate_score(input);
}

void monte_carlo(Board& board, Input& input)
{
  init_rule_actions(input);

  int iter = 0;

  //――― 100 ターン回す ―――//
  for (int turn = 0; turn < N; ++turn) {
    int empty_index = input.get_empty_index();
    board.set_candy(empty_index, input.flavor_at_turn[turn]);

    Board backup_board = board;

    /*---------- モンテカルロ探索 ----------*/
    start_timer();
    const double turn_time_limit = TL / N;

    int    trial_cnt = 0;
    double dir_score_avg[4] = {};
    int    dir_trial_cnt[4] = {};
    int sim_turn = 100;

    while (true) {
      if (trial_cnt % 8 == 0) {
        if (get_elapsed_time() > turn_time_limit) {
          break;
        }
      }

      iter++;

      board = backup_board;

      int dir = trial_cnt % 4;
      board.tilt_board(dir);

      double score = simulate_rollout_with_rules(turn + 1, board, input);

      dir_score_avg[dir] = (dir_score_avg[dir] * dir_trial_cnt[dir] + score) / (dir_trial_cnt[dir] + 1);
      dir_trial_cnt[dir]++;

      ++trial_cnt;
    }

    int best_dir = -1;
    int best_score = -1;
    for (int dir = 0; dir < 4; ++dir) {
      if (dir_score_avg[dir] > best_score) {
        best_score = static_cast<int>(dir_score_avg[dir]);
        best_dir = dir;
      }
    }

    if (exec_mode == 0) {
      std::cout << DC[best_dir] << '\n';
      std::fflush(stdout);
    }
    else {
      output_dir_list[turn] = DC[best_dir];
    }

    board = backup_board;
    board.tilt_board(best_dir);
  }

  if (exec_mode != 0) {
    cout << "iter = " << iter << endl;
  }
}

int run_solver(int case_num)
{
  start_timer();

  Input input(exec_mode, case_num);

  Board board;
  board.init();

  monte_carlo(board, input);  

  int score = 0;
  if (exec_mode != 0) {
    save_answer(case_num);
    score = board.calculate_score(input);
    std::cout << "Score = " << score << '\n';
  }

  return score;
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
