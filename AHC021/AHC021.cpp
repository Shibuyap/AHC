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

// 乱数
namespace
{
  static uint32_t rand_xorshift()
  {
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

  static double rand_01()
  {
    return (rand_xorshift() + 0.5) * (1.0 / UINT_MAX);
  }

  static double rand_range(double l, double r)
  {
    return l + (r - l) * rand_01();
  }

  static uint32_t rand_range(uint32_t l, uint32_t r)
  {
    return l + rand_xorshift() % (r - l + 1); // [l, r]
  }

  void shuffle_array(int* arr, int n)
  {
    for (int i = n - 1; i >= 0; i--) {
      int j = rand_xorshift() % (i + 1);
      int swa = arr[i];
      arr[i] = arr[j];
      arr[j] = swa;
    }
  }
}

const ll INF = 1001001001001001001LL;
const int INT_INF = 1001001001;

const int DX[4] = { -1, 0, 1, 0 };
const int DY[4] = { 0, -1, 0, 1 };

const double TIME_LIMIT = 1.8;
int exec_mode;

constexpr int BOARD_SIZE = 30;
constexpr int MAX_LOOP = 10000;
constexpr int BALL_COUNT = BOARD_SIZE * (BOARD_SIZE + 1) / 2; // ＝BALL_COUNT

class State
{
public:
  int board[BOARD_SIZE][BOARD_SIZE];
  int ball_pos[1000][2];
  int move_cnt;
  int moves[MAX_LOOP][4];

  State()
  {
    move_cnt = 0;
  }

  int get_score()
  {
    return 100000 - 5 * move_cnt;
  }

  inline void swap_ball(int x1, int y1, int x2, int y2)
  {
    int ball1 = board[x1][y1];
    int ball2 = board[x2][y2];
    std::swap(ball_pos[ball1][0], ball_pos[ball2][0]);
    std::swap(ball_pos[ball1][1], ball_pos[ball2][1]);
    std::swap(board[x1][y1], board[x2][y2]);
  }

  inline void push_move(int x1, int y1, int x2, int y2)
  {
    moves[move_cnt][0] = x1;
    moves[move_cnt][1] = y1;
    moves[move_cnt][2] = x2;
    moves[move_cnt][3] = y2;
    swap_ball(x1, y1, x2, y2);
    ++move_cnt;
  }

  inline int get_diff(int x1, int y1, int x2, int y2)
  {
    return board[x1][y1] - board[x2][y2];
  }
};

bool is_out_of_range(int x, int y)
{
  //if (x < 0 || n <= x || y < 0 || n <= y) return true;
  return false;
}

State input_data(int case_num)
{
  State state;

  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  if (!ifs.is_open()) {
    // 標準入力
    for (int i = 0; i < BOARD_SIZE; i++) {
      for (int j = 0; j < i + 1; j++) {
        cin >> state.board[i][j];
      }
    }
  }
  else {
    // ファイル入力
    for (int i = 0; i < BOARD_SIZE; i++) {
      for (int j = 0; j < i + 1; j++) {
        ifs >> state.board[i][j];
      }
    }
  }

  return state;
}

void output_data(int case_num, const State& state)
{
  if (exec_mode == 0) {
    // 標準出力
    cout << state.move_cnt << endl;
    for (int i = 0; i < state.move_cnt; i++) {
      for (int j = 0; j < 4; j++) {
        cout << state.moves[i][j] << ' ';
      }
      cout << endl;
    }
  }
  else {
    // ファイル出力
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofstream ofs(oss.str());

    ofs << state.move_cnt << endl;
    for (int i = 0; i < state.move_cnt; i++) {
      for (int j = 0; j < 4; j++) {
        ofs << state.moves[i][j] << ' ';
      }
      ofs << endl;
    }

    if (ofs.is_open()) {
      ofs.close();
    }
  }
}

State greedy_swap_max_delta_with_tie(const State& initial_state, bool with_tie)
{
  State state = initial_state;

  while (true) {
    int x1, y1, x2, y2;
    int diff = 0;
    for (int i = 0; i < BOARD_SIZE - 1; i++) {
      for (int j = 0; j < i + 1; j++) {
        int diff1 = state.get_diff(i, j, i + 1, j);
        if (diff1 > diff || (diff1 == diff && with_tie)) {
          diff = diff1;
          x1 = i;
          y1 = j;
          x2 = i + 1;
          y2 = j;
        }
        int diff2 = state.get_diff(i, j, i + 1, j + 1);
        if (diff2 > diff || (diff2 == diff && with_tie)) {
          diff = diff2;
          x1 = i;
          y1 = j;
          x2 = i + 1;
          y2 = j + 1;
        }
      }
    }
    if (diff == 0) {
      break;
    }
    state.push_move(x1, y1, x2, y2);
  }

  return state;
}

ll solve_case(int case_num)
{
  start_timer();

  State initial_state = input_data(case_num);
  State best_state;

  State greedy1_state = greedy_swap_max_delta_with_tie(initial_state, false);
  best_state = greedy1_state;
  State greedy2_state = greedy_swap_max_delta_with_tie(initial_state, true);
  if (greedy2_state.get_score() > best_state.get_score()) {
    best_state = greedy2_state;
  }

  output_data(case_num, best_state);

  ll score = 0;
  if (exec_mode != 0) {
    score = best_state.get_score();
  }
  return score;
}

int main()
{
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
