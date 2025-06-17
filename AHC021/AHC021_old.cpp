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

#define srep(i, s, t) for (int i = s; i < t; ++i)

using namespace std;

typedef long long int ll;
typedef pair<int, int> P;
typedef pair<P, P> PP;

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
  static uint32_t Rand()
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


  static double Rand01()
  {
    return (Rand() + 0.5) * (1.0 / UINT_MAX);
  }
}  // namespace

constexpr int BOARD_SIZE = 30;
constexpr int MAX_LOOP = 10000;
constexpr int BALL_COUNT = BOARD_SIZE * (BOARD_SIZE + 1) / 2; // ＝BALL_COUNT
constexpr int LOCK_PREFIX = 52;

int board[BOARD_SIZE][BOARD_SIZE];
int ball_pos[1000][2];
int move_cnt;
int moves[MAX_LOOP][4];

int initial_board[BOARD_SIZE][BOARD_SIZE];
int best_move_cnt;
int best_moves[MAX_LOOP][4];
int global_best_move_cnt;
int global_best_moves[MAX_LOOP][4];

int saved_move_cnt;
int saved_moves[MAX_LOOP][4];

void init_board()
{
  for (int i = 0; i < BOARD_SIZE; ++i)
  {
    for (int j = 0; j < i + 1; ++j)
    {
      board[i][j] = initial_board[i][j];
      ball_pos[board[i][j]][0] = i;
      ball_pos[board[i][j]][1] = j;
    }
  }
}

void update_best_moves()
{
  if (move_cnt <= best_move_cnt) {
    best_move_cnt = move_cnt;
    for (int i = 0; i < move_cnt; ++i)
    {
      for (int j = 0; j < 4; ++j)
      {
        best_moves[i][j] = moves[i][j];
      }
    }
  }
}

void update_global_best_moves()
{
  if (best_move_cnt <= global_best_move_cnt) {
    global_best_move_cnt = best_move_cnt;
    for (int i = 0; i < best_move_cnt; ++i)
    {
      for (int j = 0; j < 4; ++j)
      {
        global_best_moves[i][j] = best_moves[i][j];
      }
    }
  }
}

void save_current_best()
{
  saved_move_cnt = best_move_cnt;
  for (int i = 0; i < best_move_cnt; ++i)
  {
    for (int j = 0; j < 4; ++j)
    {
      saved_moves[i][j] = best_moves[i][j];
    }
  }
}

void restore_best()
{
  move_cnt = best_move_cnt;
  for (int i = 0; i < move_cnt; ++i)
  {
    for (int j = 0; j < 4; ++j)
    {
      moves[i][j] = best_moves[i][j];
    }
  }
}

void restore_global_best()
{
  best_move_cnt = global_best_move_cnt;
  for (int i = 0; i < best_move_cnt; ++i)
  {
    for (int j = 0; j < 4; ++j)
    {
      best_moves[i][j] = global_best_moves[i][j];
    }
  }
}

inline void search_parent(const int x, const int y, int& nx, int& ny, int& diff, bool with_tie = false)
{
  if (y != 0 && board[x - 1][y - 1] - board[x][y] > diff || (board[x - 1][y - 1] - board[x][y] == diff && with_tie)) {
    diff = board[x - 1][y - 1] - board[x][y];
    nx = x - 1;
    ny = y - 1;
  }
  if (y != x && board[x - 1][y] - board[x][y] > diff || (board[x - 1][y] - board[x][y] == diff && with_tie)) {
    diff = board[x - 1][y] - board[x][y];
    nx = x - 1;
    ny = y;
  }
}

inline void search_parent_2(const int x, const int y, int& nx, int& ny, int& diff1, int& diff2)
{
  if (y != 0 && board[x - 1][y - 1] - board[x][y] > diff1) {
    diff1 = board[x - 1][y - 1] - board[x][y];
    nx = x - 1;
    ny = y - 1;
  }
  if (y != x && board[x - 1][y] - board[x][y] > diff2) {
    diff2 = board[x - 1][y] - board[x][y];
    nx = x - 1;
    ny = y;
  }
}

inline void swap_ball(int x1, int y1, int x2, int y2)
{
  int ball1 = board[x1][y1];
  int ball2 = board[x2][y2];
  std::swap(ball_pos[ball1][0], ball_pos[ball2][0]);
  std::swap(ball_pos[ball1][1], ball_pos[ball2][1]);
  std::swap(board[x1][y1], board[x2][y2]);
}

inline void push_move(int& loop, int x, int y, int nx, int ny)
{
  moves[loop][0] = x;
  moves[loop][1] = y;
  moves[loop][2] = nx;
  moves[loop][3] = ny;
  swap_ball(x, y, nx, ny);
  ++loop;
}

inline bool one_move(int& x, int& y, int& loop)
{
  int diff = 0;
  int nx = -1;
  int ny = -1;
  search_parent(x, y, nx, ny, diff);
  if (diff == 0) {
    return false;
  }
  push_move(loop, x, y, nx, ny);
  x = nx;
  y = ny;
  return true;
}

inline void some_operation(int& loop)
{
  for (int ball = 0; ball < BALL_COUNT; ++ball)
  {
    int x = -1, y = -1;
    x = ball_pos[ball][0];
    y = ball_pos[ball][1];

    while (loop < MAX_LOOP) {
      if (x == 0) {
        break;
      }
      bool result = one_move(x, y, loop);
      if (!result) {
        break;
      }
    }
  }
}

void greedy_swap_max_delta_with_tie(bool with_tie)
{
  init_board();

  int loop = 0;
  while (loop < MAX_LOOP) {
    int tmp[4] = {};
    int diff = 0;
    for (int i = 0; i < BOARD_SIZE - 1; ++i)
    {
      for (int j = 0; j < i + 1; ++j)
      {
        if (board[i][j] - board[i + 1][j] > diff || (board[i][j] - board[i + 1][j] == diff && with_tie)) {
          diff = board[i][j] - board[i + 1][j];
          tmp[0] = i;
          tmp[1] = j;
          tmp[2] = i + 1;
          tmp[3] = j;
        }
        if (board[i][j] - board[i + 1][j + 1] > diff || (board[i][j] - board[i + 1][j + 1] == diff && with_tie)) {
          diff = board[i][j] - board[i + 1][j + 1];
          tmp[0] = i;
          tmp[1] = j;
          tmp[2] = i + 1;
          tmp[3] = j + 1;
        }
      }
    }
    if (diff == 0) {
      break;
    }
    for (int j = 0; j < 4; ++j)
    {
      moves[loop][j] = tmp[j];
    }
    swap(board[tmp[0]][tmp[1]], board[tmp[2]][tmp[3]]);
    loop++;
  }

  move_cnt = loop;

  update_best_moves();
}

void ball_wise_ascent_greedy()
{
  init_board();

  int loop = 0;
  for (int ball = 0; ball < BALL_COUNT; ++ball)
  {
    int x = -1;
    int y = -1;
    for (int i = BOARD_SIZE - 1; i >= 0; --i)
    {
      for (int j = i; j >= 0; --j)
      {
        if (board[i][j] == ball) {
          x = i;
          y = j;
          break;
        }
      }
      if (x != -1) {
        break;
      }
    }
    while (loop < MAX_LOOP) {
      if (x == 0) {
        break;
      }
      bool result = one_move(x, y, loop);
      if (!result) {
        break;
      }
    }
  }
  move_cnt = loop;

  update_best_moves();
}

int random_prefix_len = 0;
void random_prefix_greedy(double time_limit)
{
  start_timer();

  int loopCount = 0;
  while (true) {
    if (get_elapsed_time() > time_limit) {
      break;
    }
    loopCount++;

    init_board();

    int loop = 0;

    int random = Rand() % 20 + 1;
    for (int _ = 0; _ < random; ++_)
    {
      while (true) {
        int x = Rand() % (BOARD_SIZE - 1);
        int y = Rand() % (x + 1);
        bool result = one_move(x, y, loop);
        if (result) {
          break;
        }
      }
    }

    some_operation(loop);

    move_cnt = loop;

    if (move_cnt <= best_move_cnt) {
      random_prefix_len = random;
    }
    update_best_moves();
  }
}

void adaptive_prefix_greedy(double time_limit)
{
  start_timer();

  int loopCount = 0;

  while (random_prefix_len < 50) {
    if (get_elapsed_time() > time_limit) {
      break;
    }
    loopCount++;

    init_board();

    int loop = 0;

    for (int _ = 0; _ < random_prefix_len; ++_)
    {
      int x = best_moves[loop][0];
      int y = best_moves[loop][1];
      int nx = best_moves[loop][2];
      int ny = best_moves[loop][3];
      push_move(loop, x, y, nx, ny);
    }

    while (true) {
      int x = Rand() % (BOARD_SIZE - 1);
      int y = Rand() % (x + 1);
      bool result = one_move(x, y, loop);
      if (result) {
        break;
      }
    }

    some_operation(loop);

    move_cnt = loop;

    if (move_cnt < best_move_cnt) {
      random_prefix_len++;
    }

    update_best_moves();
  }
}

inline int apply_prefix(int len)
{
  int loop = 0;
  for (int i = 0; i < len; ++i) {
    push_move(loop, best_moves[i][0], best_moves[i][1], best_moves[i][2], best_moves[i][3]);
  }
  return loop;
}

void prefix_lock_local_search(double time_limit)
{
  start_timer();

  best_move_cnt = saved_move_cnt;
  for (int i = 0; i < best_move_cnt; ++i)
  {
    for (int j = 0; j < 4; ++j)
    {
      best_moves[i][j] = saved_moves[i][j];
    }
  }

  int cnt = 0;
  while (true) {
    if (get_elapsed_time() > time_limit) {
      break;
    }
    cnt++;

    init_board();

    int ra = Rand() % (best_move_cnt - LOCK_PREFIX) + LOCK_PREFIX;
    int randomOpe = 0;
    int loop = apply_prefix(LOCK_PREFIX);

    for (int ball = 0; ball < BALL_COUNT; ++ball)
    {
      int x = -1, y = -1;
      x = ball_pos[ball][0];
      y = ball_pos[ball][1];

      while (loop < MAX_LOOP) {
        if (loop < ra) {
          if (best_moves[loop][0] == x && best_moves[loop][1] == y) {
            x = best_moves[loop][2];
            y = best_moves[loop][3];
            push_move(loop, best_moves[loop][0], best_moves[loop][1], best_moves[loop][2], best_moves[loop][3]);
            continue;
          }
          break;
        }

        if (x == 0) {
          break;
        }

        int diff1 = 0;
        int diff2 = 0;
        int nx = -1;
        int ny = -1;
        search_parent_2(x, y, nx, ny, diff1, diff2);
        if (diff1 == 0 && diff2 == 0) {
          break;
        }
        nx = x - 1;
        if (diff2 == 0) {
          ny = y - 1;
        }
        else if (diff1 == 0) {
          ny = y;
        }
        else {
          if (loop >= ra && randomOpe == 0) {
            if (diff1 > diff2) {
              ny = y;
            }
            else {
              ny = y - 1;
            }
            randomOpe = 1;
          }
          else {
            if (diff1 < diff2) {
              ny = y;
            }
            else {
              ny = y - 1;
            }
          }
        }
        push_move(loop, x, y, nx, ny);
        x = nx;
        y = ny;
      }
    }
    move_cnt = loop;

    update_best_moves();
  }

  update_global_best_moves();
}

bool random_local_search_v1_inner(int& loop, int& x, int& y, const int ball, const int ra, int& randomOpe)
{
  int diff1 = 0;
  int diff2 = 0;
  int nx = -1;
  int ny = -1;
  search_parent_2(x, y, nx, ny, diff1, diff2);
  if (diff1 == 0 && diff2 == 0) {
    return false;
  }
  nx = x - 1;
  if (diff2 == 0) {
    ny = y - 1;
  }
  else if (diff1 == 0) {
    ny = y;
  }
  else {
    if (loop >= ra && randomOpe == 0) {
      if (diff1 > diff2) {
        ny = y;
      }
      else {
        ny = y - 1;
      }
      randomOpe = 1;
    }
    else {
      if (diff1 < diff2) {
        ny = y;
      }
      else {
        ny = y - 1;
      }
    }
  }
  push_move(loop, x, y, nx, ny);
  x = nx;
  y = ny;

  return true;
}

bool random_local_search_v2_inner(int& loop, int& x, int& y, const int ball, const int ra, int& randomOpe)
{
  int nx = -1;
  int ny = -1;

  if (loop >= ra && randomOpe == 0) {
    nx = x;
    int dir = Rand() % 2;
    if (dir == 0) {
      if (y != 0) {
        ny = y - 1;
      }
      else {
        ny = y + 1;
      }
    }
    else {
      if (y != x) {
        ny = y + 1;
      }
      else {
        ny = y - 1;
      }
    }

    int ball2 = board[nx][ny];
    randomOpe = 1;
    if (ball2 > ball) {
      push_move(loop, x, y, nx, ny);
      x = nx;
      y = ny;
      return true;
    }
  }

  int diff1 = 0;
  int diff2 = 0;
  search_parent_2(x, y, nx, ny, diff1, diff2);
  if (diff1 == 0 && diff2 == 0) {
    return false;
  }
  nx = x - 1;
  if (diff1 < diff2) {
    ny = y;
  }
  else {
    ny = y - 1;
  }
  push_move(loop, x, y, nx, ny);
  x = nx;
  y = ny;

  return true;
}

void random_local_search(int version)
{
  init_board();

  int ra = Rand() % (best_move_cnt - LOCK_PREFIX) + LOCK_PREFIX;
  int randomOpe = 0;
  int loop = apply_prefix(LOCK_PREFIX);

  for (int ball = 0; ball < BALL_COUNT; ++ball)
  {
    int x = ball_pos[ball][0];
    int y = ball_pos[ball][1];
    while (loop < MAX_LOOP) {
      if (loop < ra) {
        if (best_moves[loop][0] == x && best_moves[loop][1] == y) {
          x = best_moves[loop][2];
          y = best_moves[loop][3];
          push_move(loop, best_moves[loop][0], best_moves[loop][1], best_moves[loop][2], best_moves[loop][3]);
          continue;
        }
        break;
      }

      if (x == 0) {
        break;
      }

      if (version == 1) {
        bool b = random_local_search_v1_inner(loop, x, y, ball, ra, randomOpe);
        if (!b) {
          break;
        }
      }
      else if (version == 2) {
        bool b = random_local_search_v2_inner(loop, x, y, ball, ra, randomOpe);
        if (!b) {
          break;
        }
      }
    }
  }
  move_cnt = loop;

  update_best_moves();
}

void load_input(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  // 標準入力する
  if (!ifs.is_open()) {
    for (int i = 0; i < BOARD_SIZE; ++i)
    {
      for (int j = 0; j < i + 1; ++j)
      {
        cin >> board[i][j];
        initial_board[i][j] = board[i][j];
      }
    }
  }
  // ファイル入力する
  else {
    for (int i = 0; i < BOARD_SIZE; ++i)
    {
      for (int j = 0; j < i + 1; ++j)
      {
        ifs >> board[i][j];
        initial_board[i][j] = board[i][j];
      }
    }
  }
}

void write_output(int mode, int problemNum)
{
  if (mode == 0) {
    cout << move_cnt << endl;
    for (int i = 0; i < move_cnt; ++i)
    {
      for (int j = 0; j < 4; ++j)
      {
        cout << moves[i][j] << ' ';
      }
      cout << endl;
    }
  }

  // ファイル出力
  if (mode != 0) {
    string fileNameOfs = "./out/";
    string strNum;
    for (int i = 0; i < 4; ++i)
    {
      strNum += (char)(problemNum % 10 + '0');
      problemNum /= 10;
    }
    reverse(strNum.begin(), strNum.end());
    fileNameOfs += strNum + ".txt";

    ofstream ofs(fileNameOfs);

    ofs << move_cnt << endl;
    for (int i = 0; i < move_cnt; ++i)
    {
      for (int j = 0; j < 4; ++j)
      {
        ofs << moves[i][j] << ' ';
      }
      ofs << endl;
    }
    ofs.close();
  }
}

int solve_single_problem(int mode, int probNum)
{
  load_input(probNum);

  move_cnt = MAX_LOOP;
  best_move_cnt = MAX_LOOP;
  global_best_move_cnt = MAX_LOOP;

  greedy_swap_max_delta_with_tie(false);
  greedy_swap_max_delta_with_tie(true);

  ball_wise_ascent_greedy();

  random_prefix_greedy(0.2);
  adaptive_prefix_greedy(0.3);

  save_current_best();

  for (int _ = 0; _ < 50; ++_)
  {
    prefix_lock_local_search(0.01);
  }

  restore_global_best();
  restore_best();

  start_timer();
  int loopCount = 0;
  while (true) {
    if (get_elapsed_time() > 0.8) {
      break;
    }
    loopCount++;
    if (Rand() % 2 == 0) {
      random_local_search(1);
    }
    else {
      random_local_search(2);
    }
  }

  // 戻して出力
  restore_best();

  if (mode != 0) {
    cout << probNum << ' ' << loopCount << ' ' << move_cnt << endl;
  }

  write_output(mode, probNum);
  return 100000 - 5 * move_cnt;
}

int main_old()
{
  int mode = 2;

  if (mode == 0) {
    solve_single_problem(mode, 0);
  }
  else if (mode == 1) {
    int probNum;
    cin >> probNum;
    solve_single_problem(mode, probNum);
  }
  else {
    int sum = 0;
    for (int _ = 0; _ < 1; ++_)
    {
      sum += solve_single_problem(1, _);
    }
    cout << sum << endl;
  }
  return 0;
}
