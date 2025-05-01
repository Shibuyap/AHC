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
// #include <atcoder/all>
#define rep(i, n) for (int i = 0; i < (n); ++i)
#define srep(i, s, t) for (int i = s; i < t; ++i)
#define drep(i, n) for (int i = (n)-1; i >= 0; --i)
using namespace std;
// using namespace atcoder;
typedef long long int ll;
typedef pair<int, int> P;
#define MAX_N 200005

namespace /* �������C�u���� */
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

// 2�����L���[�̃N���X
class Queue2D
{
private:
  static const int MAX_SIZE = 10000;
  int arr[MAX_SIZE][2];
  int head;
  int tail;

public:
  // �R���X�g���N�^
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

const int dir_row[4] = { -1, 0, 1, 0 };
const int dir_col[4] = { 0, -1, 0, 1 };
const char dir_char[4] = { 'F', 'L', 'B', 'R' };

const int n = 100;
const int H = 10;
const int W = 10;
int flavor_at_turn[n];
int input_index_list[n];
char output_dir_list[n];
int flavor_total[4];

using board_t = std::array<std::array<int, W>, H>;

board_t board{};            
board_t work_board{};
board_t backup_board_buf{};

const double TL = 1.8;

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

inline int random_empty_index(int turn) { return xor_shift32() % (n - turn) + 1; }

// ���͎󂯎��i���s����x�����Ă΂�Ȃ����Ƃ�z��j
void load_testcase(int problemNum)
{
  string fileNameIfs = "./in/";
  string strNum;
  rep(i, 4)
  {
    strNum += (char)(problemNum % 10 + '0');
    problemNum /= 10;
  }
  reverse(strNum.begin(), strNum.end());
  fileNameIfs += strNum + ".txt";

  ifstream ifs(fileNameIfs);

  rep(i, n) ifs >> flavor_at_turn[i];
  rep(i, n) ifs >> input_index_list[i];
}

// �𓚏o��
void save_answer(int problemNum)
{
  string fileNameOfs = "./out/";
  string strNum;
  rep(i, 4)
  {
    strNum += (char)(problemNum % 10 + '0');
    problemNum /= 10;
  }
  reverse(strNum.begin(), strNum.end());
  fileNameOfs += strNum + ".txt";

  ofstream ofs(fileNameOfs);

  rep(i, n) ofs << output_dir_list[i] << endl;

  ofs.close();
}

inline bool is_out_of_board(int x, int y)
{
  if (x < 0 || x >= 10 || y < 0 || y >= 10) return true;
  return false;
}

P find_empty_cell(int ite)
{
  int cnt = 1;
  rep(i, H)
  {
    rep(j, W)
    {
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

  cout << "NG" << endl;
  return P(-1, -1);
}

int visited[H][W];
int score_work_board()
{
  // �X�R�A�v�Z
  que2d.clear_queue();
  double res = 0;
  rep(i, H) rep(j, W) visited[i][j] = 0;
  rep(i, H)
  {
    rep(j, W)
    {
      if (work_board[i][j] == 0) continue;
      if (visited[i][j]) continue;
      que2d.push(i, j);
      visited[i][j] = 1;
      int sz = 1;
      while (que2d.size()) {
        int x = que2d.front_x();
        int y = que2d.front_y();
        que2d.pop();
        rep(k, 4)
        {
          int nx = x + dir_row[k];
          int ny = y + dir_col[k];
          if (is_out_of_board(nx, ny)) continue;
          if (work_board[nx][ny] == work_board[i][j] && visited[nx][ny] == 0) {
            visited[nx][ny] = 1;
            que2d.push(nx, ny);
            sz++;
          }
        }
      }

      res += sz * sz;
    }
  }

  res *= 1000000;
  int bo = 0;
  srep(i, 1, 4) { bo += flavor_total[i] * flavor_total[i]; }
  res /= bo;

  return round(res);
}

int score_board()
{ 
  // �X�R�A�v�Z
  que2d.clear_queue();
  double res = 0;
  rep(i, H) rep(j, W) visited[i][j] = 0;
  rep(i, H)
  {
    rep(j, W)
    {
      if (board[i][j] == 0) continue;
      if (visited[i][j]) continue;
      que2d.push(i, j);
      visited[i][j] = 1;
      int sz = 1;
      while (que2d.size()) {
        int x = que2d.front_x();
        int y = que2d.front_y();
        que2d.pop();
        rep(k, 4)
        {
          int nx = x + dir_row[k];
          int ny = y + dir_col[k];
          if (is_out_of_board(nx, ny)) continue;
          if (board[nx][ny] == board[i][j] && visited[nx][ny] == 0) {
            visited[nx][ny] = 1;
            que2d.push(nx, ny);
            sz++;
          }
        }
      }

      res += sz * sz;
    }
  }

  res *= 1000000;
  int bo = 0;
  srep(i, 1, 4) { bo += flavor_total[i] * flavor_total[i]; }
  res /= bo;

  return round(res);
}

void tilt_board(int dir, board_t& board)
{
  // �ꏊ�ړ�
  if (dir == 0) {
    rep(j, W)
    {
      int now = 0;
      rep(i, H)
      {
        if (board[i][j] != 0) {
          int tmp = board[i][j];
          board[i][j] = 0;
          board[now][j] = tmp;
          now++;
        }
      }
    }
  }
  else if (dir == 1) {
    rep(i, H)
    {
      int now = 0;
      rep(j, W)
      {
        if (board[i][j] != 0) {
          int tmp = board[i][j];
          board[i][j] = 0;
          board[i][now] = tmp;
          now++;
        }
      }
    }
  }
  else if (dir == 2) {
    rep(j, W)
    {
      int now = 9;
      drep(i, H)
      {
        if (board[i][j] != 0) {
          int tmp = board[i][j];
          board[i][j] = 0;
          board[now][j] = tmp;
          now--;
        }
      }
    }
  }
  else if (dir == 3) {
    rep(i, H)
    {
      int now = 9;
      drep(j, W)
      {
        if (board[i][j] != 0) {
          int tmp = board[i][j];
          board[i][j] = 0;
          board[i][now] = tmp;
          now--;
        }
      }
    }
  }
}

int tilt_and_score_work(int dir)
{
  tilt_board(dir, work_board);
  return score_work_board();
}

int simulate_rollout(int start_turn, int sim_bias)
{
  int sim_turn = 10 + sim_bias;
  srep(turn, start_turn, start_turn + sim_turn)
  {
    if (turn == n) break;

    //�\�\�\ �L�����f�B�[�z�u �\�\�\//
    int empty_index = random_empty_index(turn);
    P   cell = find_empty_cell(empty_index);
    int row = cell.first;
    int col = cell.second;

    if (board[row][col] != 0) {
      std::cout << "NG" << std::endl;
    }
    board[row][col] = flavor_at_turn[turn];

    //�\�\�\ 4 �����V�~�����[�g���čŗǂ�I�� �\�\�\//
    int best_score = -1;
    int best_dir = -1;
    rep(dir, 4)
    {
      copy_board_to_work();
      int score = tilt_and_score_work(dir);
      if (score > best_score) {
        best_score = score;
        best_dir = dir;
      }
    }

        tilt_board(best_dir, board);   // �{�Ֆʂ��X�V
  }

  return score_board();
}

int run_solver(int mode)
{
  std::ofstream debug_ofs("debug.txt");

  clock_t clk_start = clock();
  clock_t clk_end = clk_start;

  //�\�\�\ ���� �\�\�\//
  if (mode == 0) {
    rep(i, n) std::cin >> flavor_at_turn[i];
  }
  else {
    load_testcase(mode - 1);
  }
  rep(i, n) flavor_total[flavor_at_turn[i]]++;

  int iter = 0;

  //�\�\�\ 100 �^�[���� �\�\�\//
  rep(turn, n)
  {
    /*---------- ��}�X�C���f�b�N�X�擾 ----------*/
    int empty_index;
    if (mode == 0) {
      std::cin >> empty_index;
    }
    else {
      empty_index = input_index_list[turn];
    }

    P   cell = find_empty_cell(empty_index);   // (row,col)
    int row = cell.first;
    int col = cell.second;

    if (board[row][col] != 0) std::cout << "NG\n";
    board[row][col] = flavor_at_turn[turn];

    backup_board();            // ���Ֆʂ�ۑ�

    /*---------- �����e�J�����T�� ----------*/
    clk_start = clock();
    clk_end = clk_start;
    const double turn_time_limit = TL / 100.0;

    int    trial_cnt = 0;
    double dir_score_avg[4] = {};
    int    dir_trial_cnt[4] = {};
    int    sim_bias = 0;

    while (true) {
      if (trial_cnt % 8 == 0) {
        clk_end = clock();
        double elapsed =
          static_cast<double>(clk_end - clk_start) / CLOCKS_PER_SEC;
        if (elapsed > turn_time_limit) break;
      }
      iter++;
      if (trial_cnt % 4 == 0) {
        sim_bias = static_cast<int>(xor_shift32() % 3) - 1;
      }

      restore_board();

      int dir = trial_cnt % 4;
      tilt_board(dir, board);

      double score =
        simulate_rollout(turn + 1, sim_bias);

      dir_score_avg[dir] =
        (dir_score_avg[dir] * dir_trial_cnt[dir] + score) /
        (dir_trial_cnt[dir] + 1);
      dir_trial_cnt[dir]++;

      ++trial_cnt;
    }

    /*---------- �x�X�g�����̌��� ----------*/
    int best_dir = -1;
    int best_score = -1;
    rep(dir, 4)
    {
      if (dir_score_avg[dir] > best_score) {
        best_score = static_cast<int>(dir_score_avg[dir]);
        best_dir = dir;
      }
    }

    /*---------- �o�� ----------*/
    if (mode == 0) {
      std::cout << dir_char[best_dir] << '\n';
      std::fflush(stdout);
    }
    else {
      output_dir_list[turn] = dir_char[best_dir];
    }

    /*---------- �Ֆʂ�i�߂� ----------*/
    restore_board();
    tilt_board(best_dir, board);

    if (mode != 0) {
      debug_ofs << std::right << std::setw(7) << trial_cnt << ' ' << score_board() << '\n';
    }
  }

  /*---------- �t�@�C���ۑ����ŏI�X�R�A�\�� ----------*/
  if (mode != 0) {
    save_answer(mode - 1);
    cout << "iter = " << iter << endl;
    std::cout << "Score = " << score_board() << '\n';
  }
  return 0;
}

int main()
{
  int mode = 1;
  run_solver(mode);

  return 0;
}
