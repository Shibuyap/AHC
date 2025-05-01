#include <algorithm>
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

int board[H][W];
int work_board[H][W];
int backup_board_buf[H][W];
const double TL = 1.8;

inline void copy_board_to_work()
{
  rep(i, H)
  {
    rep(j, W)
    {
      work_board[i][j] = board[i][j];
    }
  }
}
inline void copy_work_to_board()
{
  rep(i, H)
  {
    rep(j, W)
    {
      board[i][j] = work_board[i][j];
    }
  }
}
inline void backup_board()
{
  rep(i, H)
  {
    rep(j, W)
    {
      backup_board_buf[i][j] = board[i][j];
    }
  }
}
inline void restore_board()
{
  rep(i, H)
  {
    rep(j, W)
    {
      board[i][j] = backup_board_buf[i][j];
    }
  }
}

inline int random_empty_index(int turn) { return xor_shift32() % (n - turn) + 1; }

// 入力受け取り（実行中一度しか呼ばれないことを想定）
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

// 解答出力
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
  // スコア計算
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
{  // スコア計算
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

void tilt_work_board(int dir)
{
  // 場所移動
  if (dir == 0) {
    rep(j, W)
    {
      int now = 0;
      rep(i, H)
      {
        if (work_board[i][j] != 0) {
          int tmp = work_board[i][j];
          work_board[i][j] = 0;
          work_board[now][j] = tmp;
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
        if (work_board[i][j] != 0) {
          int tmp = work_board[i][j];
          work_board[i][j] = 0;
          work_board[i][now] = tmp;
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
        if (work_board[i][j] != 0) {
          int tmp = work_board[i][j];
          work_board[i][j] = 0;
          work_board[now][j] = tmp;
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
        if (work_board[i][j] != 0) {
          int tmp = work_board[i][j];
          work_board[i][j] = 0;
          work_board[i][now] = tmp;
          now--;
        }
      }
    }
  }
}

void tilt_board(int dir)
{
  // 場所移動
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
  tilt_work_board(dir);
  return score_work_board();
}

int simulate_rollout(int startTurn, int simBias)
{
  int simTurn = 10 + simBias;
  srep(_, startTurn, startTurn + simTurn)
  {
    if (_ == n) break;
    // 場所決め
    int po = random_empty_index(_);
    P p = find_empty_cell(po);
    int wx = p.first;
    int wy = p.second;
    if (board[wx][wy] != 0) {
      cout << "NG" << endl;
    }
    board[wx][wy] = flavor_at_turn[_];

    int tmpMax = -1;
    int tmpAns = -1;
    rep(i, 4)
    {
      copy_board_to_work();
      int tmp = tilt_and_score_work(i);
      if (tmp > tmpMax) {
        tmpMax = tmp;
        tmpAns = i;
      }
    }

    tilt_board(tmpAns);
  }

  return score_board();
}

int run_solver(int mode)
{
  ofstream ofs2("debug.txt");
  clock_t startTime, endTime;
  startTime = clock();
  endTime = clock();

  if (mode == 0) {
    rep(i, n) { cin >> flavor_at_turn[i]; }
  }
  else {
    load_testcase(mode - 1);
  }

  rep(i, n) { flavor_total[flavor_at_turn[i]]++; }

  rep(_, n)
  {
    int ite;
    if (mode == 0) {
      cin >> ite;
    }
    else {
      ite = input_index_list[_];
    }

    P p = find_empty_cell(ite);
    int wx = p.first;
    int wy = p.second;
    if (board[wx][wy] != 0) {
      cout << "NG" << endl;
    }
    board[wx][wy] = flavor_at_turn[_];

    backup_board();

    // ここで探索
    startTime = clock();
    endTime = clock();
    double tl = TL / 100;
    int loop = 0;

    double ave[4] = {};
    int cnt[4] = {};
    int simBias = 0;
    while (true) {
      if (loop % 8 == 0) {
        endTime = clock();
        double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
        if (nowTime > tl) break;
      }

      if (loop % 4 == 0) {
        simBias = xor_shift32() % 3 - 1;
      }

      restore_board();

      int dir = loop % 4;
      tilt_board(dir);
      ave[dir] = (ave[dir] * cnt[dir] + simulate_rollout(_ + 1, simBias)) / (cnt[dir] + 1);
      cnt[dir]++;

      loop++;
    }

    int tmpMax = -1;
    int tmpAns = -1;
    rep(i, 4)
    {
      if (ave[i] > tmpMax) {
        tmpMax = ave[i];
        tmpAns = i;
      }
    }

    if (mode == 0) {
      std::cout << dir_char[tmpAns] << endl;
      std::fflush(stdout);
    }
    else {
      output_dir_list[_] = dir_char[tmpAns];
    }

    restore_board();
    tilt_board(tmpAns);

    if (mode != 0) {
      ofs2 << std::right << std::setw(7) << loop << ' ' << score_board() << endl;
    }
  }

  if (mode != 0) {
    save_answer(mode - 1);
  }

  if (mode != 0) {
    cout << "Score = " << score_board() << endl;
  }

  return 0;
}

int main()
{
  int mode = 1;
  run_solver(mode);

  return 0;
}
