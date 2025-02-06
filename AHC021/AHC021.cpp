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
  static uint32_t randxor()
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

  // 0以上1未満の小数をとる乱数
  static double rand01() { return (randxor() + 0.5) * (1.0 / UINT_MAX); }
}  // namespace

int n = 30;
int b[30][30];
int argb[1000][2];
int ansSize;
int ans[11000][4];

int real_b[30][30];
int real_ansSize;
int real_ans[11000][4];
int real_real_ansSize;
int real_real_ans[11000][4];

int real_keep_ansSize;
int real_keep_ans[11000][4];

void InitB()
{
  rep(i, n)
  {
    rep(j, i + 1)
    {
      b[i][j] = real_b[i][j];
      argb[b[i][j]][0] = i;
      argb[b[i][j]][1] = j;
    }
  }
}

void UpdateAns()
{
  if (ansSize <= real_ansSize) {
    real_ansSize = ansSize;
    rep(i, ansSize)
    {
      rep(j, 4) { real_ans[i][j] = ans[i][j]; }
    }
  }
}

void UpdateRealAns()
{
  if (real_ansSize <= real_real_ansSize) {
    real_real_ansSize = real_ansSize;
    rep(i, real_ansSize)
    {
      rep(j, 4) { real_real_ans[i][j] = real_ans[i][j]; }
    }
  }
}

void KeepAns()
{
  real_keep_ansSize = real_ansSize;
  rep(i, real_ansSize)
  {
    rep(j, 4) { real_keep_ans[i][j] = real_ans[i][j]; }
  }
}

void RollBackAns()
{
  ansSize = real_ansSize;
  rep(i, ansSize)
  {
    rep(j, 4) { ans[i][j] = real_ans[i][j]; }
  }
}

void RollBackRealAns()
{
  real_ansSize = real_real_ansSize;
  rep(i, real_ansSize)
  {
    rep(j, 4) { real_ans[i][j] = real_real_ans[i][j]; }
  }
}

void Method1()
{
  int loop = 0;
  while (loop < 10000) {
    int tmp[4] = {};
    int diff = 0;
    rep(i, n - 1)
    {
      rep(j, i + 1)
      {
        if (b[i][j] - b[i + 1][j] > diff) {
          diff = b[i][j] - b[i + 1][j];
          tmp[0] = i;
          tmp[1] = j;
          tmp[2] = i + 1;
          tmp[3] = j;
        }
        if (b[i][j] - b[i + 1][j + 1] > diff) {
          diff = b[i][j] - b[i + 1][j + 1];
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
    rep(j, 4) { ans[loop][j] = tmp[j]; }
    swap(b[tmp[0]][tmp[1]], b[tmp[2]][tmp[3]]);
    loop++;
  }

  ansSize = loop;
  real_ansSize = ansSize;
  rep(i, ansSize)
  {
    rep(j, 4) { real_ans[i][j] = ans[i][j]; }
  }
}

void Method2()
{
  InitB();

  int loop = 0;
  while (loop < 10000) {
    int tmp[4] = {};
    int diff = 0;
    rep(i, n - 1)
    {
      rep(j, i + 1)
      {
        if (b[i][j] - b[i + 1][j] >= diff) {
          diff = b[i][j] - b[i + 1][j];
          tmp[0] = i;
          tmp[1] = j;
          tmp[2] = i + 1;
          tmp[3] = j;
        }
        if (b[i][j] - b[i + 1][j + 1] >= diff) {
          diff = b[i][j] - b[i + 1][j + 1];
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
    rep(j, 4) { ans[loop][j] = tmp[j]; }
    swap(b[tmp[0]][tmp[1]], b[tmp[2]][tmp[3]]);
    loop++;
  }

  ansSize = loop;

  UpdateAns();
}

void Method3()
{
  InitB();

  int loop = 0;
  rep(ball, 465)
  {
    int x = -1, y = -1;
    drep(i, n)
    {
      drep(j, i + 1)
      {
        if (b[i][j] == ball) {
          x = i;
          y = j;
          break;
        }
      }
      if (x != -1) break;
    }
    while (loop < 10000) {
      if (x == 0) break;
      int diff = 0;
      int nx = -1;
      int ny = -1;
      if (y != 0 && b[x - 1][y - 1] - b[x][y] > diff) {
        diff = b[x - 1][y - 1] - b[x][y];
        nx = x - 1;
        ny = y - 1;
      }
      if (y != x && b[x - 1][y] - b[x][y] > diff) {
        diff = b[x - 1][y] - b[x][y];
        nx = x - 1;
        ny = y;
      }
      if (diff == 0) break;
      ans[loop][0] = x;
      ans[loop][1] = y;
      ans[loop][2] = nx;
      ans[loop][3] = ny;
      swap(b[x][y], b[nx][ny]);
      x = nx;
      y = ny;
      loop++;
    }
  }
  ansSize = loop;

  UpdateAns();
}

int RandomCount = 0;
void Method3_2()
{
  clock_t startTime, endTime;
  startTime = clock();
  endTime = clock();

  int loopCount = 0;
  while (true) {
    endTime = clock();
    double nowTime = (double)(endTime - startTime) / CLOCKS_PER_SEC;
    if (nowTime > 0.2) {
      break;
    }
    loopCount++;

    InitB();

    int loop = 0;

    int random = randxor() % 20 + 1;
    rep(_, random)
    {
      while (true) {
        int x = randxor() % (n - 1);
        int y = randxor() % (x + 1);
        int diff = 0;
        int nx = -1;
        int ny = -1;
        if (y != 0 && b[x - 1][y - 1] - b[x][y] > diff) {
          diff = b[x - 1][y - 1] - b[x][y];
          nx = x - 1;
          ny = y - 1;
        }
        if (y != x && b[x - 1][y] - b[x][y] > diff) {
          diff = b[x - 1][y] - b[x][y];
          nx = x - 1;
          ny = y;
        }
        if (diff == 0) continue;
        ans[loop][0] = x;
        ans[loop][1] = y;
        ans[loop][2] = nx;
        ans[loop][3] = ny;
        int ball1 = b[x][y];
        int ball2 = b[nx][ny];
        swap(argb[ball1][0], argb[ball2][0]);
        swap(argb[ball1][1], argb[ball2][1]);
        swap(b[x][y], b[nx][ny]);
        x = nx;
        y = ny;
        loop++;
        break;
      }
    }

    rep(ball, 465)
    {
      int x = -1, y = -1;
      x = argb[ball][0];
      y = argb[ball][1];

      while (loop < 10000) {
        if (x == 0) break;
        int diff = 0;
        int nx = -1;
        int ny = -1;
        if (y != 0 && b[x - 1][y - 1] - b[x][y] > diff) {
          diff = b[x - 1][y - 1] - b[x][y];
          nx = x - 1;
          ny = y - 1;
        }
        if (y != x && b[x - 1][y] - b[x][y] > diff) {
          diff = b[x - 1][y] - b[x][y];
          nx = x - 1;
          ny = y;
        }
        if (diff == 0) break;
        ans[loop][0] = x;
        ans[loop][1] = y;
        ans[loop][2] = nx;
        ans[loop][3] = ny;
        int ball2 = b[ans[loop][2]][ans[loop][3]];
        argb[ball][0] = ans[loop][2];
        argb[ball][1] = ans[loop][3];
        argb[ball2][0] = ans[loop][0];
        argb[ball2][1] = ans[loop][1];
        swap(b[x][y], b[nx][ny]);
        x = nx;
        y = ny;
        loop++;
      }
    }
    ansSize = loop;

    if (ansSize <= real_ansSize) {
      RandomCount = random;
    }
    UpdateAns();
  }
}

void Method3_3()
{
  clock_t startTime, endTime;
  startTime = clock();
  endTime = clock();

  int loopCount = 0;

  while (RandomCount < 50) {
    endTime = clock();
    double nowTime = (double)(endTime - startTime) / CLOCKS_PER_SEC;
    if (nowTime > 0.3) {
      break;
    }
    loopCount++;

    InitB();

    int loop = 0;

    rep(_, RandomCount)
    {
      int x = real_ans[loop][0];
      int y = real_ans[loop][1];
      int nx = real_ans[loop][2];
      int ny = real_ans[loop][3];
      int ball1 = b[x][y];
      int ball2 = b[nx][ny];
      swap(argb[ball1][0], argb[ball2][0]);
      swap(argb[ball1][1], argb[ball2][1]);
      swap(b[x][y], b[nx][ny]);
      rep(j, 4) { ans[loop][j] = real_ans[loop][j]; }
      loop++;
    }

    while (true) {
      int x = randxor() % (n - 1);
      int y = randxor() % (x + 1);
      int diff = 0;
      int nx = -1;
      int ny = -1;
      if (y != 0 && b[x - 1][y - 1] - b[x][y] > diff) {
        diff = b[x - 1][y - 1] - b[x][y];
        nx = x - 1;
        ny = y - 1;
      }
      if (y != x && b[x - 1][y] - b[x][y] > diff) {
        diff = b[x - 1][y] - b[x][y];
        nx = x - 1;
        ny = y;
      }
      if (diff == 0) continue;
      ans[loop][0] = x;
      ans[loop][1] = y;
      ans[loop][2] = nx;
      ans[loop][3] = ny;
      int ball1 = b[x][y];
      int ball2 = b[nx][ny];
      swap(argb[ball1][0], argb[ball2][0]);
      swap(argb[ball1][1], argb[ball2][1]);
      swap(b[x][y], b[nx][ny]);
      x = nx;
      y = ny;
      loop++;
      break;
    }

    rep(ball, 465)
    {
      int x = -1, y = -1;
      x = argb[ball][0];
      y = argb[ball][1];

      while (loop < 10000) {
        if (x == 0) break;
        int diff = 0;
        int nx = -1;
        int ny = -1;
        if (y != 0 && b[x - 1][y - 1] - b[x][y] > diff) {
          diff = b[x - 1][y - 1] - b[x][y];
          nx = x - 1;
          ny = y - 1;
        }
        if (y != x && b[x - 1][y] - b[x][y] > diff) {
          diff = b[x - 1][y] - b[x][y];
          nx = x - 1;
          ny = y;
        }
        if (diff == 0) break;
        ans[loop][0] = x;
        ans[loop][1] = y;
        ans[loop][2] = nx;
        ans[loop][3] = ny;
        int ball2 = b[ans[loop][2]][ans[loop][3]];
        argb[ball][0] = ans[loop][2];
        argb[ball][1] = ans[loop][3];
        argb[ball2][0] = ans[loop][0];
        argb[ball2][1] = ans[loop][1];
        swap(b[x][y], b[nx][ny]);
        x = nx;
        y = ny;
        loop++;
      }
    }
    ansSize = loop;

    if (ansSize < real_ansSize) {
      RandomCount++;
    }

    UpdateAns();
  }
}

void Method4()
{
  clock_t startTime, endTime;
  startTime = clock();
  endTime = clock();

  real_ansSize = real_keep_ansSize;
  rep(i, real_ansSize)
  {
    rep(j, 4) { real_ans[i][j] = real_keep_ans[i][j]; }
  }

  int cnt = 0;
  while (true) {
    endTime = clock();
    double nowTime = (double)(endTime - startTime) / CLOCKS_PER_SEC;
    if (nowTime > 0.01) {
      break;
    }
    cnt++;

    InitB();

    int ra = rand() % (real_ansSize - 52) + 52;
    int randomOpe = 0;
    int loop = 0;

    rep(_, 52)
    {
      int x = real_ans[loop][0];
      int y = real_ans[loop][1];
      int nx = real_ans[loop][2];
      int ny = real_ans[loop][3];
      int ball1 = b[x][y];
      int ball2 = b[nx][ny];
      swap(argb[ball1][0], argb[ball2][0]);
      swap(argb[ball1][1], argb[ball2][1]);
      swap(b[x][y], b[nx][ny]);
      rep(j, 4) { ans[loop][j] = real_ans[loop][j]; }
      loop++;
    }

    rep(ball, 465)
    {
      int x = -1, y = -1;
      x = argb[ball][0];
      y = argb[ball][1];

      while (loop < 10000) {
        if (loop < ra) {
          if (real_ans[loop][0] == x && real_ans[loop][1] == y) {
            rep(j, 4) { ans[loop][j] = real_ans[loop][j]; }
            x = ans[loop][2];
            y = ans[loop][3];

            int ball2 = b[ans[loop][2]][ans[loop][3]];
            argb[ball][0] = ans[loop][2];
            argb[ball][1] = ans[loop][3];
            argb[ball2][0] = ans[loop][0];
            argb[ball2][1] = ans[loop][1];
            swap(b[ans[loop][0]][ans[loop][1]], b[ans[loop][2]][ans[loop][3]]);
            loop++;
            continue;
          }
          else {
            break;
          }
        }
        if (x == 0) break;
        int diff1 = 0;
        int diff2 = 0;
        int nx = -1;
        int ny = -1;
        if (y != 0 && b[x - 1][y - 1] - b[x][y] > diff1) {
          diff1 = b[x - 1][y - 1] - b[x][y];
          nx = x - 1;
          ny = y - 1;
        }
        if (y != x && b[x - 1][y] - b[x][y] > diff2) {
          diff2 = b[x - 1][y] - b[x][y];
          nx = x - 1;
          ny = y;
        }
        if (diff1 == 0 && diff2 == 0) break;
        if (diff2 == 0) {
          nx = x - 1;
          ny = y - 1;
        }
        else if (diff1 == 0) {
          nx = x - 1;
          ny = y;
        }
        else {
          if (loop >= ra && randomOpe == 0) {
            if (diff1 > diff2) {
              nx = x - 1;
              ny = y;
            }
            else {
              nx = x - 1;
              ny = y - 1;
            }
            randomOpe = 1;
          }
          else {
            if (diff1 < diff2) {
              nx = x - 1;
              ny = y;
            }
            else {
              nx = x - 1;
              ny = y - 1;
            }
          }
        }
        ans[loop][0] = x;
        ans[loop][1] = y;
        ans[loop][2] = nx;
        ans[loop][3] = ny;
        int ball2 = b[ans[loop][2]][ans[loop][3]];
        argb[ball][0] = ans[loop][2];
        argb[ball][1] = ans[loop][3];
        argb[ball2][0] = ans[loop][0];
        argb[ball2][1] = ans[loop][1];
        swap(b[x][y], b[nx][ny]);
        x = nx;
        y = ny;
        loop++;
      }
    }
    ansSize = loop;

    UpdateAns();
  }

  UpdateRealAns();
}

void Method5_1()
{
  InitB();

  int ra = rand() % (real_ansSize - 52) + 52;
  int randomOpe = 0;
  int loop = 0;

  rep(_, 52)
  {
    int x = real_ans[loop][0];
    int y = real_ans[loop][1];
    int nx = real_ans[loop][2];
    int ny = real_ans[loop][3];
    int ball1 = b[x][y];
    int ball2 = b[nx][ny];
    swap(argb[ball1][0], argb[ball2][0]);
    swap(argb[ball1][1], argb[ball2][1]);
    swap(b[x][y], b[nx][ny]);
    rep(j, 4) { ans[loop][j] = real_ans[loop][j]; }
    loop++;
  }

  rep(ball, 465)
  {
    int x = -1, y = -1;
    x = argb[ball][0];
    y = argb[ball][1];
    while (loop < 10000) {
      if (loop < ra) {
        if (real_ans[loop][0] == x && real_ans[loop][1] == y) {
          rep(j, 4) { ans[loop][j] = real_ans[loop][j]; }
          x = ans[loop][2];
          y = ans[loop][3];
          int ball2 = b[ans[loop][2]][ans[loop][3]];
          argb[ball][0] = ans[loop][2];
          argb[ball][1] = ans[loop][3];
          argb[ball2][0] = ans[loop][0];
          argb[ball2][1] = ans[loop][1];
          swap(b[ans[loop][0]][ans[loop][1]], b[ans[loop][2]][ans[loop][3]]);
          loop++;
          continue;
        }
        else {
          break;
        }
      }
      if (x == 0) break;
      int diff1 = 0;
      int diff2 = 0;
      int nx = -1;
      int ny = -1;
      if (y != 0 && b[x - 1][y - 1] - b[x][y] > diff1) {
        diff1 = b[x - 1][y - 1] - b[x][y];
        nx = x - 1;
        ny = y - 1;
      }
      if (y != x && b[x - 1][y] - b[x][y] > diff2) {
        diff2 = b[x - 1][y] - b[x][y];
        nx = x - 1;
        ny = y;
      }
      if (diff1 == 0 && diff2 == 0) break;
      if (diff2 == 0) {
        nx = x - 1;
        ny = y - 1;
      }
      else if (diff1 == 0) {
        nx = x - 1;
        ny = y;
      }
      else {
        if (loop >= ra && randomOpe == 0) {
          if (diff1 > diff2) {
            nx = x - 1;
            ny = y;
          }
          else {
            nx = x - 1;
            ny = y - 1;
          }
          randomOpe = 1;
        }
        else {
          if (diff1 < diff2) {
            nx = x - 1;
            ny = y;
          }
          else {
            nx = x - 1;
            ny = y - 1;
          }
        }
      }
      ans[loop][0] = x;
      ans[loop][1] = y;
      ans[loop][2] = nx;
      ans[loop][3] = ny;
      int ball2 = b[ans[loop][2]][ans[loop][3]];
      argb[ball][0] = ans[loop][2];
      argb[ball][1] = ans[loop][3];
      argb[ball2][0] = ans[loop][0];
      argb[ball2][1] = ans[loop][1];
      swap(b[x][y], b[nx][ny]);
      x = nx;
      y = ny;
      loop++;
    }
  }
  ansSize = loop;

  UpdateAns();
}

void Method5_2()
{
  InitB();

  int ra = rand() % (real_ansSize - 52) + 52;
  int randomOpe = 0;
  int loop = 0;

  rep(_, 52)
  {
    int x = real_ans[loop][0];
    int y = real_ans[loop][1];
    int nx = real_ans[loop][2];
    int ny = real_ans[loop][3];
    int ball1 = b[x][y];
    int ball2 = b[nx][ny];
    swap(argb[ball1][0], argb[ball2][0]);
    swap(argb[ball1][1], argb[ball2][1]);
    swap(b[x][y], b[nx][ny]);
    rep(j, 4) { ans[loop][j] = real_ans[loop][j]; }
    loop++;
  }

  rep(ball, 465)
  {
    int x = -1, y = -1;
    x = argb[ball][0];
    y = argb[ball][1];
    while (loop < 10000) {
      if (loop < ra) {
        if (real_ans[loop][0] == x && real_ans[loop][1] == y) {
          rep(j, 4) { ans[loop][j] = real_ans[loop][j]; }
          x = ans[loop][2];
          y = ans[loop][3];
          int ball2 = b[ans[loop][2]][ans[loop][3]];
          argb[ball][0] = ans[loop][2];
          argb[ball][1] = ans[loop][3];
          argb[ball2][0] = ans[loop][0];
          argb[ball2][1] = ans[loop][1];
          swap(b[ans[loop][0]][ans[loop][1]], b[ans[loop][2]][ans[loop][3]]);
          loop++;
          continue;
        }
        else {
          break;
        }
      }
      if (x == 0) break;

      int diff1 = 0;
      int diff2 = 0;
      int nx = -1;
      int ny = -1;

      if (loop >= ra && randomOpe == 0) {
        int dir = randxor() % 2;
        if (dir == 0) {
          if (y != 0) {
            nx = x;
            ny = y - 1;
          }
          else {
            nx = x;
            ny = y + 1;
          }
        }
        else {
          if (y != x) {
            nx = x;
            ny = y + 1;
          }
          else {
            nx = x;
            ny = y - 1;
          }
        }

        ans[loop][0] = x;
        ans[loop][1] = y;
        ans[loop][2] = nx;
        ans[loop][3] = ny;
        int ball2 = b[ans[loop][2]][ans[loop][3]];
        randomOpe = 1;
        if (ball2 > ball) {
          argb[ball][0] = ans[loop][2];
          argb[ball][1] = ans[loop][3];
          argb[ball2][0] = ans[loop][0];
          argb[ball2][1] = ans[loop][1];
          swap(b[x][y], b[nx][ny]);
          x = nx;
          y = ny;
          loop++;
          continue;
        }
      }

      if (y != 0 && b[x - 1][y - 1] - b[x][y] > diff1) {
        diff1 = b[x - 1][y - 1] - b[x][y];
        nx = x - 1;
        ny = y - 1;
      }
      if (y != x && b[x - 1][y] - b[x][y] > diff2) {
        diff2 = b[x - 1][y] - b[x][y];
        nx = x - 1;
        ny = y;
      }
      if (diff1 == 0 && diff2 == 0) break;
      if (diff2 == 0) {
        nx = x - 1;
        ny = y - 1;
      }
      else if (diff1 == 0) {
        nx = x - 1;
        ny = y;
      }
      else {
        if (diff1 < diff2) {
          nx = x - 1;
          ny = y;
        }
        else {
          nx = x - 1;
          ny = y - 1;
        }
      }
      ans[loop][0] = x;
      ans[loop][1] = y;
      ans[loop][2] = nx;
      ans[loop][3] = ny;
      int ball2 = b[ans[loop][2]][ans[loop][3]];
      argb[ball][0] = ans[loop][2];
      argb[ball][1] = ans[loop][3];
      argb[ball2][0] = ans[loop][0];
      argb[ball2][1] = ans[loop][1];
      swap(b[x][y], b[nx][ny]);
      x = nx;
      y = ny;
      loop++;
    }
  }
  ansSize = loop;

  UpdateAns();
}

void Input(int problemNum)
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

  // 標準入力する
  if (!ifs.is_open()) {
    rep(i, n)
    {
      rep(j, i + 1)
      {
        cin >> b[i][j];
        real_b[i][j] = b[i][j];
      }
    }
  }
  // ファイル入力する
  else {
    rep(i, n)
    {
      rep(j, i + 1)
      {
        ifs >> b[i][j];
        real_b[i][j] = b[i][j];
      }
    }
  }
}

void Output(int mode, int problemNum)
{
  if (mode == 0) {
    cout << ansSize << endl;
    rep(i, ansSize)
    {
      rep(j, 4) { cout << ans[i][j] << ' '; }
      cout << endl;
    }
  }

  // ファイル出力
  if (mode != 0) {
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

    ofs << ansSize << endl;
    rep(i, ansSize)
    {
      rep(j, 4) { ofs << ans[i][j] << ' '; }
      ofs << endl;
    }
    ofs.close();
  }
}

int Solve(int mode, int probNum)
{
  // 乱数調整
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    randxor();
  }

  Input(probNum);

  ansSize = 10000;
  real_ansSize = 10000;
  real_real_ansSize = 10000;

  Method1();

  Method2();

  Method3();

  Method3_2();
  Method3_3();

  KeepAns();

  rep(_, 50) { Method4(); }

  RollBackRealAns();
  RollBackAns();

  clock_t startTime, endTime;
  startTime = clock();
  endTime = clock();
  int loopCount = 0;
  while (true) {
    endTime = clock();
    double nowTime = (double)(endTime - startTime) / CLOCKS_PER_SEC;
    if (nowTime > 0.8) {
      break;
    }
    loopCount++;
    if (randxor() % 2 == 0) {
      Method5_1();
    }
    else {
      Method5_2();
    }
  }

  // 戻して出力
  RollBackAns();

  if (mode != 0) {
    cout << probNum << ' ' << loopCount << ' ' << ansSize << endl;
  }

  Output(mode, probNum);
  return 100000 - 5 * ansSize;
}

int main()
{
  int mode = 0;

  if (mode == 0) {
    Solve(mode, 0);
  }
  else if (mode == 1) {
    int probNum;
    cin >> probNum;
    Solve(mode, probNum);
  }
  else {
    int sum = 0;
    rep(_, 15) { sum += Solve(1, _); }
    cout << sum << endl;
  }
}