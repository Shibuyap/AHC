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
#define rep(i, n) for (int i = 0; i < (n); ++i)
#define srep(i, s, t) for (int i = s; i < t; ++i)
#define drep(i, n) for (int i = (n)-1; i >= 0; --i)
using namespace std;
typedef long long int ll;
typedef pair<int, int> P;

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

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

double TL = 1.8;
int mode;

const int n = 50;
const int m = 100;
int c[52][52];
int d[52][52];
int g[101][101];
int cntM[101];
double gx[101], gy[101];
int zero;

// 複数ケース回すときに内部状態を初期値に戻す
void SetUp()
{
  zero = 0;
  rep(i, 101)
  {
    rep(j, 101) { g[i][j] = 0; }
  }
  rep(i, 101)
  {
    cntM[i] = 0;
    gx[i] = 0;
    gy[i] = 0;
  }
}

// 入力受け取り
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
    int nn, mm;
    cin >> nn >> mm;
    srep(i, 1, n + 1)
    {
      srep(j, 1, n + 1) { cin >> c[i][j]; }
    }
  }
  // ファイル入力する
  else {
    int nn, mm;
    ifs >> nn >> mm;
    srep(i, 1, n + 1)
    {
      srep(j, 1, n + 1) { ifs >> c[i][j]; }
    }
  }

  rep(i, n + 2)
  {
    rep(j, n + 1)
    {
      g[c[i][j]][c[i][j + 1]]++;
      g[c[i][j + 1]][c[i][j]]++;
    }
  }
  rep(i, n + 1)
  {
    rep(j, n + 2)
    {
      g[c[i][j]][c[i + 1][j]]++;
      g[c[i + 1][j]][c[i][j]]++;
    }
  }
  rep(i, n + 2)
  {
    rep(j, n + 2)
    {
      int z = c[i][j];
      cntM[z]++;
      gx[z] += i;
      gy[z] += j;
    }
  }
}

// 出力ファイルストリームオープン
void OpenOfs(int probNum, ofstream& ofs)
{
  if (mode != 0) {
    string fileNameOfs = "./out/";
    string strNum;
    rep(i, 4)
    {
      strNum += (char)(probNum % 10 + '0');
      probNum /= 10;
    }
    reverse(strNum.begin(), strNum.end());
    fileNameOfs += strNum + ".txt";

    ofs.open(fileNameOfs);
  }
}

// スコア計算
int CalcScore() { return zero + 1; }

// 初期解生成
void Initialize()
{
  rep(i, n + 2)
  {
    rep(j, n + 2) { d[i][j] = c[i][j]; }
  }
}

// 解答出力
void Output(ofstream& ofs)
{
  if (mode == 0) {
    srep(i, 1, n + 1)
    {
      srep(j, 1, n + 1) { cout << d[i][j] << ' '; }
      cout << endl;
    }
  }
  else {
    srep(i, 1, n + 1)
    {
      srep(j, 1, n + 1) { ofs << d[i][j] << ' '; }
      ofs << endl;
    }
  }
}

int Method1(double temperature, double time)
{
  int x = randxor() % n + 1;
  int y = randxor() % n + 1;
  int z = d[x][y];
  int dir = randxor() % 4;
  int x2 = x + dx[dir];
  int y2 = y + dy[dir];
  int z2 = d[x2][y2];
  if (z == z2) {
    return 0;
  }

  if (z == 0) {
    return 0;
  }

  if (time < TL / 2) {
    double g1 = abs(gx[z] / cntM[z] - 25.5) + abs(gy[z] / cntM[z] - 25.5);
    double g2 = abs(gx[z2] / cntM[z2] - 25.5) + abs(gy[z2] / cntM[z2] - 25.5);
    if (g1 > g2) {
      return 0;
    }
  }

  // 連結チェック
  {
    int ok = 1;
    int r = 0;
    rep(i, 4)
    {
      int nx = x + dx[i];
      int ny = y + dy[i];
      if (d[nx][ny] == z) {
        r += (1 << i);
      }
    }
    switch (r) {
      case 0:
        // NG
        return 0;
      case 1:
        // OK
        break;
      case 2:
        // OK
        break;
      case 3:
        // 上左
        if (d[x - 1][y - 1] != z) {
          ok = 0;
        }
        break;
      case 4:
        // OK
        break;
      case 5:
        // 上下 NG
        return 0;
      case 6:
        // 左下
        if (d[x + 1][y - 1] != z) {
          ok = 0;
        }
        break;
      case 7:
        // 上左下
        if (d[x - 1][y - 1] != z || d[x + 1][y - 1] != z) {
          ok = 0;
        }
        break;
      case 8:
        // OK
        break;
      case 9:
        // 上右
        if (d[x - 1][y + 1] != z) {
          ok = 0;
        }
        break;
      case 10:
        // 左右 NG
        return 0;
      case 11:
        // 上左右
        if (d[x - 1][y - 1] != z || d[x - 1][y + 1] != z) {
          ok = 0;
        }
        break;
      case 12:
        // 下右
        if (d[x + 1][y + 1] != z) {
          ok = 0;
        }
        break;
      case 13:
        // 上下右
        if (d[x - 1][y + 1] != z || d[x + 1][y + 1] != z) {
          ok = 0;
        }
        break;
      case 14:
        // 左下右
        if (d[x + 1][y - 1] != z || d[x + 1][y + 1] != z) {
          ok = 0;
        }
        break;
      case 15:
        // あり得ない
        break;
    }
    if (ok == 0) {
      return 0;
    }
  }

  // 隣接チェック
  {
    int ng = 0;
    rep(i, 4)
    {
      int nx = x + dx[i];
      int ny = y + dy[i];
      int nz = d[nx][ny];
      if (nz == z) {
        continue;
      }
      g[z][nz]--;
      g[nz][z]--;
      if (g[z][nz] == 0) {
        ng = 1;
      }
    }
    if (ng) {
      rep(i, 4)
      {
        int nx = x + dx[i];
        int ny = y + dy[i];
        int nz = d[nx][ny];
        if (nz == z) {
          continue;
        }
        g[z][nz]++;
        g[nz][z]++;
      }
      return 0;
    }

    d[x][y] = z2;
    rep(i, 4)
    {
      int nx = x + dx[i];
      int ny = y + dy[i];
      int nz = d[nx][ny];
      if (nz == z2) {
        continue;
      }
      if (g[z2][nz] == 0) {
        ng = 1;
      }
      g[z2][nz]++;
      g[nz][z2]++;
    }

    if (ng) {
      rep(i, 4)
      {
        int nx = x + dx[i];
        int ny = y + dy[i];
        int nz = d[nx][ny];
        if (nz == z2) {
          continue;
        }
        g[z2][nz]--;
        g[nz][z2]--;
      }
      d[x][y] = z;
      return 0;
    }
  }

  // OK
  cntM[z]--;
  gx[z] -= x;
  gy[z] -= y;
  cntM[z2]++;
  gx[z2] += x;
  gy[z2] += y;
  int res = 0;
  if (z2 == 0) {
    zero++;
    res = 1;
  }
  return res;
}

void SA1()
{
  // realに退避

  // 焼きなまし
  clock_t startTime, endTime;
  startTime = clock();
  endTime = clock();
  double startTemperature = 0;
  double endTemperature = 0;
  int loopCount = 0;
  double nowProgress = 0;
  double nowTime;
  while (true) {
    if (loopCount % 100 == 0) {
      endTime = clock();
      nowTime = (double)(endTime - startTime) / CLOCKS_PER_SEC;
      if (nowTime > TL) {
        break;
      }
      nowProgress = nowTime / TL;
    }

    loopCount++;

    double temperature =
      startTemperature + (endTemperature - startTemperature) * nowProgress;

    int tmp = Method1(temperature, nowTime);
    // if (tmp > 0) {
    //   cout << loopCount << endl;
    // }
  }

  // 戻す
}

ll Solve(int probNum)
{
  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  // 入力受け取り
  Input(probNum);

  // 出力ファイルストリームオープン
  ofstream ofs;
  OpenOfs(probNum, ofs);

  // 初期解生成
  Initialize();

  SA1();

  // 解答を出力
  Output(ofs);

  if (ofs.is_open()) {
    ofs.close();
  }

  ll score = 0;
  if (mode != 0) {
    score = CalcScore();
  }
  return score;
}

int main()
{
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    randxor();
  }

  mode = 0;

  if (mode == 0) {
    Solve(0);
  }
  else if (mode == 1) {
    ll sum = 0;
    srep(i, 0, 15)
    {
      ll score = Solve(i);
      sum += score;
      cout << "num = " << i << ", ";
      cout << "score = " << score << ", ";
      cout << "sum = " << sum << endl;
    }
  }

  return 0;
}
