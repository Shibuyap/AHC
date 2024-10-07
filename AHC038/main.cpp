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
#define drep(i, n) for (int i = (n) - 1; i >= 0; --i)
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

  // 配列シャッフル
  void FisherYates(int* data, int n)
  {
    for (int i = n - 1; i >= 0; i--) {
      int j = randxor() % (i + 1);
      int swa = data[i];
      data[i] = data[j];
      data[j] = swa;
    }
  }
}  // namespace

// 配列シャッフル
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
// std::shuffle(v.begin(), v.end(), engine);

const ll INF = 1001001001001001001;
const int INT_INF = 1001001001;

const int dx[4] = { 1, 0, -1, 0 };
const int dy[4] = { 0, -1, 0, 1 };
const char dirChar[5] = { 'D', 'L', 'U', 'R', '.' };
const char rotChar[3] = { 'L', '.', 'R' };
const char tipChar[2] = { '.', 'P' };

const double TL = 2.8;
int mode;
clock_t startTime, endTime;

const int MAX_N = 30;
const int MAX_M = 450;
const int MAX_V = 15;
const int MAX_T = 100005;

int n, m, v;
int init_a[MAX_N][MAX_N];
int init_b[MAX_N][MAX_N];

int V;
int pa[MAX_V];
int le[MAX_V];
int ansCount;
int dir[MAX_T];
int rot[MAX_T][MAX_V];
int tip[MAX_T][MAX_V];
int sx;
int sy;

int real_V;
int real_pa[MAX_V];
int real_le[MAX_V];
int real_ansCount;
int real_dir[MAX_T];
int real_rot[MAX_T][MAX_V];
int real_tip[MAX_T][MAX_V];
int real_sx;
int real_sy;

int route[MAX_N * MAX_N * 3][2];

void CopyToReal()
{
  real_V = V;
  rep(i, V)
  {
    real_pa[i] = pa[i];
    real_le[i] = le[i];
  }
  real_ansCount = ansCount;
  rep(i, ansCount)
  {
    real_dir[i] = dir[i];
    rep(j, V)
    {
      real_rot[i][j] = rot[i][j];
      real_tip[i][j] = tip[i][j];
    }
  }
  real_sx = sx;
  real_sy = sy;
}

void CopyToAns()
{
  V = real_V;
  rep(i, V)
  {
    pa[i] = real_pa[i];
    le[i] = real_le[i];
  }
  ansCount = real_ansCount;
  rep(i, ansCount)
  {
    dir[i] = real_dir[i];
    rep(j, V)
    {
      rot[i][j] = real_rot[i][j];
      tip[i][j] = real_tip[i][j];
    }
  }
  sx = real_sx;
  sy = real_sy;
}

double GetNowTime()
{
  endTime = clock();
  double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
  return nowTime;
}

bool IsNG(int x, int y)
{
  if (x < 0 || n <= x || y < 0 || n <= y)return true;
  return false;
}

// 複数ケース回すときに内部状態を初期値に戻す
void SetUp()
{
  V = 0;
  ansCount = 0;
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
    cin >> n >> m >> v;
    rep(i, n)
    {
      string s;
      cin >> s;
      rep(j, n)
      {
        init_a[i][j] = s[j] - '0';
      }
    }
    rep(i, n)
    {
      string t;
      cin >> t;
      rep(j, n)
      {
        init_b[i][j] = t[j] - '0';
      }
    }
  }
  else {
    // ファイル入力する
    ifs >> n >> m >> v;
    rep(i, n)
    {
      string s;
      ifs >> s;
      rep(j, n)
      {
        init_a[i][j] = s[j] - '0';
      }
    }
    rep(i, n)
    {
      string t;
      ifs >> t;
      rep(j, n)
      {
        init_b[i][j] = t[j] - '0';
      }
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
int CalcScore()
{
  return ansCount;
}

// 初期解生成
void Initialize()
{
  V = v;
  srep(i, 1, V)
  {
    pa[i] = 0;
    le[i] = i;
  }
  sx = 0;
  sy = 0;

  {
    int idx = 0;
    rep(i, n)
    {
      if (i % 2 == 0) {
        rep(j, n)
        {
          route[idx][0] = i;
          route[idx][1] = j;
          idx++;
        }
      }
      else {
        drep(j, n)
        {
          route[idx][0] = i;
          route[idx][1] = j;
          idx++;
        }
      }
    }
    drep(i, n * n)
    {
      route[idx][0] = route[i][0];
      route[idx][1] = route[i][1];
      idx++;
    }
  }
}

// 解答出力
void Output(ofstream& ofs)
{
  if (mode == 0) {
    cout << V << endl;
    srep(i, 1, V)
    {
      cout << pa[i] << ' ' << le[i] << endl;
    }
    cout << sx << ' ' << sy << endl;
    rep(i, ansCount)
    {
      cout << dirChar[dir[i]];
      srep(j, 1, v)cout << rotChar[rot[i][j]];
      rep(j, V)cout << tipChar[tip[i][j]];
      cout << endl;
    }
  }
  else {
    ofs << V << endl;
    srep(i, 1, V)
    {
      ofs << pa[i] << ' ' << le[i] << endl;
    }
    ofs << sx << ' ' << sy << endl;
    rep(i, ansCount)
    {
      ofs << dirChar[dir[i]];
      srep(j, 1, v)ofs << rotChar[rot[i][j]];
      rep(j, V)ofs << tipChar[tip[i][j]];
      ofs << endl;
    }
  }
}

void Method1()
{
  int nn2 = n * n * 2;

  int x = sx;
  int y = sy;
  int t = 0;

  int nowRot[MAX_V];
  int nowTip[MAX_V];
  rep(i, V)
  {
    nowRot[i] = 3;
    nowTip[i] = 0;
  }

  int a[MAX_N][MAX_N];
  int b[MAX_N][MAX_N];
  int mCount = 0;
  rep(i, n)
  {
    rep(j, n)
    {
      a[i][j] = init_a[i][j];
      b[i][j] = init_b[i][j];
      if (a[i][j] == 1 && b[i][j] == 1) {
        mCount++;
        a[i][j] = 0;
        b[i][j] = 0;
      }
    }
  }

  int action[MAX_V] = {};
  while (mCount < m) {
    // dir
    int nx = route[t % nn2][0];
    int ny = route[t % nn2][1];
    dir[t] = 4;
    rep(i, 4)
    {
      if (nx == x + dx[i] && ny == y + dy[i])dir[t] = i;
    }

    // rot, tip
    rep(i, V)
    {
      rot[t][i] = 1;
      tip[t][i] = 0;
      action[i] = 0;
    }

    // このターン
    srep(i, 1, V)
    {
      srep(j, -1, 2)
      {
        int nRot = (nowRot[i] + j + 4) % 4;
        int nrx = nx + le[i] * dx[nRot];
        int nry = ny + le[i] * dy[nRot];

        if (IsNG(nrx, nry))continue;

        if (nowTip[i] == 0) {
          if (a[nrx][nry] == 1) {
            a[nrx][nry] = 0;
            rot[t][i] = j + 1;
            nowRot[i] = nRot;
            tip[t][i] = 1;
            nowTip[i] = 1;
            action[i] = 1;
            break;
          }
        }
        else {
          if (b[nrx][nry] == 1) {
            b[nrx][nry] = 0;
            rot[t][i] = j + 1;
            nowRot[i] = nRot;
            tip[t][i] = 1;
            nowTip[i] = 0;
            action[i] = 1;
            mCount++;
            break;
          }
        }
      }
    }

    // 次のターン
    int nnx = route[(t + 1) % nn2][0];
    int nny = route[(t + 1) % nn2][1];
    srep(i, 1, V)
    {
      if (action[i])continue;

      srep(j, -1, 3)
      {
        int nRot = (nowRot[i] + j + 4) % 4;
        int nnrx = nnx + le[i] * dx[nRot];
        int nnry = nny + le[i] * dy[nRot];

        if (IsNG(nnrx, nnry))continue;

        if (nowTip[i] == 0) {
          if (a[nnrx][nnry] == 1) {
            if (j == 2) {
              rot[t][i] = 2;
              nowRot[i] = (nowRot[i] + 1) % 4;
            }
            else {
              rot[t][i] = j + 1;
              nowRot[i] = nRot;
            }
            action[i] = 1;
            break;
          }
        }
        else {
          if (b[nnrx][nnry] == 1) {
            if (j == 2) {
              rot[t][i] = 2;
              nowRot[i] = (nowRot[i] + 1) % 4;
            }
            else {
              rot[t][i] = j + 1;
              nowRot[i] = nRot;
            }
            action[i] = 1;
            break;
          }
        }
      }
    }

    t++;
    x = nx;
    y = ny;
  }

  ansCount = t;
}

void Method2()
{
  int loop = 0;
  while (true) {

    if (GetNowTime() > TL)break;

    loop++;

    int nn2 = n * n * 2;

    V = v;
    srep(i, 1, V)
    {
      pa[i] = 0;
      le[i] = randxor() % (n - 1) + 1;
    }

    int startT = randxor() % nn2;
    sx = route[startT][0];
    sy = route[startT][1];

    int x = sx;
    int y = sy;
    int _t = 0;

    int nowRot[MAX_V];
    int nowTip[MAX_V];
    rep(i, V)
    {
      nowRot[i] = 3;
      nowTip[i] = 0;
    }

    int a[MAX_N][MAX_N];
    int b[MAX_N][MAX_N];
    int mCount = 0;
    rep(i, n)
    {
      rep(j, n)
      {
        a[i][j] = init_a[i][j];
        b[i][j] = init_b[i][j];
        if (a[i][j] == 1 && b[i][j] == 1) {
          mCount++;
          a[i][j] = 0;
          b[i][j] = 0;
        }
      }
    }

    int action[MAX_V] = {};
    while (mCount < m && _t < real_ansCount) {
      int t = (_t + startT) % nn2;

      // dir
      int nx = route[t][0];
      int ny = route[t][1];
      dir[_t] = 4;
      rep(i, 4)
      {
        if (nx == x + dx[i] && ny == y + dy[i])dir[_t] = i;
      }

      // rot, tip
      rep(i, V)
      {
        rot[_t][i] = 1;
        tip[_t][i] = 0;
        action[i] = 0;
      }

      // このターン
      srep(i, 1, V)
      {
        srep(j, -1, 2)
        {
          int nRot = (nowRot[i] + j + 4) % 4;
          int nrx = nx + le[i] * dx[nRot];
          int nry = ny + le[i] * dy[nRot];

          if (IsNG(nrx, nry))continue;

          if (nowTip[i] == 0) {
            if (a[nrx][nry] == 1) {
              a[nrx][nry] = 0;
              rot[_t][i] = j + 1;
              nowRot[i] = nRot;
              tip[_t][i] = 1;
              nowTip[i] = 1;
              action[i] = 1;
              break;
            }
          }
          else {
            if (b[nrx][nry] == 1) {
              b[nrx][nry] = 0;
              rot[_t][i] = j + 1;
              nowRot[i] = nRot;
              tip[_t][i] = 1;
              nowTip[i] = 0;
              action[i] = 1;
              mCount++;
              break;
            }
          }
        }
      }

      // 次のターン
      int nnx = route[(t + 1) % nn2][0];
      int nny = route[(t + 1) % nn2][1];
      srep(i, 1, V)
      {
        if (action[i])continue;

        srep(j, -1, 3)
        {
          int nRot = (nowRot[i] + j + 4) % 4;
          int nnrx = nnx + le[i] * dx[nRot];
          int nnry = nny + le[i] * dy[nRot];

          if (IsNG(nnrx, nnry))continue;

          if (nowTip[i] == 0) {
            if (a[nnrx][nnry] == 1) {
              if (j == 2) {
                rot[_t][i] = 2;
                nowRot[i] = (nowRot[i] + 1) % 4;
              }
              else {
                rot[_t][i] = j + 1;
                nowRot[i] = nRot;
              }
              action[i] = 1;
              break;
            }
          }
          else {
            if (b[nnrx][nnry] == 1) {
              if (j == 2) {
                rot[_t][i] = 2;
                nowRot[i] = (nowRot[i] + 1) % 4;
              }
              else {
                rot[_t][i] = j + 1;
                nowRot[i] = nRot;
              }
              action[i] = 1;
              break;
            }
          }
        }
      }

      _t++;
      x = nx;
      y = ny;
    }

    ansCount = _t;
    if (mCount == m && ansCount < real_ansCount) {
      CopyToReal();
    }
  }

  if (mode != 0) {
    cout << "Method2 loop = " << loop << endl;
  }
}

void Method3()
{
  int loop = 0;
  while (true) {

    if (GetNowTime() > TL)break;

    loop++;

    int nn2 = n * n * 2;

    V = v;
    srep(i, 1, V)
    {
      pa[i] = 0;
      le[i] = randxor() % (n - 1) + 1;
    }

    int startT = randxor() % nn2;
    sx = route[startT][0];
    sy = route[startT][1];

    int t = startT;

    int x = sx;
    int y = sy;
    int _t = 0;

    int nowRot[MAX_V];
    int nowTip[MAX_V];
    rep(i, V)
    {
      nowRot[i] = 3;
      nowTip[i] = 0;
    }

    int a[MAX_N][MAX_N];
    int b[MAX_N][MAX_N];
    int mCount = 0;
    rep(i, n)
    {
      rep(j, n)
      {
        a[i][j] = init_a[i][j];
        b[i][j] = init_b[i][j];
        if (a[i][j] == 1 && b[i][j] == 1) {
          mCount++;
          a[i][j] = 0;
          b[i][j] = 0;
        }
      }
    }

    int action[MAX_V] = {};
    int lastT = -1;
    int lastX = sx;
    int lastY = sy;
    while (mCount < m && lastT < real_ansCount) {
      // dir
      int nx = route[t][0];
      int ny = route[t][1];
      dir[_t] = 4;
      rep(i, 4)
      {
        if (nx == x + dx[i] && ny == y + dy[i])dir[_t] = i;
      }

      // rot, tip
      rep(i, V)
      {
        rot[_t][i] = 1;
        tip[_t][i] = 0;
        action[i] = 0;
      }

      // このターン
      srep(i, 1, V)
      {
        srep(j, -1, 2)
        {
          int nRot = (nowRot[i] + j + 4) % 4;
          int nrx = nx + le[i] * dx[nRot];
          int nry = ny + le[i] * dy[nRot];

          if (IsNG(nrx, nry))continue;

          if (nowTip[i] == 0) {
            if (a[nrx][nry] == 1) {
              a[nrx][nry] = 0;
              rot[_t][i] = j + 1;
              nowRot[i] = nRot;
              tip[_t][i] = 1;
              nowTip[i] = 1;
              action[i] = 1;
              break;
            }
          }
          else {
            if (b[nrx][nry] == 1) {
              b[nrx][nry] = 0;
              rot[_t][i] = j + 1;
              nowRot[i] = nRot;
              tip[_t][i] = 1;
              nowTip[i] = 0;
              action[i] = 1;
              mCount++;
              break;
            }
          }
        }
      }

      // 次のターン
      int nnx = route[(t + 1) % nn2][0];
      int nny = route[(t + 1) % nn2][1];
      srep(i, 1, V)
      {
        if (action[i])continue;

        srep(j, -1, 3)
        {
          int nRot = (nowRot[i] + j + 4) % 4;
          int nnrx = nnx + le[i] * dx[nRot];
          int nnry = nny + le[i] * dy[nRot];

          if (IsNG(nnrx, nnry))continue;

          if (nowTip[i] == 0) {
            if (a[nnrx][nnry] == 1) {
              if (j == 2) {
                rot[_t][i] = 2;
                nowRot[i] = (nowRot[i] + 1) % 4;
              }
              else {
                rot[_t][i] = j + 1;
                nowRot[i] = nRot;
              }
              action[i] = 1;
              break;
            }
          }
          else {
            if (b[nnrx][nnry] == 1) {
              if (j == 2) {
                rot[_t][i] = 2;
                nowRot[i] = (nowRot[i] + 1) % 4;
              }
              else {
                rot[_t][i] = j + 1;
                nowRot[i] = nRot;
              }
              action[i] = 1;
              break;
            }
          }
        }
      }

      int isAction = 0;
      rep(i, V)
      {
        if (action[i] == 1) {
          isAction = 1;
          break;
        }
      }

      if (isAction) {
        int keepT = _t;
        _t = lastT + 1;
        x = lastX;
        y = lastY;
        while (x != nx || y != ny) {
          rep(i, V)
          {
            rot[_t][i] = 1;
            tip[_t][i] = 0;
            action[i] = 0;
          }
          if (x < nx) {
            dir[_t] = 0;
            x++;
          }
          else if (x > nx) {
            dir[_t] = 2;
            x--;
          }
          else if (y < ny) {
            dir[_t] = 3;
            y++;
          }
          else {
            y--;
          }
          if (x == nx && y == ny) {
            break;
          }
          else {
            _t++;
          }
        }

        srep(i, 1, V)
        {
          rot[_t][i] = rot[keepT][i];
          tip[_t][i] = tip[keepT][i];
        }

        lastT = _t;
        lastX = x;
        lastY = y;
        _t++;
      }
      else {
        _t++;
        x = nx;
        y = ny;
      }
      t = (t + 1) % nn2;
    }

    ansCount = _t;
    if (mCount == m && ansCount < real_ansCount) {
      CopyToReal();
    }
  }

  if (mode != 0) {
    cout << "Method2 loop = " << loop << endl;
  }
}


ll Solve(int probNum)
{
  startTime = clock();

  // 複数ケース回すときに内部状態を初期値に戻す
  SetUp();

  // 入力受け取り
  Input(probNum);

  // 出力ファイルストリームオープン
  ofstream ofs;
  OpenOfs(probNum, ofs);

  // 初期解生成
  Initialize();

  Method1();

  CopyToReal();

  //Method2();
  Method3();

  CopyToAns();

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

  mode = 1;

  if (mode == 0) {
    Solve(0);
  }
  else if (mode == 1) {
    ll sum = 0;
    srep(i, 0, 100)
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
