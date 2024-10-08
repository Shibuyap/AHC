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

const int dx[5] = { 1, 0, -1, 0, 0 };
const int dy[5] = { 0, -1, 0, 1, 0 };
const char dirChar[5] = { 'D', 'L', 'U', 'R', '.' };
const char rotChar[3] = { 'L', '.', 'R' };
const char tipChar[2] = { '.', 'P' };
const int BASE_DIR = 3;

const double TL = 2.8;
int mode;
clock_t startTime, endTime;

const int MAX_N = 30;
const int MAX_M = 450;
const int MAX_V = 15;
const int MAX_T = 100005;

const int ACTION_RATIO = 5;

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
int startT;
int armCount;
int leafs1;
int Method;

int real_V;
int real_pa[MAX_V];
int real_le[MAX_V];
int real_ansCount;
int real_dir[MAX_T];
int real_rot[MAX_T][MAX_V];
int real_tip[MAX_T][MAX_V];
int real_sx;
int real_sy;
int real_startT;
int real_armCount;
int real_leafs1;
int real_Method;

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
  real_startT = startT;
  real_armCount = armCount;
  real_leafs1 = leafs1;
  real_Method = Method;
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
  startT = real_startT;
  armCount = real_armCount;
  leafs1 = real_leafs1;
  Method = real_Method;
}

void ResetTime()
{
  startTime = clock();
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
    //if (V < 10) {
    //  assert(false);
    //}
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

//初期位置(0, 0)
//辺の長さは1~V - 1
//行動できるときは行動する
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
            action[i] = ACTION_RATIO;
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
            action[i] = ACTION_RATIO;
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

//初期位置はランダム
//辺の長さはランダム
//行動できるときは行動する
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
              action[i] = ACTION_RATIO;
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
              action[i] = ACTION_RATIO;
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
    //cout << "Method2 loop = " << loop << endl;
  }
}

//初期位置はランダム
//辺の長さはランダム
//行動できるときは行動する
//不要な移動はショートカットする
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
    while (mCount < m && _t < real_ansCount + 20) {
      // dir
      int nx = route[(t + 1) % nn2][0];
      int ny = route[(t + 1) % nn2][1];
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
              action[i] = ACTION_RATIO;
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
              action[i] = ACTION_RATIO;
              mCount++;
              break;
            }
          }
        }
      }

      // 次のターン
      int nnx = route[(t + 2) % nn2][0];
      int nny = route[(t + 2) % nn2][1];
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
            dir[_t] = 1;
            y--;
          }
          if (x == nx && y == ny) {
            break;
          }
          else {
            _t++;
          }
        }

        srep(i, 0, V)
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
    //cout << "Method3 loop = " << loop << endl;
  }
}

//初期位置はランダム
//関節一つ
//辺の長さはランダム
//行動できるときは行動する
//不要な移動はショートカットする
void Method4(double timeLimit)
{
  ResetTime();
  Method = 4;

  real_startT = -1;
  int loop = 0;
  while (true) {

    double time = GetNowTime();
    if (time > timeLimit)break;

    loop++;

    int nn2 = n * n * 2;

    V = v;
    srep(i, 1, V)
    {
      if (i == 1) {
        pa[i] = 0;
      }
      else {
        pa[i] = 1;
      }
      le[i] = randxor() % (n - 1) + 1;
    }

    startT = randxor() % nn2;
    sx = route[startT][0];
    sy = route[startT][1];
    if (randxor() % 2 == 0 && real_startT != -1 && time > timeLimit * 0.75) {
      V = real_V;
      srep(i, 1, V)
      {
        pa[i] = real_pa[i];
        le[i] = real_le[i];
      }
      startT = real_startT;
      sx = real_sx;
      sy = real_sy;

      int ra = randxor() % 2;
      if (ra == 0) {
        startT = randxor() % nn2;
        sx = route[startT][0];
        sy = route[startT][1];
      }
      else {
        while (true) {
          int raV = randxor() % (V - 2) + 2;
          int newLe = randxor() % (n - 1) + 1;
          if (newLe == le[raV])continue;
          le[raV] = newLe;
          break;
        }
      }
    }

    int t = startT;

    int x = sx;
    int y = sy;
    int _t = 0;

    int nowRot[MAX_V];
    int nowTip[MAX_V];
    rep(i, V)
    {
      nowRot[i] = 0;
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

    int maxActionScore;
    int maxAction;
    int maxRot[MAX_V];
    int maxTip[MAX_V];
    int maxNowRot[MAX_V];
    int maxNowTip[MAX_V];

    int tmpRot[MAX_V];
    int tmpTip[MAX_V];
    int tmpNowRot[MAX_V];
    int tmpNowTip[MAX_V];

    int maxA[MAX_V][3];
    int maxB[MAX_V][3];
    int maxACount = 0;
    int maxBCount = 0;
    int keepA[MAX_V][3];
    int keepB[MAX_V][3];
    int keepACount = 0;
    int keepBCount = 0;

    while (mCount < m && _t < real_ansCount + 20 && lastT < real_ansCount) {
      // dir
      int nx = route[(t + 1) % nn2][0];
      int ny = route[(t + 1) % nn2][1];
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


      maxAction = -1;
      maxActionScore = -1;
      rep(i, V)
      {
        maxRot[i] = 1;
        maxTip[i] = 0;
        maxNowRot[i] = nowRot[i];
        maxNowTip[i] = nowTip[i];
      }

      srep(iii, 0, 3)
      {
        int ii = iii;
        if (ii == 2) ii -= 3;
        srep(jjj, 0, 3)
        {
          int jj = jjj;
          if (jj == 2)jj -= 3;

          int nRot1 = (BASE_DIR + nowRot[1] + ii + 4) % 4;
          int nRot2 = (nRot1 + jj + 4) % 4;

          rep(i, V)
          {
            action[i] = 0;
            tmpRot[i] = 1;
            tmpTip[i] = 0;
            tmpNowRot[i] = nowRot[i];
            tmpNowTip[i] = nowTip[i];
          }

          keepACount = 0;
          keepBCount = 0;

          // このターン
          srep(i, 2, V)
          {
            srep(j, -1, 2)
            {
              int nRot = (nRot1 + nowRot[i] + j + 4) % 4;
              int nrx = nx + le[i] * dx[nRot] + le[1] * dx[nRot1];
              int nry = ny + le[i] * dy[nRot] + le[1] * dy[nRot1];

              if (IsNG(nrx, nry))continue;

              if (nowTip[i] == 0) {
                if (a[nrx][nry] == 1) {
                  keepA[keepACount][0] = nrx;
                  keepA[keepACount][1] = nry;
                  keepA[keepACount][2] = a[nrx][nry];
                  keepACount++;
                  a[nrx][nry] = 0;
                  tmpRot[i] = j + 1;
                  tmpNowRot[i] = (nowRot[i] + j) % 4;
                  tmpTip[i] = 1;
                  tmpNowTip[i] = 1;
                  action[i] = ACTION_RATIO;
                  break;
                }
              }
              else {
                if (b[nrx][nry] == 1) {
                  keepB[keepBCount][0] = nrx;
                  keepB[keepBCount][1] = nry;
                  keepB[keepBCount][2] = b[nrx][nry];
                  keepBCount++;
                  b[nrx][nry] = 0;
                  tmpRot[i] = j + 1;
                  tmpNowRot[i] = (nowRot[i] + j) % 4;
                  tmpTip[i] = 1;
                  tmpNowTip[i] = 0;
                  action[i] = ACTION_RATIO;
                  break;
                }
              }
            }
          }

          // 次のターン
          int nnx = route[(t + 2) % nn2][0];
          int nny = route[(t + 2) % nn2][1];
          srep(i, 2, V)
          {
            if (action[i])continue;

            srep(j, -1, 3)
            {
              int nRot = (nRot2 + nowRot[i] + j + 4) % 4;
              int nnrx = nnx + le[i] * dx[nRot] + le[1] * dx[nRot2];
              int nnry = nny + le[i] * dy[nRot] + le[1] * dy[nRot2];

              if (IsNG(nnrx, nnry))continue;

              if (nowTip[i] == 0) {
                if (a[nnrx][nnry] == 1) {
                  if (j == 2) {
                    tmpRot[i] = 2;
                    tmpNowRot[i] = (nowRot[i] + 1) % 4;
                  }
                  else {
                    tmpRot[i] = j + 1;
                    tmpNowRot[i] = (nowRot[i] + j) % 4;
                  }
                  action[i] = 1;
                  break;
                }
              }
              else {
                if (b[nnrx][nnry] == 1) {
                  if (j == 2) {
                    tmpRot[i] = 2;
                    tmpNowRot[i] = (nowRot[i] + 1) % 4;
                  }
                  else {
                    tmpRot[i] = j + 1;
                    tmpNowRot[i] = (nowRot[i] + j) % 4;
                  }
                  action[i] = 1;
                  break;
                }
              }
            }
          }

          int tmpActionScore = 0;
          srep(i, 2, V)
          {
            tmpActionScore += action[i];
          }

          if (tmpActionScore > maxActionScore) {
            maxAction = ii;
            maxActionScore = tmpActionScore;
            rep(i, V)
            {
              maxRot[i] = tmpRot[i];
              maxTip[i] = tmpTip[i];
              maxNowRot[i] = tmpNowRot[i];
              maxNowTip[i] = tmpNowTip[i];
            }

            maxACount = keepACount;
            rep(i, keepACount)
            {
              rep(j, 3)
              {
                maxA[i][j] = keepA[i][j];
              }
            }
            maxBCount = keepBCount;
            rep(i, keepBCount)
            {
              rep(j, 3)
              {
                maxB[i][j] = keepB[i][j];
              }
            }
          }

          rep(i, keepACount)
          {
            a[keepA[i][0]][keepA[i][1]] = keepA[i][2];
          }
          rep(i, keepBCount)
          {
            b[keepB[i][0]][keepB[i][1]] = keepB[i][2];
          }
        }
      }

      srep(i, 2, V)
      {
        rot[_t][i] = maxRot[i];
        tip[_t][i] = maxTip[i];
        nowRot[i] = maxNowRot[i];
        nowTip[i] = maxNowTip[i];
      }

      rot[_t][1] = maxAction + 1;
      nowRot[1] = (nowRot[1] + maxAction + 4) % 4;

      mCount += maxBCount;
      rep(i, maxACount)
      {
        a[maxA[i][0]][maxA[i][1]] = 0;
      }
      rep(i, maxBCount)
      {
        b[maxB[i][0]][maxB[i][1]] = 0;
      }

      int isAction = 0;
      if (maxActionScore > 0) {
        isAction = 1;
      }

      if (isAction) {
        int keepT = _t;
        _t = lastT + 1;
        x = lastX;
        y = lastY;
        while (x != nx || y != ny) {
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
            dir[_t] = 1;
            y--;
          }
          if (x == nx && y == ny) {
            break;
          }
          else {
            _t++;
          }
        }

        rep(i, V)
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

  if (mode == 2) {
    cout << "Method4 loop = " << loop << " " << endl;
  }
}

//初期位置はランダム
//関節二つ
//辺の長さはランダム
//行動できるときは行動する
//不要な移動はショートカットする
void Method5(double timeLimit)
{
  ResetTime();
  Method = 5;

  real_startT = -1;
  int loop = 0;
  while (true) {

    double time = GetNowTime();
    if (time > timeLimit)break;

    loop++;

    int nn2 = n * n * 2;

    V = v;
    srep(i, 1, V)
    {
      if (i == 1) {
        pa[i] = 0;
      }
      else if (i == 2) {
        pa[i] = 1;
      }
      else {
        pa[i] = 2;
      }
      le[i] = randxor() % (n - 1) + 1;
    }

    startT = randxor() % nn2;
    sx = route[startT][0];
    sy = route[startT][1];
    if (randxor() % 2 == 0 && real_startT != -1 && time > timeLimit * 0.75) {
      V = real_V;
      srep(i, 1, V)
      {
        pa[i] = real_pa[i];
        le[i] = real_le[i];
      }
      startT = real_startT;
      sx = real_sx;
      sy = real_sy;

      int ra = randxor() % 2;
      if (ra == 0) {
        startT = randxor() % nn2;
        sx = route[startT][0];
        sy = route[startT][1];
      }
      else {
        while (true) {
          int raV = randxor() % (V - 2) + 2;
          int newLe = randxor() % (n - 1) + 1;
          if (newLe == le[raV])continue;
          le[raV] = newLe;
          break;
        }
      }
    }

    int t = startT;

    int x = sx;
    int y = sy;
    int _t = 0;

    int nowRot[MAX_V];
    int nowTip[MAX_V];
    rep(i, V)
    {
      nowRot[i] = 0;
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

    int maxActionScore;
    int maxAction;
    int maxAction2;
    int maxRot[MAX_V];
    int maxTip[MAX_V];
    int maxNowRot[MAX_V];
    int maxNowTip[MAX_V];

    int tmpRot[MAX_V];
    int tmpTip[MAX_V];
    int tmpNowRot[MAX_V];
    int tmpNowTip[MAX_V];

    int maxA[MAX_V][3];
    int maxB[MAX_V][3];
    int maxACount = 0;
    int maxBCount = 0;
    int keepA[MAX_V][3];
    int keepB[MAX_V][3];
    int keepACount = 0;
    int keepBCount = 0;

    while (mCount < m && _t < real_ansCount + 20 && lastT < real_ansCount) {
      // dir
      int nx = route[(t + 1) % nn2][0];
      int ny = route[(t + 1) % nn2][1];
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


      maxAction = -1;
      maxAction2 = -1;
      maxActionScore = -1;
      rep(i, V)
      {
        maxRot[i] = 1;
        maxTip[i] = 0;
        maxNowRot[i] = nowRot[i];
        maxNowTip[i] = nowTip[i];
      }

      srep(ii2, 0, 3)
      {
        int ii = ii2;
        if (ii == 2) ii -= 3;
        srep(jj2, 0, 3)
        {
          int jj = jj2;
          if (jj == 2)jj -= 3;

          int nRot1 = (BASE_DIR + nowRot[1] + ii + 4) % 4;
          int nRot2 = (nRot1 + jj + 4) % 4;

          srep(iii2, 0, 3)
          {
            int iii = iii2;
            if (iii == 2) iii -= 3;
            srep(jjj2, 0, 3)
            {
              int jjj = jjj2;
              if (jjj == 2)jjj -= 3;

              int nRot3 = (nRot1 + nowRot[2] + iii + 4) % 4;
              int nRot4 = (nRot2 + nowRot[2] + iii + jjj + 4) % 4;

              rep(i, V)
              {
                action[i] = 0;
                tmpRot[i] = 1;
                tmpTip[i] = 0;
                tmpNowRot[i] = nowRot[i];
                tmpNowTip[i] = nowTip[i];
              }

              keepACount = 0;
              keepBCount = 0;

              // このターン
              srep(i, 3, V)
              {
                srep(j, -1, 2)
                {
                  int nRot = (nRot3 + nowRot[i] + j + 4) % 4;
                  int nrx = nx + le[i] * dx[nRot] + le[1] * dx[nRot1] + le[2] * dx[nRot3];
                  int nry = ny + le[i] * dy[nRot] + le[1] * dy[nRot1] + le[2] * dy[nRot3];

                  if (IsNG(nrx, nry))continue;

                  if (nowTip[i] == 0) {
                    if (a[nrx][nry] == 1) {
                      keepA[keepACount][0] = nrx;
                      keepA[keepACount][1] = nry;
                      keepA[keepACount][2] = a[nrx][nry];
                      keepACount++;
                      a[nrx][nry] = 0;
                      tmpRot[i] = j + 1;
                      tmpNowRot[i] = (nowRot[i] + j) % 4;
                      tmpTip[i] = 1;
                      tmpNowTip[i] = 1;
                      action[i] = ACTION_RATIO;
                      break;
                    }
                  }
                  else {
                    if (b[nrx][nry] == 1) {
                      keepB[keepBCount][0] = nrx;
                      keepB[keepBCount][1] = nry;
                      keepB[keepBCount][2] = b[nrx][nry];
                      keepBCount++;
                      b[nrx][nry] = 0;
                      tmpRot[i] = j + 1;
                      tmpNowRot[i] = (nowRot[i] + j) % 4;
                      tmpTip[i] = 1;
                      tmpNowTip[i] = 0;
                      action[i] = ACTION_RATIO;
                      break;
                    }
                  }
                }
              }

              // 次のターン
              int nnx = route[(t + 2) % nn2][0];
              int nny = route[(t + 2) % nn2][1];
              srep(i, 3, V)
              {
                if (action[i])continue;

                srep(j, -1, 3)
                {
                  int nRot = (nRot4 + nowRot[i] + j + 4) % 4;
                  int nnrx = nnx + le[i] * dx[nRot] + le[1] * dx[nRot2] + le[2] * dx[nRot4];
                  int nnry = nny + le[i] * dy[nRot] + le[1] * dy[nRot2] + le[2] * dy[nRot4];

                  if (IsNG(nnrx, nnry))continue;

                  if (nowTip[i] == 0) {
                    if (a[nnrx][nnry] == 1) {
                      if (j == 2) {
                        tmpRot[i] = 2;
                        tmpNowRot[i] = (nowRot[i] + 1) % 4;
                      }
                      else {
                        tmpRot[i] = j + 1;
                        tmpNowRot[i] = (nowRot[i] + j) % 4;
                      }
                      action[i] = 1;
                      break;
                    }
                  }
                  else {
                    if (b[nnrx][nnry] == 1) {
                      if (j == 2) {
                        tmpRot[i] = 2;
                        tmpNowRot[i] = (nowRot[i] + 1) % 4;
                      }
                      else {
                        tmpRot[i] = j + 1;
                        tmpNowRot[i] = (nowRot[i] + j) % 4;
                      }
                      action[i] = 1;
                      break;
                    }
                  }
                }
              }


              int tmpActionScore = 0;
              srep(i, 3, V)
              {
                tmpActionScore += action[i];
              }

              if (tmpActionScore > maxActionScore) {
                maxAction = ii;
                maxAction2 = iii;
                maxActionScore = tmpActionScore;
                rep(i, V)
                {
                  maxRot[i] = tmpRot[i];
                  maxTip[i] = tmpTip[i];
                  maxNowRot[i] = tmpNowRot[i];
                  maxNowTip[i] = tmpNowTip[i];
                }

                maxACount = keepACount;
                rep(i, keepACount)
                {
                  rep(j, 3)
                  {
                    maxA[i][j] = keepA[i][j];
                  }
                }
                maxBCount = keepBCount;
                rep(i, keepBCount)
                {
                  rep(j, 3)
                  {
                    maxB[i][j] = keepB[i][j];
                  }
                }
              }


              rep(i, keepACount)
              {
                a[keepA[i][0]][keepA[i][1]] = keepA[i][2];
              }
              rep(i, keepBCount)
              {
                b[keepB[i][0]][keepB[i][1]] = keepB[i][2];
              }
            }
          }
        }
      }

      srep(i, 3, V)
      {
        rot[_t][i] = maxRot[i];
        tip[_t][i] = maxTip[i];
        nowRot[i] = maxNowRot[i];
        nowTip[i] = maxNowTip[i];
      }

      rot[_t][1] = maxAction + 1;
      nowRot[1] = (nowRot[1] + maxAction + 4) % 4;
      rot[_t][2] = maxAction2 + 1;
      nowRot[2] = (nowRot[2] + maxAction2 + 4) % 4;

      mCount += maxBCount;
      rep(i, maxACount)
      {
        a[maxA[i][0]][maxA[i][1]] = 0;
      }
      rep(i, maxBCount)
      {
        b[maxB[i][0]][maxB[i][1]] = 0;
      }

      int isAction = 0;
      if (maxActionScore > 0) {
        isAction = 1;
      }

      if (isAction) {
        int keepT = _t;
        _t = lastT + 1;
        x = lastX;
        y = lastY;
        while (x != nx || y != ny) {
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
            dir[_t] = 1;
            y--;
          }
          if (x == nx && y == ny) {
            break;
          }
          else {
            _t++;
          }
        }

        rep(i, V)
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
      //cout << real_ansCount << ' ' << le[1] << endl;
    }
  }

  if (mode == 2) {
    cout << "Method5 loop = " << loop << " " << endl;
  }
}

//初期位置はランダム
//関節二つ
//辺の長さはランダム
//行動できるときは行動する
//不要な移動はショートカットする
// ランダムウォーク
void Method52(double timeLimit)
{
  ResetTime();
  Method = 52;

  real_startT = -1;
  int loop = 0;
  while (true) {

    double time = GetNowTime();
    if (time > timeLimit)break;

    loop++;

    int nn2 = n * n * 2;

    V = v;
    srep(i, 1, V)
    {
      if (i == 1) {
        pa[i] = 0;
      }
      else if (i == 2) {
        pa[i] = 1;
      }
      else {
        pa[i] = 2;
      }
      le[i] = randxor() % (n - 1) + 1;
    }

    startT = randxor() % nn2;
    sx = route[startT][0];
    sy = route[startT][1];
    if (randxor() % 2 == 0 && real_startT != -1 && time > timeLimit * 0.75) {
      V = real_V;
      srep(i, 1, V)
      {
        pa[i] = real_pa[i];
        le[i] = real_le[i];
      }
      startT = real_startT;
      sx = real_sx;
      sy = real_sy;

      int ra = randxor() % 2;
      if (ra == 0) {
        startT = randxor() % nn2;
        sx = route[startT][0];
        sy = route[startT][1];
      }
      else {
        while (true) {
          int raV = randxor() % (V - 2) + 2;
          int newLe = randxor() % (n - 1) + 1;
          if (newLe == le[raV])continue;
          le[raV] = newLe;
          break;
        }
      }
    }

    int t = startT;

    int x = sx;
    int y = sy;
    int _t = 0;

    int nowRot[MAX_V];
    int nowTip[MAX_V];
    rep(i, V)
    {
      nowRot[i] = 0;
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

    int maxDir;

    int maxActionScore;
    int maxAction;
    int maxAction2;
    int maxRot[MAX_V];
    int maxTip[MAX_V];
    int maxNowRot[MAX_V];
    int maxNowTip[MAX_V];

    int tmpRot[MAX_V];
    int tmpTip[MAX_V];
    int tmpNowRot[MAX_V];
    int tmpNowTip[MAX_V];

    int maxA[MAX_V][3];
    int maxB[MAX_V][3];
    int maxACount = 0;
    int maxBCount = 0;
    int keepA[MAX_V][3];
    int keepB[MAX_V][3];
    int keepACount = 0;
    int keepBCount = 0;

    int order[5] = { 0,1,2,3,4 };

    while (mCount < m && _t < real_ansCount + 20 && lastT < real_ansCount) {
      FisherYates(order, 5);

      maxDir = -1;
      maxAction = -1;
      maxAction2 = -1;
      maxActionScore = -1;
      rep(i, V)
      {
        maxRot[i] = 1;
        maxTip[i] = 0;
        maxNowRot[i] = nowRot[i];
        maxNowTip[i] = nowTip[i];
      }

      rep(ord, 5)
      {
        dir[_t] = order[ord];
        int nx = x + dx[order[ord]];
        int ny = y + dy[order[ord]];

        // rot, tip
        rep(i, V)
        {
          rot[_t][i] = 1;
          tip[_t][i] = 0;
          action[i] = 0;
        }

        srep(ii2, 0, 3)
        {
          int ii = ii2;
          if (ii == 2) ii -= 3;
          srep(jj2, 0, 3)
          {
            int jj = jj2;
            if (jj == 2)jj -= 3;

            int nRot1 = (BASE_DIR + nowRot[1] + ii + 4) % 4;
            int nRot2 = (nRot1 + jj + 4) % 4;

            srep(iii2, 0, 3)
            {
              int iii = iii2;
              if (iii == 2) iii -= 3;
              srep(jjj2, 0, 3)
              {
                int jjj = jjj2;
                if (jjj == 2)jjj -= 3;

                int nRot3 = (nRot1 + nowRot[2] + iii + 4) % 4;
                int nRot4 = (nRot2 + nowRot[2] + iii + jjj + 4) % 4;

                rep(i, V)
                {
                  action[i] = 0;
                  tmpRot[i] = 1;
                  tmpTip[i] = 0;
                  tmpNowRot[i] = nowRot[i];
                  tmpNowTip[i] = nowTip[i];
                }

                keepACount = 0;
                keepBCount = 0;

                // このターン
                srep(i, 3, V)
                {
                  srep(j, -1, 2)
                  {
                    int nRot = (nRot3 + nowRot[i] + j + 4) % 4;
                    int nrx = nx + le[i] * dx[nRot] + le[1] * dx[nRot1] + le[2] * dx[nRot3];
                    int nry = ny + le[i] * dy[nRot] + le[1] * dy[nRot1] + le[2] * dy[nRot3];

                    if (IsNG(nrx, nry))continue;

                    if (nowTip[i] == 0) {
                      if (a[nrx][nry] == 1) {
                        keepA[keepACount][0] = nrx;
                        keepA[keepACount][1] = nry;
                        keepA[keepACount][2] = a[nrx][nry];
                        keepACount++;
                        a[nrx][nry] = 0;
                        tmpRot[i] = j + 1;
                        tmpNowRot[i] = (nowRot[i] + j) % 4;
                        tmpTip[i] = 1;
                        tmpNowTip[i] = 1;
                        action[i] = ACTION_RATIO;
                        break;
                      }
                    }
                    else {
                      if (b[nrx][nry] == 1) {
                        keepB[keepBCount][0] = nrx;
                        keepB[keepBCount][1] = nry;
                        keepB[keepBCount][2] = b[nrx][nry];
                        keepBCount++;
                        b[nrx][nry] = 0;
                        tmpRot[i] = j + 1;
                        tmpNowRot[i] = (nowRot[i] + j) % 4;
                        tmpTip[i] = 1;
                        tmpNowTip[i] = 0;
                        action[i] = ACTION_RATIO;
                        break;
                      }
                    }
                  }
                }

                // 次のターン
                //int nnx = route[(t + 2) % nn2][0];
                //int nny = route[(t + 2) % nn2][1];
                //srep(i, 3, V)
                //{
                //  if (action[i])continue;

                //  srep(j, -1, 3)
                //  {
                //    int nRot = (nRot4 + nowRot[i] + j + 4) % 4;
                //    int nnrx = nnx + le[i] * dx[nRot] + le[1] * dx[nRot2] + le[2] * dx[nRot4];
                //    int nnry = nny + le[i] * dy[nRot] + le[1] * dy[nRot2] + le[2] * dy[nRot4];

                //    if (IsNG(nnrx, nnry))continue;

                //    if (nowTip[i] == 0) {
                //      if (a[nnrx][nnry] == 1) {
                //        if (j == 2) {
                //          tmpRot[i] = 2;
                //          tmpNowRot[i] = (nowRot[i] + 1) % 4;
                //        }
                //        else {
                //          tmpRot[i] = j + 1;
                //          tmpNowRot[i] = (nowRot[i] + j) % 4;
                //        }
                //        action[i] = 1;
                //        break;
                //      }
                //    }
                //    else {
                //      if (b[nnrx][nnry] == 1) {
                //        if (j == 2) {
                //          tmpRot[i] = 2;
                //          tmpNowRot[i] = (nowRot[i] + 1) % 4;
                //        }
                //        else {
                //          tmpRot[i] = j + 1;
                //          tmpNowRot[i] = (nowRot[i] + j) % 4;
                //        }
                //        action[i] = 1;
                //        break;
                //      }
                //    }
                //  }
                //}


                int tmpActionScore = 0;
                srep(i, 3, V)
                {
                  tmpActionScore += action[i];
                }

                if (tmpActionScore > maxActionScore) {
                  maxDir = order[ord];
                  maxAction = ii;
                  maxAction2 = iii;
                  maxActionScore = tmpActionScore;
                  rep(i, V)
                  {
                    maxRot[i] = tmpRot[i];
                    maxTip[i] = tmpTip[i];
                    maxNowRot[i] = tmpNowRot[i];
                    maxNowTip[i] = tmpNowTip[i];
                  }

                  maxACount = keepACount;
                  rep(i, keepACount)
                  {
                    rep(j, 3)
                    {
                      maxA[i][j] = keepA[i][j];
                    }
                  }
                  maxBCount = keepBCount;
                  rep(i, keepBCount)
                  {
                    rep(j, 3)
                    {
                      maxB[i][j] = keepB[i][j];
                    }
                  }
                }


                rep(i, keepACount)
                {
                  a[keepA[i][0]][keepA[i][1]] = keepA[i][2];
                }
                rep(i, keepBCount)
                {
                  b[keepB[i][0]][keepB[i][1]] = keepB[i][2];
                }
              }
            }
          }
        }
      }


      srep(i, 3, V)
      {
        rot[_t][i] = maxRot[i];
        tip[_t][i] = maxTip[i];
        nowRot[i] = maxNowRot[i];
        nowTip[i] = maxNowTip[i];
      }

      dir[_t] = maxDir;
      rot[_t][1] = maxAction + 1;
      nowRot[1] = (nowRot[1] + maxAction + 4) % 4;
      rot[_t][2] = maxAction2 + 1;
      nowRot[2] = (nowRot[2] + maxAction2 + 4) % 4;

      mCount += maxBCount;
      rep(i, maxACount)
      {
        a[maxA[i][0]][maxA[i][1]] = 0;
      }
      rep(i, maxBCount)
      {
        b[maxB[i][0]][maxB[i][1]] = 0;
      }

      int isAction = 0;
      if (maxActionScore > 0) {
        isAction = 1;
      }

      int nx = x + dx[maxDir];
      int ny = y + dy[maxDir];

      if (isAction) {
        int keepT = _t;
        _t = lastT + 1;
        x = lastX;
        y = lastY;
        while (x != nx || y != ny) {
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
            dir[_t] = 1;
            y--;
          }
          if (x == nx && y == ny) {
            break;
          }
          else {
            _t++;
          }
        }

        rep(i, V)
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
    }

    ansCount = _t;
    if (mCount == m && ansCount < real_ansCount) {
      CopyToReal();
      //cout << real_ansCount << ' ' << le[1] << endl;
    }
  }

  if (mode == 2) {
    cout << "Method52 loop = " << loop << " " << endl;
  }
}

//初期位置はランダム
//関節三つ
//辺の長さはランダム
//行動できるときは行動する
//不要な移動はショートカットする
void Method6(double timeLimit)
{
  ResetTime();
  Method = 6;

  real_startT = -1;
  int loop = 0;
  while (true) {

    double time = GetNowTime();
    if (time > timeLimit)break;

    loop++;

    int nn2 = n * n * 2;

    V = v;
    srep(i, 1, V)
    {
      if (i == 1) {
        pa[i] = 0;
      }
      else if (i == 2) {
        pa[i] = 1;
      }
      else if (i == 3) {
        pa[i] = 2;
      }
      else {
        pa[i] = 3;
      }
      le[i] = randxor() % (n - 1) + 1;
    }

    startT = randxor() % nn2;
    sx = route[startT][0];
    sy = route[startT][1];
    if (randxor() % 2 == 0 && real_startT != -1 && time > timeLimit * 0.75) {
      V = real_V;
      srep(i, 1, V)
      {
        pa[i] = real_pa[i];
        le[i] = real_le[i];
      }
      startT = real_startT;
      sx = real_sx;
      sy = real_sy;

      int ra = randxor() % 2;
      if (ra == 0) {
        startT = randxor() % nn2;
        sx = route[startT][0];
        sy = route[startT][1];
      }
      else {
        while (true) {
          int raV = randxor() % (V - 2) + 2;
          int newLe = randxor() % (n - 1) + 1;
          if (newLe == le[raV])continue;
          le[raV] = newLe;
          break;
        }
      }
    }

    int t = startT;

    int x = sx;
    int y = sy;
    int _t = 0;

    int nowRot[MAX_V];
    int nowTip[MAX_V];
    rep(i, V)
    {
      nowRot[i] = 0;
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

    int maxActionScore;
    int maxAction;
    int maxAction2;
    int maxAction3;
    int maxRot[MAX_V];
    int maxTip[MAX_V];
    int maxNowRot[MAX_V];
    int maxNowTip[MAX_V];

    int tmpRot[MAX_V];
    int tmpTip[MAX_V];
    int tmpNowRot[MAX_V];
    int tmpNowTip[MAX_V];

    int maxA[MAX_V][3];
    int maxB[MAX_V][3];
    int maxACount = 0;
    int maxBCount = 0;
    int keepA[MAX_V][3];
    int keepB[MAX_V][3];
    int keepACount = 0;
    int keepBCount = 0;

    while (mCount < m && _t < real_ansCount + 20 && lastT < real_ansCount) {
      // dir
      int nx = route[(t + 1) % nn2][0];
      int ny = route[(t + 1) % nn2][1];
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


      maxAction = -1;
      maxAction2 = -1;
      maxAction3 = -1;
      maxActionScore = -1;
      rep(i, V)
      {
        maxRot[i] = 1;
        maxTip[i] = 0;
        maxNowRot[i] = nowRot[i];
        maxNowTip[i] = nowTip[i];
      }

      srep(ii2, 0, 3)
      {
        int ii = ii2;
        if (ii == 2) ii -= 3;
        srep(jj2, 0, 3)
        {
          int jj = jj2;
          if (jj == 2)jj -= 3;

          int nRot1 = (BASE_DIR + nowRot[1] + ii + 4) % 4;
          int nRot2 = (BASE_DIR + nowRot[1] + ii + jj + 4) % 4;

          srep(iii2, 0, 3)
          {
            int iii = iii2;
            if (iii == 2) iii -= 3;
            srep(jjj2, 0, 3)
            {
              int jjj = jjj2;
              if (jjj == 2)jjj -= 3;

              int nRot3 = (nRot1 + nowRot[2] + iii + 4) % 4;
              int nRot4 = (nRot2 + nowRot[2] + iii + jjj + 4) % 4;

              srep(iiii2, 0, 3)
              {
                int iiii = iiii2;
                if (iiii == 2) iiii -= 3;
                srep(jjjj2, 0, 3)
                {
                  int jjjj = jjjj2;
                  if (jjjj == 2)jjjj -= 3;

                  int nRot5 = (nRot3 + nowRot[3] + iiii + 4) % 4;
                  int nRot6 = (nRot4 + nowRot[3] + iiii + jjjj + 4) % 4;


                  rep(i, V)
                  {
                    action[i] = 0;
                    tmpRot[i] = 1;
                    tmpTip[i] = 0;
                    tmpNowRot[i] = nowRot[i];
                    tmpNowTip[i] = nowTip[i];
                  }

                  keepACount = 0;
                  keepBCount = 0;

                  // このターン
                  srep(i, 4, V)
                  {
                    srep(j, -1, 2)
                    {
                      int nRot = (nRot5 + nowRot[i] + j + 4) % 4;
                      int nrx = nx + le[i] * dx[nRot] + le[1] * dx[nRot1] + le[2] * dx[nRot3] + le[3] * dx[nRot5];
                      int nry = ny + le[i] * dy[nRot] + le[1] * dy[nRot1] + le[2] * dy[nRot3] + le[3] * dy[nRot5];

                      if (IsNG(nrx, nry))continue;

                      if (nowTip[i] == 0) {
                        if (a[nrx][nry] == 1) {
                          keepA[keepACount][0] = nrx;
                          keepA[keepACount][1] = nry;
                          keepA[keepACount][2] = a[nrx][nry];
                          keepACount++;
                          a[nrx][nry] = 0;
                          tmpRot[i] = j + 1;
                          tmpNowRot[i] = (nowRot[i] + j) % 4;
                          tmpTip[i] = 1;
                          tmpNowTip[i] = 1;
                          action[i] = ACTION_RATIO;
                          break;
                        }
                      }
                      else {
                        if (b[nrx][nry] == 1) {
                          keepB[keepBCount][0] = nrx;
                          keepB[keepBCount][1] = nry;
                          keepB[keepBCount][2] = b[nrx][nry];
                          keepBCount++;
                          b[nrx][nry] = 0;
                          tmpRot[i] = j + 1;
                          tmpNowRot[i] = (nowRot[i] + j) % 4;
                          tmpTip[i] = 1;
                          tmpNowTip[i] = 0;
                          action[i] = ACTION_RATIO;
                          break;
                        }
                      }
                    }
                  }

                  // 次のターン
                  int nnx = route[(t + 2) % nn2][0];
                  int nny = route[(t + 2) % nn2][1];
                  srep(i, 4, V)
                  {
                    if (action[i])continue;

                    srep(j, -1, 3)
                    {
                      int nRot = (nRot6 + nowRot[i] + j + 4) % 4;
                      int nnrx = nnx + le[i] * dx[nRot] + le[1] * dx[nRot2] + le[2] * dx[nRot4] + le[3] * dx[nRot6];
                      int nnry = nny + le[i] * dy[nRot] + le[1] * dy[nRot2] + le[2] * dy[nRot4] + le[3] * dy[nRot6];

                      if (IsNG(nnrx, nnry))continue;

                      if (nowTip[i] == 0) {
                        if (a[nnrx][nnry] == 1) {
                          if (j == 2) {
                            tmpRot[i] = 2;
                            tmpNowRot[i] = (nowRot[i] + 1) % 4;
                          }
                          else {
                            tmpRot[i] = j + 1;
                            tmpNowRot[i] = (nowRot[i] + j) % 4;
                          }
                          action[i] = 1;
                          break;
                        }
                      }
                      else {
                        if (b[nnrx][nnry] == 1) {
                          if (j == 2) {
                            tmpRot[i] = 2;
                            tmpNowRot[i] = (nowRot[i] + 1) % 4;
                          }
                          else {
                            tmpRot[i] = j + 1;
                            tmpNowRot[i] = (nowRot[i] + j) % 4;
                          }
                          action[i] = 1;
                          break;
                        }
                      }
                    }
                  }

                  int tmpActionScore = 0;
                  srep(i, 4, V)
                  {
                    tmpActionScore += action[i];
                  }

                  if (tmpActionScore > maxActionScore) {
                    maxAction = ii;
                    maxAction2 = iii;
                    maxAction3 = iiii;
                    maxActionScore = tmpActionScore;
                    rep(i, V)
                    {
                      maxRot[i] = tmpRot[i];
                      maxTip[i] = tmpTip[i];
                      maxNowRot[i] = tmpNowRot[i];
                      maxNowTip[i] = tmpNowTip[i];
                    }

                    maxACount = keepACount;
                    rep(i, keepACount)
                    {
                      rep(j, 3)
                      {
                        maxA[i][j] = keepA[i][j];
                      }
                    }
                    maxBCount = keepBCount;
                    rep(i, keepBCount)
                    {
                      rep(j, 3)
                      {
                        maxB[i][j] = keepB[i][j];
                      }
                    }
                  }


                  rep(i, keepACount)
                  {
                    a[keepA[i][0]][keepA[i][1]] = keepA[i][2];
                  }
                  rep(i, keepBCount)
                  {
                    b[keepB[i][0]][keepB[i][1]] = keepB[i][2];
                  }

                }
              }
            }
          }
        }
      }

      srep(i, 4, V)
      {
        rot[_t][i] = maxRot[i];
        tip[_t][i] = maxTip[i];
        nowRot[i] = maxNowRot[i];
        nowTip[i] = maxNowTip[i];
      }

      rot[_t][1] = maxAction + 1;
      nowRot[1] = (nowRot[1] + maxAction + 4) % 4;
      rot[_t][2] = maxAction2 + 1;
      nowRot[2] = (nowRot[2] + maxAction2 + 4) % 4;
      rot[_t][3] = maxAction3 + 1;
      nowRot[3] = (nowRot[3] + maxAction3 + 4) % 4;

      mCount += maxBCount;
      rep(i, maxACount)
      {
        a[maxA[i][0]][maxA[i][1]] = 0;
      }
      rep(i, maxBCount)
      {
        b[maxB[i][0]][maxB[i][1]] = 0;
      }

      int isAction = 0;
      if (maxActionScore > 0) {
        isAction = 1;
      }

      if (isAction) {
        int keepT = _t;
        _t = lastT + 1;
        x = lastX;
        y = lastY;
        while (x != nx || y != ny) {
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
            dir[_t] = 1;
            y--;
          }
          if (x == nx && y == ny) {
            break;
          }
          else {
            _t++;
          }
        }

        rep(i, V)
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
      //cout << real_ansCount << ' ' << le[1] << endl;
    }
  }

  if (mode == 2) {
    cout << "Method6 loop = " << loop << " " << endl;
  }
}

//初期位置はランダム
//関節三つ
//辺の長さはランダム
//行動できるときは行動する
//不要な移動はショートカットする
// ランダムウォーク
void Method62(double timeLimit)
{
  ResetTime();
  Method = 62;

  real_startT = -1;
  int loop = 0;
  while (true) {

    double time = GetNowTime();
    if (time > timeLimit)break;

    loop++;

    int nn2 = n * n * 2;

    V = v;
    srep(i, 1, V)
    {
      if (i == 1) {
        pa[i] = 0;
      }
      else if (i == 2) {
        pa[i] = 1;
      }
      else if (i == 3) {
        pa[i] = 2;
      }
      else {
        pa[i] = 3;
      }
      le[i] = randxor() % (n - 1) + 1;
    }

    startT = randxor() % nn2;
    sx = route[startT][0];
    sy = route[startT][1];
    if (randxor() % 2 == 0 && real_startT != -1 && time > timeLimit * 0.75) {
      V = real_V;
      srep(i, 1, V)
      {
        pa[i] = real_pa[i];
        le[i] = real_le[i];
      }
      startT = real_startT;
      sx = real_sx;
      sy = real_sy;

      int ra = randxor() % 2;
      if (ra == 0) {
        startT = randxor() % nn2;
        sx = route[startT][0];
        sy = route[startT][1];
      }
      else {
        while (true) {
          int raV = randxor() % (V - 2) + 2;
          int newLe = randxor() % (n - 1) + 1;
          if (newLe == le[raV])continue;
          le[raV] = newLe;
          break;
        }
      }
    }

    int t = startT;

    int x = sx;
    int y = sy;
    int _t = 0;

    int nowRot[MAX_V];
    int nowTip[MAX_V];
    rep(i, V)
    {
      nowRot[i] = 0;
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

    int maxDir;

    int maxActionScore;
    int maxAction;
    int maxAction2;
    int maxAction3;
    int maxRot[MAX_V];
    int maxTip[MAX_V];
    int maxNowRot[MAX_V];
    int maxNowTip[MAX_V];

    int tmpRot[MAX_V];
    int tmpTip[MAX_V];
    int tmpNowRot[MAX_V];
    int tmpNowTip[MAX_V];

    int maxA[MAX_V][3];
    int maxB[MAX_V][3];
    int maxACount = 0;
    int maxBCount = 0;
    int keepA[MAX_V][3];
    int keepB[MAX_V][3];
    int keepACount = 0;
    int keepBCount = 0;

    int order[5] = { 0,1,2,3,4 };

    while (mCount < m && _t < real_ansCount + 20 && lastT < real_ansCount) {
      FisherYates(order, 5);

      maxDir = -1;
      maxAction = -1;
      maxAction2 = -1;
      maxAction3 = -1;
      maxActionScore = -1;
      rep(i, V)
      {
        maxRot[i] = 1;
        maxTip[i] = 0;
        maxNowRot[i] = nowRot[i];
        maxNowTip[i] = nowTip[i];
      }

      rep(ord, 5)
      {
        // dir
        dir[_t] = order[ord];
        int nx = x + dx[order[ord]];
        int ny = y + dy[order[ord]];

        // rot, tip
        rep(i, V)
        {
          rot[_t][i] = 1;
          tip[_t][i] = 0;
          action[i] = 0;
        }

        srep(ii2, 0, 3)
        {
          int ii = ii2;
          if (ii == 2) ii -= 3;
          srep(jj2, 0, 3)
          {
            int jj = jj2;
            if (jj == 2)jj -= 3;

            int nRot1 = (BASE_DIR + nowRot[1] + ii + 4) % 4;
            int nRot2 = (BASE_DIR + nowRot[1] + ii + jj + 4) % 4;

            srep(iii2, 0, 3)
            {
              int iii = iii2;
              if (iii == 2) iii -= 3;
              srep(jjj2, 0, 3)
              {
                int jjj = jjj2;
                if (jjj == 2)jjj -= 3;

                int nRot3 = (nRot1 + nowRot[2] + iii + 4) % 4;
                int nRot4 = (nRot2 + nowRot[2] + iii + jjj + 4) % 4;

                srep(iiii2, 0, 3)
                {
                  int iiii = iiii2;
                  if (iiii == 2) iiii -= 3;
                  srep(jjjj2, 0, 3)
                  {
                    int jjjj = jjjj2;
                    if (jjjj == 2)jjjj -= 3;

                    int nRot5 = (nRot3 + nowRot[3] + iiii + 4) % 4;
                    int nRot6 = (nRot4 + nowRot[3] + iiii + jjjj + 4) % 4;


                    rep(i, V)
                    {
                      action[i] = 0;
                      tmpRot[i] = 1;
                      tmpTip[i] = 0;
                      tmpNowRot[i] = nowRot[i];
                      tmpNowTip[i] = nowTip[i];
                    }

                    keepACount = 0;
                    keepBCount = 0;

                    // このターン
                    srep(i, 4, V)
                    {
                      srep(j, -1, 2)
                      {
                        int nRot = (nRot5 + nowRot[i] + j + 4) % 4;
                        int nrx = nx + le[i] * dx[nRot] + le[1] * dx[nRot1] + le[2] * dx[nRot3] + le[3] * dx[nRot5];
                        int nry = ny + le[i] * dy[nRot] + le[1] * dy[nRot1] + le[2] * dy[nRot3] + le[3] * dy[nRot5];

                        if (IsNG(nrx, nry))continue;

                        if (nowTip[i] == 0) {
                          if (a[nrx][nry] == 1) {
                            keepA[keepACount][0] = nrx;
                            keepA[keepACount][1] = nry;
                            keepA[keepACount][2] = a[nrx][nry];
                            keepACount++;
                            a[nrx][nry] = 0;
                            tmpRot[i] = j + 1;
                            tmpNowRot[i] = (nowRot[i] + j) % 4;
                            tmpTip[i] = 1;
                            tmpNowTip[i] = 1;
                            action[i] = ACTION_RATIO;
                            break;
                          }
                        }
                        else {
                          if (b[nrx][nry] == 1) {
                            keepB[keepBCount][0] = nrx;
                            keepB[keepBCount][1] = nry;
                            keepB[keepBCount][2] = b[nrx][nry];
                            keepBCount++;
                            b[nrx][nry] = 0;
                            tmpRot[i] = j + 1;
                            tmpNowRot[i] = (nowRot[i] + j) % 4;
                            tmpTip[i] = 1;
                            tmpNowTip[i] = 0;
                            action[i] = ACTION_RATIO;
                            break;
                          }
                        }
                      }
                    }

                    // 次のターン
                    //int nnx = route[(t + 2) % nn2][0];
                    //int nny = route[(t + 2) % nn2][1];
                    //srep(i, 4, V)
                    //{
                    //  if (action[i])continue;

                    //  srep(j, -1, 3)
                    //  {
                    //    int nRot = (nRot6 + nowRot[i] + j + 4) % 4;
                    //    int nnrx = nnx + le[i] * dx[nRot] + le[1] * dx[nRot2] + le[2] * dx[nRot4] + le[3] * dx[nRot6];
                    //    int nnry = nny + le[i] * dy[nRot] + le[1] * dy[nRot2] + le[2] * dy[nRot4] + le[3] * dy[nRot6];

                    //    if (IsNG(nnrx, nnry))continue;

                    //    if (nowTip[i] == 0) {
                    //      if (a[nnrx][nnry] == 1) {
                    //        if (j == 2) {
                    //          tmpRot[i] = 2;
                    //          tmpNowRot[i] = (nowRot[i] + 1) % 4;
                    //        }
                    //        else {
                    //          tmpRot[i] = j + 1;
                    //          tmpNowRot[i] = (nowRot[i] + j) % 4;
                    //        }
                    //        action[i] = 1;
                    //        break;
                    //      }
                    //    }
                    //    else {
                    //      if (b[nnrx][nnry] == 1) {
                    //        if (j == 2) {
                    //          tmpRot[i] = 2;
                    //          tmpNowRot[i] = (nowRot[i] + 1) % 4;
                    //        }
                    //        else {
                    //          tmpRot[i] = j + 1;
                    //          tmpNowRot[i] = (nowRot[i] + j) % 4;
                    //        }
                    //        action[i] = 1;
                    //        break;
                    //      }
                    //    }
                    //  }
                    //}

                    int tmpActionScore = 0;
                    srep(i, 4, V)
                    {
                      tmpActionScore += action[i];
                    }

                    if (tmpActionScore > maxActionScore) {
                      maxDir = order[ord];
                      maxAction = ii;
                      maxAction2 = iii;
                      maxAction3 = iiii;
                      maxActionScore = tmpActionScore;
                      rep(i, V)
                      {
                        maxRot[i] = tmpRot[i];
                        maxTip[i] = tmpTip[i];
                        maxNowRot[i] = tmpNowRot[i];
                        maxNowTip[i] = tmpNowTip[i];
                      }

                      maxACount = keepACount;
                      rep(i, keepACount)
                      {
                        rep(j, 3)
                        {
                          maxA[i][j] = keepA[i][j];
                        }
                      }
                      maxBCount = keepBCount;
                      rep(i, keepBCount)
                      {
                        rep(j, 3)
                        {
                          maxB[i][j] = keepB[i][j];
                        }
                      }
                    }


                    rep(i, keepACount)
                    {
                      a[keepA[i][0]][keepA[i][1]] = keepA[i][2];
                    }
                    rep(i, keepBCount)
                    {
                      b[keepB[i][0]][keepB[i][1]] = keepB[i][2];
                    }

                  }
                }
              }
            }
          }
        }
      }



      srep(i, 4, V)
      {
        rot[_t][i] = maxRot[i];
        tip[_t][i] = maxTip[i];
        nowRot[i] = maxNowRot[i];
        nowTip[i] = maxNowTip[i];
      }

      dir[_t] = maxDir;
      rot[_t][1] = maxAction + 1;
      nowRot[1] = (nowRot[1] + maxAction + 4) % 4;
      rot[_t][2] = maxAction2 + 1;
      nowRot[2] = (nowRot[2] + maxAction2 + 4) % 4;
      rot[_t][3] = maxAction3 + 1;
      nowRot[3] = (nowRot[3] + maxAction3 + 4) % 4;

      mCount += maxBCount;
      rep(i, maxACount)
      {
        a[maxA[i][0]][maxA[i][1]] = 0;
      }
      rep(i, maxBCount)
      {
        b[maxB[i][0]][maxB[i][1]] = 0;
      }

      int isAction = 0;
      if (maxActionScore > 0) {
        isAction = 1;
      }

      int nx = x + dx[maxDir];
      int ny = y + dy[maxDir];

      if (isAction) {
        int keepT = _t;
        _t = lastT + 1;
        x = lastX;
        y = lastY;
        while (x != nx || y != ny) {
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
            dir[_t] = 1;
            y--;
          }
          if (x == nx && y == ny) {
            break;
          }
          else {
            _t++;
          }
        }

        rep(i, V)
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
    }

    ansCount = _t;
    if (mCount == m && ansCount < real_ansCount) {
      CopyToReal();
      //cout << real_ansCount << ' ' << le[1] << endl;
    }
  }

  if (mode == 2) {
    cout << "Method62 loop = " << loop << " " << endl;
  }
}

//初期位置はランダム
//関節二つ
//腕二本
//辺の長さはランダム
//行動できるときは行動する
//不要な移動はショートカットする
void Method7(double timeLimit)
{
  ResetTime();
  Method = 7;

  armCount = 1;
  real_armCount = 1;

  real_startT = -1;
  int loop = 0;
  while (true) {

    double time = GetNowTime();
    if (time > timeLimit)break;

    loop++;

    int nn2 = n * n * 2;

    V = v;

    armCount = 2;
    leafs1 = randxor() % (V - 6) + 1;

    srep(i, 1, V)
    {
      if (i == 1) {
        pa[i] = 0;
      }
      else if (i == 2) {
        pa[i] = 1;
      }
      else if (i == 3) {
        pa[i] = 0;
      }
      else if (i == 4) {
        pa[i] = 3;
      }
      else if (i < 5 + leafs1) {
        pa[i] = 2;
      }
      else {
        pa[i] = 4;
      }
      le[i] = randxor() % (n - 1) + 1;
    }

    startT = randxor() % nn2;
    sx = route[startT][0];
    sy = route[startT][1];
    if (randxor() % 2 == 0 && real_startT != -1 && time > timeLimit * 0.75) {
      V = real_V;
      srep(i, 1, V)
      {
        pa[i] = real_pa[i];
        le[i] = real_le[i];
      }
      leafs1 = real_leafs1;
      startT = real_startT;
      sx = real_sx;
      sy = real_sy;

      int ra = randxor() % 2;
      if (ra == 0) {
        startT = randxor() % nn2;
        sx = route[startT][0];
        sy = route[startT][1];
      }
      else {
        while (true) {
          int raV = randxor() % (V - 2) + 2;
          int newLe = randxor() % (n - 1) + 1;
          if (newLe == le[raV])continue;
          le[raV] = newLe;
          break;
        }
      }
    }

    int t = startT;

    int x = sx;
    int y = sy;
    int _t = 0;

    int nowRot[MAX_V];
    int nowTip[MAX_V];
    rep(i, V)
    {
      nowRot[i] = 0;
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

    int maxActionScore;
    int maxAction;
    int maxAction2;
    int maxRot[MAX_V];
    int maxTip[MAX_V];
    int maxNowRot[MAX_V];
    int maxNowTip[MAX_V];

    int tmpRot[MAX_V];
    int tmpTip[MAX_V];
    int tmpNowRot[MAX_V];
    int tmpNowTip[MAX_V];

    int maxA[MAX_V][3];
    int maxB[MAX_V][3];
    int maxACount = 0;
    int maxBCount = 0;
    int keepA[MAX_V][3];
    int keepB[MAX_V][3];
    int keepACount = 0;
    int keepBCount = 0;

    while (mCount < m && _t < real_ansCount + 20 && lastT < real_ansCount) {
      // dir
      int nx = route[(t + 1) % nn2][0];
      int ny = route[(t + 1) % nn2][1];
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


      maxAction = -1;
      maxAction2 = -1;
      maxActionScore = -1;
      rep(i, V)
      {
        maxRot[i] = 1;
        maxTip[i] = 0;
        maxNowRot[i] = nowRot[i];
        maxNowTip[i] = nowTip[i];
      }

      srep(ii2, 0, 3)
      {
        int ii = ii2;
        if (ii == 2) ii -= 3;
        srep(jj2, 0, 3)
        {
          int jj = jj2;
          if (jj == 2)jj -= 3;

          int nRot1 = (BASE_DIR + nowRot[1] + ii + 4) % 4;
          int nRot2 = (nRot1 + jj + 4) % 4;

          srep(iii2, 0, 3)
          {
            int iii = iii2;
            if (iii == 2) iii -= 3;
            srep(jjj2, 0, 3)
            {
              int jjj = jjj2;
              if (jjj == 2)jjj -= 3;

              int nRot3 = (nRot1 + nowRot[2] + iii + 4) % 4;
              int nRot4 = (nRot2 + nowRot[2] + iii + jjj + 4) % 4;

              rep(i, V)
              {
                action[i] = 0;
                tmpRot[i] = 1;
                tmpTip[i] = 0;
                tmpNowRot[i] = nowRot[i];
                tmpNowTip[i] = nowTip[i];
              }

              keepACount = 0;
              keepBCount = 0;

              // このターン
              srep(i, 5, 5 + leafs1)
              {
                srep(j, -1, 2)
                {
                  int nRot = (nRot3 + nowRot[i] + j + 4) % 4;
                  int nrx = nx + le[i] * dx[nRot] + le[1] * dx[nRot1] + le[2] * dx[nRot3];
                  int nry = ny + le[i] * dy[nRot] + le[1] * dy[nRot1] + le[2] * dy[nRot3];

                  if (IsNG(nrx, nry))continue;

                  if (nowTip[i] == 0) {
                    if (a[nrx][nry] == 1) {
                      keepA[keepACount][0] = nrx;
                      keepA[keepACount][1] = nry;
                      keepA[keepACount][2] = a[nrx][nry];
                      keepACount++;
                      a[nrx][nry] = 0;
                      tmpRot[i] = j + 1;
                      tmpNowRot[i] = (nowRot[i] + j) % 4;
                      tmpTip[i] = 1;
                      tmpNowTip[i] = 1;
                      action[i] = ACTION_RATIO;
                      break;
                    }
                  }
                  else {
                    if (b[nrx][nry] == 1) {
                      keepB[keepBCount][0] = nrx;
                      keepB[keepBCount][1] = nry;
                      keepB[keepBCount][2] = b[nrx][nry];
                      keepBCount++;
                      b[nrx][nry] = 0;
                      tmpRot[i] = j + 1;
                      tmpNowRot[i] = (nowRot[i] + j) % 4;
                      tmpTip[i] = 1;
                      tmpNowTip[i] = 0;
                      action[i] = ACTION_RATIO;
                      break;
                    }
                  }
                }
              }

              // 次のターン
              int nnx = route[(t + 2) % nn2][0];
              int nny = route[(t + 2) % nn2][1];
              srep(i, 5, 5 + leafs1)
              {
                if (action[i])continue;

                srep(j, -1, 3)
                {
                  int nRot = (nRot4 + nowRot[i] + j + 4) % 4;
                  int nnrx = nnx + le[i] * dx[nRot] + le[1] * dx[nRot2] + le[2] * dx[nRot4];
                  int nnry = nny + le[i] * dy[nRot] + le[1] * dy[nRot2] + le[2] * dy[nRot4];

                  if (IsNG(nnrx, nnry))continue;

                  if (nowTip[i] == 0) {
                    if (a[nnrx][nnry] == 1) {
                      if (j == 2) {
                        tmpRot[i] = 2;
                        tmpNowRot[i] = (nowRot[i] + 1) % 4;
                      }
                      else {
                        tmpRot[i] = j + 1;
                        tmpNowRot[i] = (nowRot[i] + j) % 4;
                      }
                      action[i] = 1;
                      break;
                    }
                  }
                  else {
                    if (b[nnrx][nnry] == 1) {
                      if (j == 2) {
                        tmpRot[i] = 2;
                        tmpNowRot[i] = (nowRot[i] + 1) % 4;
                      }
                      else {
                        tmpRot[i] = j + 1;
                        tmpNowRot[i] = (nowRot[i] + j) % 4;
                      }
                      action[i] = 1;
                      break;
                    }
                  }
                }
              }


              int tmpActionScore = 0;
              srep(i, 3, V)
              {
                tmpActionScore += action[i];
              }

              if (tmpActionScore > maxActionScore) {
                maxAction = ii;
                maxAction2 = iii;
                maxActionScore = tmpActionScore;
                rep(i, V)
                {
                  maxRot[i] = tmpRot[i];
                  maxTip[i] = tmpTip[i];
                  maxNowRot[i] = tmpNowRot[i];
                  maxNowTip[i] = tmpNowTip[i];
                }

                maxACount = keepACount;
                rep(i, keepACount)
                {
                  rep(j, 3)
                  {
                    maxA[i][j] = keepA[i][j];
                  }
                }
                maxBCount = keepBCount;
                rep(i, keepBCount)
                {
                  rep(j, 3)
                  {
                    maxB[i][j] = keepB[i][j];
                  }
                }
              }


              rep(i, keepACount)
              {
                a[keepA[i][0]][keepA[i][1]] = keepA[i][2];
              }
              rep(i, keepBCount)
              {
                b[keepB[i][0]][keepB[i][1]] = keepB[i][2];
              }
            }
          }
        }
      }

      srep(i, 5, 5 + leafs1)
      {
        rot[_t][i] = maxRot[i];
        tip[_t][i] = maxTip[i];
        nowRot[i] = maxNowRot[i];
        nowTip[i] = maxNowTip[i];
      }

      rot[_t][1] = maxAction + 1;
      nowRot[1] = (nowRot[1] + maxAction + 4) % 4;
      rot[_t][2] = maxAction2 + 1;
      nowRot[2] = (nowRot[2] + maxAction2 + 4) % 4;

      mCount += maxBCount;
      rep(i, maxACount)
      {
        a[maxA[i][0]][maxA[i][1]] = 0;
      }
      rep(i, maxBCount)
      {
        b[maxB[i][0]][maxB[i][1]] = 0;
      }

      int isAction = 0;
      if (maxActionScore > 0) {
        isAction = 1;
      }

      // 2本目の腕
      maxAction = -1;
      maxAction2 = -1;
      maxActionScore = -1;
      rep(i, V)
      {
        maxRot[i] = 1;
        maxTip[i] = 0;
        maxNowRot[i] = nowRot[i];
        maxNowTip[i] = nowTip[i];
      }

      srep(ii2, 0, 3)
      {
        int ii = ii2;
        if (ii == 2) ii -= 3;
        srep(jj2, 0, 3)
        {
          int jj = jj2;
          if (jj == 2)jj -= 3;

          int nRot1 = (BASE_DIR + nowRot[3] + ii + 4) % 4;
          int nRot2 = (nRot1 + jj + 4) % 4;

          srep(iii2, 0, 3)
          {
            int iii = iii2;
            if (iii == 2) iii -= 3;
            srep(jjj2, 0, 3)
            {
              int jjj = jjj2;
              if (jjj == 2)jjj -= 3;

              int nRot3 = (nRot1 + nowRot[4] + iii + 4) % 4;
              int nRot4 = (nRot2 + nowRot[4] + iii + jjj + 4) % 4;

              rep(i, V)
              {
                action[i] = 0;
                tmpRot[i] = 1;
                tmpTip[i] = 0;
                tmpNowRot[i] = nowRot[i];
                tmpNowTip[i] = nowTip[i];
              }

              keepACount = 0;
              keepBCount = 0;

              // このターン
              srep(i, 5 + leafs1, V)
              {
                srep(j, -1, 2)
                {
                  int nRot = (nRot3 + nowRot[i] + j + 4) % 4;
                  int nrx = nx + le[i] * dx[nRot] + le[3] * dx[nRot1] + le[4] * dx[nRot3];
                  int nry = ny + le[i] * dy[nRot] + le[3] * dy[nRot1] + le[4] * dy[nRot3];

                  if (IsNG(nrx, nry))continue;

                  if (nowTip[i] == 0) {
                    if (a[nrx][nry] == 1) {
                      keepA[keepACount][0] = nrx;
                      keepA[keepACount][1] = nry;
                      keepA[keepACount][2] = a[nrx][nry];
                      keepACount++;
                      a[nrx][nry] = 0;
                      tmpRot[i] = j + 1;
                      tmpNowRot[i] = (nowRot[i] + j) % 4;
                      tmpTip[i] = 1;
                      tmpNowTip[i] = 1;
                      action[i] = ACTION_RATIO;
                      break;
                    }
                  }
                  else {
                    if (b[nrx][nry] == 1) {
                      keepB[keepBCount][0] = nrx;
                      keepB[keepBCount][1] = nry;
                      keepB[keepBCount][2] = b[nrx][nry];
                      keepBCount++;
                      b[nrx][nry] = 0;
                      tmpRot[i] = j + 1;
                      tmpNowRot[i] = (nowRot[i] + j) % 4;
                      tmpTip[i] = 1;
                      tmpNowTip[i] = 0;
                      action[i] = ACTION_RATIO;
                      break;
                    }
                  }
                }
              }

              // 次のターン
              int nnx = route[(t + 2) % nn2][0];
              int nny = route[(t + 2) % nn2][1];
              srep(i, 5 + leafs1, V)
              {
                if (action[i])continue;

                srep(j, -1, 3)
                {
                  int nRot = (nRot4 + nowRot[i] + j + 4) % 4;
                  int nnrx = nnx + le[i] * dx[nRot] + le[3] * dx[nRot2] + le[4] * dx[nRot4];
                  int nnry = nny + le[i] * dy[nRot] + le[3] * dy[nRot2] + le[4] * dy[nRot4];

                  if (IsNG(nnrx, nnry))continue;

                  if (nowTip[i] == 0) {
                    if (a[nnrx][nnry] == 1) {
                      if (j == 2) {
                        tmpRot[i] = 2;
                        tmpNowRot[i] = (nowRot[i] + 1) % 4;
                      }
                      else {
                        tmpRot[i] = j + 1;
                        tmpNowRot[i] = (nowRot[i] + j) % 4;
                      }
                      action[i] = 1;
                      break;
                    }
                  }
                  else {
                    if (b[nnrx][nnry] == 1) {
                      if (j == 2) {
                        tmpRot[i] = 2;
                        tmpNowRot[i] = (nowRot[i] + 1) % 4;
                      }
                      else {
                        tmpRot[i] = j + 1;
                        tmpNowRot[i] = (nowRot[i] + j) % 4;
                      }
                      action[i] = 1;
                      break;
                    }
                  }
                }
              }


              int tmpActionScore = 0;
              srep(i, 5 + leafs1, V)
              {
                tmpActionScore += action[i];
              }

              if (tmpActionScore > maxActionScore) {
                maxAction = ii;
                maxAction2 = iii;
                maxActionScore = tmpActionScore;
                rep(i, V)
                {
                  maxRot[i] = tmpRot[i];
                  maxTip[i] = tmpTip[i];
                  maxNowRot[i] = tmpNowRot[i];
                  maxNowTip[i] = tmpNowTip[i];
                }

                maxACount = keepACount;
                rep(i, keepACount)
                {
                  rep(j, 3)
                  {
                    maxA[i][j] = keepA[i][j];
                  }
                }
                maxBCount = keepBCount;
                rep(i, keepBCount)
                {
                  rep(j, 3)
                  {
                    maxB[i][j] = keepB[i][j];
                  }
                }
              }


              rep(i, keepACount)
              {
                a[keepA[i][0]][keepA[i][1]] = keepA[i][2];
              }
              rep(i, keepBCount)
              {
                b[keepB[i][0]][keepB[i][1]] = keepB[i][2];
              }
            }
          }
        }
      }

      srep(i, 5 + leafs1, V)
      {
        rot[_t][i] = maxRot[i];
        tip[_t][i] = maxTip[i];
        nowRot[i] = maxNowRot[i];
        nowTip[i] = maxNowTip[i];
      }

      rot[_t][3] = maxAction + 1;
      nowRot[3] = (nowRot[3] + maxAction + 4) % 4;
      rot[_t][4] = maxAction2 + 1;
      nowRot[4] = (nowRot[4] + maxAction2 + 4) % 4;

      mCount += maxBCount;
      rep(i, maxACount)
      {
        a[maxA[i][0]][maxA[i][1]] = 0;
      }
      rep(i, maxBCount)
      {
        b[maxB[i][0]][maxB[i][1]] = 0;
      }

      if (maxActionScore > 0) {
        isAction = 1;
      }

      if (isAction) {
        int keepT = _t;
        _t = lastT + 1;
        x = lastX;
        y = lastY;
        while (x != nx || y != ny) {
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
            dir[_t] = 1;
            y--;
          }
          if (x == nx && y == ny) {
            break;
          }
          else {
            _t++;
          }
        }

        rep(i, V)
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
      //cout << real_ansCount << ' ' << le[1] << endl;
    }
  }

  if (mode == 2) {
    cout << "Method7 loop = " << loop << " " << endl;
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
  //Method3();

  if (v < 7) {
    Method4(TL * 0.1);
    Method52(TL * 0.45);
    Method62(TL * 0.45);
  }
  else {
    Method4(TL * 0.1);
    Method52(TL * 0.3);
    Method62(TL * 0.3);
    Method7(TL * 0.3);
  }

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

/////////////////////////////////////////////////////////////////////////
/*
 TODO

 ・腕の長さをランダムにしない

 ・枝刈り高速化

 ・ビームサーチ化

 ・近傍探索

 ・アドホックな手法考える

リファクタ

 ・関数の共通化
   ・ひとまず共通化する

*/
/////////////////////////////////////////////////////////////////////////
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
  else if (mode <= 2) {
    ll sum = 0;
    srep(i, 0, 100)
    {
      ll score = Solve(i);
      sum += score;
      if (mode == 1) {
        cout << score << endl;
      }
      else {
        cout << "num = " << setw(2) << i << ", ";
        cout << "score = " << setw(4) << score << ", ";
        cout << "sum = " << setw(5) << sum << ", ";
        cout << "N = " << setw(2) << n << ", ";
        cout << "M = " << setw(3) << m << ", ";
        cout << "V = " << setw(2) << v << ", ";
        cout << "Method = " << Method << ", ";
        cout << endl;
      }
    }
  }

  return 0;
}
