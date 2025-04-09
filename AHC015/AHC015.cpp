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

  // 0以上1未満の小数をとる乱数
  static double rand01() { return (Rand() + 0.5) * (1.0 / UINT_MAX); }
}  // namespace

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };
const char cc[4] = { 'F', 'L', 'B', 'R' };

int n = 100;
int f[100];
int iteLocal[100];
char outLocal[100];
int countAll[4];

int a[10][10];
int aa[10][10];
int keepA[10][10];
const double TL = 1.8;

inline void CopyAtoAA()
{
  rep(i, 10)
  {
    rep(j, 10) { aa[i][j] = a[i][j]; }
  }
}
inline void CopyAAtoA()
{
  rep(i, 10)
  {
    rep(j, 10) { a[i][j] = aa[i][j]; }
  }
}
inline void CopyAtoKeepA()
{
  rep(i, 10)
  {
    rep(j, 10) { keepA[i][j] = a[i][j]; }
  }
}
inline void CopyKeepAtoA()
{
  rep(i, 10)
  {
    rep(j, 10) { a[i][j] = keepA[i][j]; }
  }
}

inline int RandomPoint(int turn) { return Rand() % (100 - turn) + 1; }

// 入力受け取り（実行中一度しか呼ばれないことを想定）
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

  rep(i, 100) ifs >> f[i];
  rep(i, 100) ifs >> iteLocal[i];
}

// 解答出力
void Output(int problemNum)
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

  rep(i, 100) ofs << outLocal[i] << endl;

  ofs.close();
}

inline bool IsNG(int x, int y)
{
  if (x < 0 || x >= 10 || y < 0 || y >= 10) return true;
  return false;
}

P Where(int ite)
{
  int cnt = 1;
  rep(i, 10)
  {
    rep(j, 10)
    {
      if (a[i][j] == 0) {
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

int calcVisit[10][10];
//queue<P> calcQue;
int qX[1000], qY[1000];
int CalcAA()
{  // スコア計算
  int qI = 0;
  int qSz = 0;
  double res = 0;
  rep(i, 10) rep(j, 10) calcVisit[i][j] = 0;
  rep(i, 10)
  {
    rep(j, 10)
    {
      if (aa[i][j] == 0) continue;
      if (calcVisit[i][j]) continue;
      //calcQue.push(P(i, j));
      qX[qSz] = i;
      qY[qSz] = j;
      qSz++;
      calcVisit[i][j] = 1;
      int sz = 1;
      while (qI < qSz) {
        //P p = calcQue.front();
        //calcQue.pop();
        int x = qX[qI];
        int y = qY[qI];
        qI++;
        rep(k, 4)
        {
          int nx = x + dx[k];
          int ny = y + dy[k];
          if (IsNG(nx, ny)) continue;
          if (aa[nx][ny] == aa[i][j] && calcVisit[nx][ny] == 0) {
            calcVisit[nx][ny] = 1;
            //calcQue.push(P(nx, ny));
            qX[qSz] = nx;
            qY[qSz] = ny;
            qSz++;
            sz++;
          }
        }
      }

      res += sz * sz;
    }
  }

  res *= 1000000;
  int bo = 0;
  srep(i, 1, 4) { bo += countAll[i] * countAll[i]; }
  res /= bo;

  return round(res);
}


int CalcA()
{  // スコア計算
  int qI = 0;
  int qSz = 0;
  double res = 0;
  rep(i, 10) rep(j, 10) calcVisit[i][j] = 0;
  rep(i, 10)
  {
    rep(j, 10)
    {
      if (a[i][j] == 0) continue;
      if (calcVisit[i][j]) continue;
      //calcQue.push(P(i, j));
      qX[qSz] = i;
      qY[qSz] = j;
      qSz++;
      calcVisit[i][j] = 1;
      int sz = 1;
      while (qI < qSz) {
        //P p = calcQue.front();
        //calcQue.pop();
        int x = qX[qI];
        int y = qY[qI];
        qI++;
        rep(k, 4)
        {
          int nx = x + dx[k];
          int ny = y + dy[k];
          if (IsNG(nx, ny)) continue;
          if (a[nx][ny] == a[i][j] && calcVisit[nx][ny] == 0) {
            calcVisit[nx][ny] = 1;
            //calcQue.push(P(nx, ny));
            qX[qSz] = nx;
            qY[qSz] = ny;
            qSz++;
            sz++;
          }
        }
      }

      res += sz * sz;
    }
  }

  res *= 1000000;
  int bo = 0;
  srep(i, 1, 4) { bo += countAll[i] * countAll[i]; }
  res /= bo;

  return round(res);
}

void MoveAA(int dir)
{
  // 場所移動
  if (dir == 0) {
    rep(j, 10)
    {
      int now = 0;
      rep(i, 10)
      {
        if (aa[i][j] != 0) {
          int tmp = aa[i][j];
          aa[i][j] = 0;
          aa[now][j] = tmp;
          now++;
        }
      }
    }
  }
  else if (dir == 1) {
    rep(i, 10)
    {
      int now = 0;
      rep(j, 10)
      {
        if (aa[i][j] != 0) {
          int tmp = aa[i][j];
          aa[i][j] = 0;
          aa[i][now] = tmp;
          now++;
        }
      }
    }
  }
  else if (dir == 2) {
    rep(j, 10)
    {
      int now = 9;
      drep(i, 10)
      {
        if (aa[i][j] != 0) {
          int tmp = aa[i][j];
          aa[i][j] = 0;
          aa[now][j] = tmp;
          now--;
        }
      }
    }
  }
  else if (dir == 3) {
    rep(i, 10)
    {
      int now = 9;
      drep(j, 10)
      {
        if (aa[i][j] != 0) {
          int tmp = aa[i][j];
          aa[i][j] = 0;
          aa[i][now] = tmp;
          now--;
        }
      }
    }
  }
}

void MoveA(int dir)
{
  // 場所移動
  if (dir == 0) {
    rep(j, 10)
    {
      int now = 0;
      rep(i, 10)
      {
        if (a[i][j] != 0) {
          int tmp = a[i][j];
          a[i][j] = 0;
          a[now][j] = tmp;
          now++;
        }
      }
    }
  }
  else if (dir == 1) {
    rep(i, 10)
    {
      int now = 0;
      rep(j, 10)
      {
        if (a[i][j] != 0) {
          int tmp = a[i][j];
          a[i][j] = 0;
          a[i][now] = tmp;
          now++;
        }
      }
    }
  }
  else if (dir == 2) {
    rep(j, 10)
    {
      int now = 9;
      drep(i, 10)
      {
        if (a[i][j] != 0) {
          int tmp = a[i][j];
          a[i][j] = 0;
          a[now][j] = tmp;
          now--;
        }
      }
    }
  }
  else if (dir == 3) {
    rep(i, 10)
    {
      int now = 9;
      drep(j, 10)
      {
        if (a[i][j] != 0) {
          int tmp = a[i][j];
          a[i][j] = 0;
          a[i][now] = tmp;
          now--;
        }
      }
    }
  }
}

int MoveAndCalcAA(int dir)
{
  MoveAA(dir);
  return CalcAA();
}

int Sim(int startTurn, int simBias)
{
  int simTurn = 10 + simBias;
  srep(_, startTurn, startTurn + simTurn)
  {
    if (_ == 100) break;
    // 場所決め
    int po = RandomPoint(_);
    P p = Where(po);
    int wx = p.first;
    int wy = p.second;
    if (a[wx][wy] != 0) {
      cout << "NG" << endl;
    }
    a[wx][wy] = f[_];

    int tmpMax = -1;
    int tmpAns = -1;
    rep(i, 4)
    {
      CopyAtoAA();
      int tmp = MoveAndCalcAA(i);
      if (tmp > tmpMax) {
        tmpMax = tmp;
        tmpAns = i;
      }
    }

    MoveA(tmpAns);
  }

  return CalcA();
}

int Solve(int mode)
{
  ofstream ofs2("debug.txt");
  clock_t startTime, endTime;
  startTime = clock();
  endTime = clock();

  if (mode == 0) {
    rep(i, 100) { cin >> f[i]; }
  }
  else {
    Input(mode - 1);
  }

  rep(i, 100) { countAll[f[i]]++; }

  rep(_, 100)
  {
    int ite;
    if (mode == 0) {
      cin >> ite;
    }
    else {
      ite = iteLocal[_];
    }

    P p = Where(ite);
    int wx = p.first;
    int wy = p.second;
    if (a[wx][wy] != 0) {
      cout << "NG" << endl;
    }
    a[wx][wy] = f[_];

    CopyAtoKeepA();

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
        simBias = Rand() % 3 - 1;
      }

      CopyKeepAtoA();

      int dir = loop % 4;
      MoveA(dir);
      // ave[dir] = max(ave[dir], (double)Sim(_ + 1));
      ave[dir] = (ave[dir] * cnt[dir] + Sim(_ + 1, simBias)) / (cnt[dir] + 1);
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
      std::cout << cc[tmpAns] << endl;
      std::fflush(stdout);
    }
    else {
      outLocal[_] = cc[tmpAns];
    }

    CopyKeepAtoA();
    MoveA(tmpAns);

    if (mode != 0) {
      ofs2 << std::right << std::setw(7) << loop << ' ' << CalcA() << endl;
    }
  }

  if (mode != 0) {
    Output(mode - 1);
  }

  if (mode != 0) {
    cout << "Score = " << CalcA() << endl;
    // rep(i, 10) {
    //   rep(j, 10) { cout << a[i][j]; }
    //   cout << endl;
    // }
  }

  return 0;
}

int main()
{
  // 乱数調整
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    Rand();
  }

  int mode = 0;
  Solve(mode);

  return 0;
}
