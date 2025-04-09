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

  // 配列シャッフル
  void FisherYates(int* data, int n)
  {
    for (int i = n - 1; i >= 0; i--) {
      int j = Rand() % (i + 1);
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

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };
const char dc[4] = { 'U', 'L', 'D', 'R' };

double TL = 1.8;
int mode;
clock_t startTime, endTime;

const int n = 20;
int baseH[20][20];
int h[20][20];
int ans[100000][2];
int ansSize;

int real_ans[100000][2];
int real_ansSize;
int real_ansScore;

double GetNowTime()
{
  endTime = clock();
  double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
  return nowTime;
}

void InitH()
{
  rep(i, n)
  {
    rep(j, n) { h[i][j] = baseH[i][j]; }
  }
}

void CopyToAns()
{
  ansSize = real_ansSize;
  rep(i, ansSize)
  {
    rep(j, 2) { ans[i][j] = real_ans[i][j]; }
  }
}

void CopyToReal()
{
  real_ansSize = ansSize;
  rep(i, ansSize)
  {
    rep(j, 2) { real_ans[i][j] = ans[i][j]; }
  }
}

// 複数ケース回すときに内部状態を初期値に戻す
void SetUp() { ansSize = 0; }

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
    int nn;
    cin >> nn;
    rep(i, n)
    {
      rep(j, n) { cin >> h[i][j]; }
    }
  }
  // ファイル入力する
  else {
    int nn;
    ifs >> nn;
    rep(i, n)
    {
      rep(j, n) { ifs >> h[i][j]; }
    }
  }

  rep(i, n)
  {
    rep(j, n) { baseH[i][j] = h[i][j]; }
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

int CalcCost()
{
  int d = 0;
  int sum = 0;
  rep(i, ansSize)
  {
    if (ans[i][0] < 4) {
      sum += 100 + d;
    }
    else {
      sum += abs(ans[i][1]);
      d += ans[i][1];
    }
  }
  return sum;
}

// スコア計算
ll CalcScore()
{
  ll base = 0;
  rep(i, n)
  {
    rep(j, n) { base += abs(baseH[i][j]); }
  }
  ll res = round(1e9 * base / CalcCost());
  return res;
}

bool IsOK(int x)
{
  rep(j, n)
  {
    if (h[x][j] != 0) {
      return false;
    }
  }
  return true;
}

// 初期解生成
void Initialize()
{
  int x = 0, y = 0;
  int d = 0;
  int dir = 1;
  while (true) {
    // cout << h[x][y] << endl;
    if (h[x][y] < 0 && d > 0) {
      if (-h[x][y] <= d) {
        ans[ansSize][0] = 4;
        ans[ansSize][1] = h[x][y];
        ansSize++;
        d += h[x][y];
        h[x][y] = 0;
      }
      else {
        ans[ansSize][0] = 4;
        ans[ansSize][1] = -d;
        ansSize++;
        h[x][y] += d;
        d = 0;
      }
    }
    else if (h[x][y] > 0) {
      ans[ansSize][0] = 4;
      ans[ansSize][1] = h[x][y];
      ansSize++;
      d += h[x][y];
      h[x][y] = 0;
    }

    if (dir == 1) {
      if (y == n - 1) {
        if (x == n - 1) break;
        ans[ansSize][0] = 2;
        ansSize++;
        x++;
        dir *= -1;
      }
      else {
        y++;
        ans[ansSize][0] = 3;
        ansSize++;
      }
    }
    else {
      if (y == 0) {
        if (x == n - 1) break;
        x++;
        ans[ansSize][0] = 2;
        ansSize++;
        dir *= -1;
      }
      else {
        y--;
        ans[ansSize][0] = 1;
        ansSize++;
      }
    }
  }

  dir = 1;
  if (y == n - 1) dir = -1;

  while (x >= 0) {
    if (IsOK(x)) {
      if (x == 0) break;
      ans[ansSize][0] = 0;
      ansSize++;
      x--;
      continue;
    }
    if (h[x][y] < 0 && d > 0) {
      if (-h[x][y] <= d) {
        ans[ansSize][0] = 4;
        ans[ansSize][1] = h[x][y];
        ansSize++;
        d += h[x][y];
        h[x][y] = 0;
      }
      else {
        ans[ansSize][0] = 4;
        ans[ansSize][1] = -d;
        h[x][y] += d;
        d = 0;
      }
    }
    else if (h[x][y] > 0) {
      ans[ansSize][0] = 4;
      ans[ansSize][1] = h[x][y];
      ansSize++;
      d += h[x][y];
      h[x][y] = 0;
    }

    if (dir == 1) {
      if (y == n - 1) {
        if (!IsOK(x)) {
          dir *= -1;
          continue;
        }
        if (x == 0) break;
        x--;
        ans[ansSize][0] = 0;
        ansSize++;
        dir *= -1;
      }
      else {
        y++;
        ans[ansSize][0] = 3;
        ansSize++;
      }
    }
    else {
      if (y == 0) {
        if (!IsOK(x)) {
          dir *= -1;
          continue;
        }
        if (x == 0) break;
        x--;
        ans[ansSize][0] = 0;
        ansSize++;
        dir *= -1;
      }
      else {
        y--;
        ans[ansSize][0] = 1;
        ansSize++;
      }
    }
  }

  CopyToReal();
  real_ansScore = CalcScore();
}

void GetLoad(int x, int y, int diff, int& d)
{
  ans[ansSize][0] = 4;
  ans[ansSize][1] = diff;
  ansSize++;
  d += diff;
  h[x][y] -= diff;
}

void PutLoad(int x, int y, int diff, int& d)
{
  ans[ansSize][0] = 4;
  ans[ansSize][1] = -diff;
  ansSize++;
  d -= diff;
  h[x][y] += diff;
}

void Up()
{
  ans[ansSize][0] = 0;
  ansSize++;
}
void Down()
{
  ans[ansSize][0] = 2;
  ansSize++;
}
void Left()
{
  ans[ansSize][0] = 1;
  ansSize++;
}
void Right()
{
  ans[ansSize][0] = 3;
  ansSize++;
}

vector<P> route;
void InitRoute1()
{
  route.clear();
  rep(i, n)
  {
    if (i % 2 == 0) {
      rep(j, n) { route.emplace_back(i, j); }
    }
    else {
      drep(j, n) { route.emplace_back(i, j); }
    }
  }
}

void InitRoute2()
{
  route.clear();
  rep(j, n)
  {
    if (j % 2 == 0) {
      rep(i, n) { route.emplace_back(i, j); }
    }
    else {
      drep(i, n) { route.emplace_back(i, j); }
    }
  }
}

void InitRoute3()
{
  route.clear();
  int ra = Rand() % 9 + 1;
  int a[n][n];
  rep(i, n) rep(j, n) a[i][j] = 0;
  int x = 0, y = 0;
  a[x][y] = 1;
  route.emplace_back(x, y);
  rep(i, ra)
  {
    while (y < n - 1 - i) {
      y++;
      route.emplace_back(x, y);
      a[x][y] = 1;
    }
    while (x < n - 1 - i) {
      x++;
      route.emplace_back(x, y);
      a[x][y] = 1;
    }
    while (y > i) {
      y--;
      route.emplace_back(x, y);
      a[x][y] = 1;
    }
    while (x > i + 1) {
      x--;
      route.emplace_back(x, y);
      a[x][y] = 1;
    }
  }

  // rep(i, n) {
  //   rep(j, n) { cout << a[i][j]; }
  //   cout << endl;
  // }

  int dir = 1;
  while (true) {
    if (dir == 1) {
      if (a[x][y + 1] == 0) {
        y++;
        route.emplace_back(x, y);
        a[x][y] = 1;
      }
      else {
        if (a[x + 1][y] == 0) {
          x++;
          route.emplace_back(x, y);
          a[x][y] = 1;
          dir = -1;
        }
        else {
          break;
        }
      }
    }
    else {
      if (a[x][y - 1] == 0) {
        y--;
        route.emplace_back(x, y);
        a[x][y] = 1;
      }
      else {
        if (a[x + 1][y] == 0) {
          x++;
          route.emplace_back(x, y);
          a[x][y] = 1;
          dir = 1;
        }
        else {
          break;
        }
      }
    }
  }
}

void InitRoute4()
{
  route.clear();
  int ra = Rand() % 9 + 1;
  int a[n][n];
  rep(i, n) rep(j, n) a[i][j] = 0;
  int x = 0, y = 0;
  a[x][y] = 1;
  route.emplace_back(x, y);
  rep(i, ra)
  {
    while (y < n - 1 - i) {
      y++;
      route.emplace_back(x, y);
      a[x][y] = 1;
    }
    while (x < n - 1 - i) {
      x++;
      route.emplace_back(x, y);
      a[x][y] = 1;
    }
    while (y > i) {
      y--;
      route.emplace_back(x, y);
      a[x][y] = 1;
    }
    while (x > i + 1) {
      x--;
      route.emplace_back(x, y);
      a[x][y] = 1;
    }
  }

  // rep(i, n) {
  //   rep(j, n) { cout << a[i][j]; }
  //   cout << endl;
  // }

  int dir = 1;
  if (a[x][y + 1] == 0) {
    y++;
    route.emplace_back(x, y);
    a[x][y] = 1;
  }
  while (true) {
    if (dir == 1) {
      if (a[x + 1][y] == 0) {
        x++;
        route.emplace_back(x, y);
        a[x][y] = 1;
      }
      else {
        if (a[x][y + 1] == 0) {
          y++;
          route.emplace_back(x, y);
          a[x][y] = 1;
          dir = -1;
        }
        else {
          break;
        }
      }
    }
    else {
      if (a[x - 1][y] == 0) {
        x--;
        route.emplace_back(x, y);
        a[x][y] = 1;
      }
      else {
        if (a[x][y + 1] == 0) {
          y++;
          route.emplace_back(x, y);
          a[x][y] = 1;
          dir = 1;
        }
        else {
          break;
        }
      }
    }
  }
}

struct Amount
{
  int idx;
  int d;
};

vector<Amount> spot[n][n];
void Method1(int ikichi1 = 300, int ikichi2 = 300)
{
  rep(i, n)
  {
    rep(j, n) { spot[i][j].clear(); }
  }
  InitH();
  ansSize = 0;
  {
    Amount a;
    int now = 0;
    int nowX = route[now].first;
    int nowY = route[now].second;
    int nokori = max(0, -h[nowX][nowY]);
    for (auto p : route) {
      int i = p.first;
      int j = p.second;
      if (h[i][j] > 0) {
        int hh = h[i][j];
        while (hh > 0) {
          if (nokori >= hh) {
            a.idx = now;
            a.d = hh;
            spot[i][j].emplace_back(a);
            nokori -= hh;
            hh = 0;
          }
          else {
            if (nokori > 0) {
              a.idx = now;
              a.d = nokori;
              spot[i][j].emplace_back(a);
              hh -= nokori;
              nokori = 0;
            }

            now++;
            nowX = route[now].first;
            nowY = route[now].second;
            nokori = max(0, -h[nowX][nowY]);
          }
        }
      }
    }
  }

  int d = 0;
  rep(i, route.size())
  {
    int x = route[i].first;
    int y = route[i].second;
    if (spot[x][y].size() > 0) {
      for (auto am : spot[x][y]) {
        if (am.idx > i) {
          GetLoad(x, y, am.d, d);
        }
      }
    }
    if (h[x][y] < 0 && d > 0) {
      if (d >= abs(h[x][y])) {
        PutLoad(x, y, abs(h[x][y]), d);
      }
      else {
        PutLoad(x, y, d, d);
      }
    }

    if (d >= ikichi1) {
      int iii = i;
      while (d > 0) {
        int nidx = iii + 1;
        while (nidx <= route.size() - 1) {
          int skip = 1;
          int ii = nidx;
          int xx = route[ii].first;
          int yy = route[ii].second;
          if (h[xx][yy] < 0 && d > 0) {
            skip = 0;
          }
          if (skip) {
            nidx++;
          }
          else {
            break;
          }
        }

        int nx = route[nidx].first;
        int ny = route[nidx].second;
        while (nx > x) {
          Down();
          x++;
        };
        while (nx < x) {
          Up();
          x--;
        }
        while (ny > y) {
          Right();
          y++;
        }
        while (ny < y) {
          Left();
          y--;
        }

        if (h[x][y] < 0 && d > 0) {
          if (d >= abs(h[x][y])) {
            PutLoad(x, y, abs(h[x][y]), d);
          }
          else {
            PutLoad(x, y, d, d);
          }
        }
      }
    }

    if (i < route.size() - 1) {
      int nidx = i + 1;
      while (nidx < route.size() - 1) {
        int skip = 1;
        int ii = nidx;
        int xx = route[ii].first;
        int yy = route[ii].second;
        if (spot[xx][yy].size() > 0) {
          for (auto am : spot[xx][yy]) {
            if (am.idx > ii) {
              skip = 0;
            }
          }
        }
        if (h[xx][yy] < 0 && d > 0) {
          skip = 0;
        }
        if (skip) {
          nidx++;
        }
        else {
          break;
        }
      }

      int nx = route[nidx].first;
      int ny = route[nidx].second;
      while (nx > x) {
        Down();
        x++;
      };
      while (nx < x) {
        Up();
        x--;
      }
      while (ny > y) {
        Right();
        y++;
      }
      while (ny < y) {
        Left();
        y--;
      }

      i = nidx - 1;
    }
  }

  drep(i, route.size())
  {
    int x = route[i].first;
    int y = route[i].second;
    if (spot[x][y].size() > 0) {
      for (auto am : spot[x][y]) {
        if (am.idx < i) {
          GetLoad(x, y, am.d, d);
        }
      }
    }
    if (h[x][y] < 0 && d > 0) {
      if (d >= abs(h[x][y])) {
        PutLoad(x, y, abs(h[x][y]), d);
      }
      else {
        PutLoad(x, y, d, d);
      }
    }

    if (d >= ikichi1) {
      int iii = i;
      while (d > 0) {
        int nidx = iii - 1;
        while (nidx >= 0) {
          int skip = 1;
          int ii = nidx;
          int xx = route[ii].first;
          int yy = route[ii].second;
          if (h[xx][yy] < 0 && d > 0) {
            skip = 0;
          }
          if (skip) {
            nidx--;
          }
          else {
            break;
          }
        }

        int nx = route[nidx].first;
        int ny = route[nidx].second;
        while (nx > x) {
          Down();
          x++;
        };
        while (nx < x) {
          Up();
          x--;
        }
        while (ny > y) {
          Right();
          y++;
        }
        while (ny < y) {
          Left();
          y--;
        }

        if (h[x][y] < 0 && d > 0) {
          if (d >= abs(h[x][y])) {
            PutLoad(x, y, abs(h[x][y]), d);
          }
          else {
            PutLoad(x, y, d, d);
          }
        }
      }
    }

    if (d == 0) {
      int ok = 1;
      drep(j, i + 1)
      {
        int x = route[j].first;
        int y = route[j].second;
        if (h[x][y] != 0) {
          ok = 0;
          break;
        }
      }
      if (ok) {
        break;
      }
    }

    if (i > 0) {
      int nidx = i - 1;
      while (nidx > 0) {
        int skip = 1;
        int ii = nidx;
        int xx = route[ii].first;
        int yy = route[ii].second;
        if (spot[xx][yy].size() > 0) {
          for (auto am : spot[xx][yy]) {
            if (am.idx < ii) {
              skip = 0;
            }
          }
        }
        if (h[xx][yy] < 0 && d > 0) {
          skip = 0;
        }
        if (skip) {
          nidx--;
        }
        else {
          break;
        }
      }

      int nx = route[nidx].first;
      int ny = route[nidx].second;
      while (nx > x) {
        Down();
        x++;
      };
      while (nx < x) {
        Up();
        x--;
      }
      while (ny > y) {
        Right();
        y++;
      }
      while (ny < y) {
        Left();
        y--;
      }

      i = nidx + 1;
    }
  }

  int score = CalcScore();
  if (score > real_ansScore) {
    CopyToReal();
    real_ansScore = score;
  }
}

// 解答出力
void Output(ofstream& ofs)
{
  if (mode == 0) {
    rep(i, ansSize)
    {
      if (ans[i][0] < 4) {
        cout << dc[ans[i][0]] << endl;
      }
      else {
        if (ans[i][1] > 0) cout << '+';
        cout << ans[i][1] << endl;
      }
    }
  }
  else {
    rep(i, ansSize)
    {
      if (ans[i][0] < 4) {
        ofs << dc[ans[i][0]] << endl;
      }
      else {
        if (ans[i][1] > 0) ofs << '+';
        ofs << ans[i][1] << endl;
      }
    }
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

  InitRoute1();
  Method1();

  InitRoute2();
  Method1();

  int loop = 0;
  while (true) {
    if (GetNowTime() > 1.8) {
      break;
    }
    int ra1 = Rand() % 1000;
    int ra2 = Rand() % 1000;
    int ra3 = Rand() % 100;
    if (ra3 < 10) {
      InitRoute1();
    }
    else if (ra3 < 20) {
      InitRoute2();
    }
    else if (ra3 < 60) {
      InitRoute3();
    }
    else {
      InitRoute4();
    }
    Method1(ra1, ra2);
    loop++;
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

int main()
{
  srand((unsigned)time(NULL));
  while (rand() % 100) {
    Rand();
  }

  mode = 0;

  if (mode == 0) {
    Solve(0);
  }
  else if (mode == 1) {
    ll sum = 0;
    srep(i, 0, 150)
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
