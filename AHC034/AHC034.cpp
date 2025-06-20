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

using namespace std;
typedef long long int ll;
typedef pair<int, int> P;

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

  // 配列シャッフル
  void FisherYates(int* data, int n)
  {
    for (int i = n - 1; i >= 0; i--) {
      int j = Rand() % (i + 1);
      int tmp = data[i];
      data[i] = data[j];
      data[j] = tmp;
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

const int n = 20;
int baseH[20][20];
int h[20][20];
int ans[100000][2];
int ansSize;

int best_ans[100000][2];
int best_ansSize;
int best_ansScore;

void InitH()
{
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) { h[i][j] = baseH[i][j]; }
  }
}

void CopyToAns()
{
  ansSize = best_ansSize;
  for (int i = 0; i < ansSize; ++i) {
    for (int j = 0; j < 2; ++j) { ans[i][j] = best_ans[i][j]; }
  }
}

void CopyToBest()
{
  best_ansSize = ansSize;
  for (int i = 0; i < ansSize; ++i) {
    for (int j = 0; j < 2; ++j) { best_ans[i][j] = ans[i][j]; }
  }
}

// 複数ケース回すときに内部状態を初期値に戻す
void SetUp() { ansSize = 0; }

// 入力受け取り
void Input(int problemNum)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << problemNum << ".txt";
  ifstream ifs(oss.str());

  // 標準入力する
  if (!ifs.is_open()) {
    int _n;
    cin >> _n;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) { cin >> h[i][j]; }
    }
  }
  // ファイル入力する
  else {
    int _n;
    ifs >> _n;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) { ifs >> h[i][j]; }
    }
  }

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) { baseH[i][j] = h[i][j]; }
  }
}

// 出力ファイルストリームオープン
void OpenOfs(int probNum, ofstream& ofs)
{
  if (mode != 0) {
    string fileNameOfs = "./out/";
    string strNum;
    for (int i = 0; i < 4; ++i) {
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
  for (int i = 0; i < ansSize; ++i) {
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
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) { base += abs(baseH[i][j]); }
  }
  ll result = round(1e9 * base / CalcCost());
  return result;
}

bool IsRowEmpty(int x)
{
  for (int j = 0; j < n; ++j) {
    if (h[x][j] != 0) {
      return false;
    }
  }
  return true;
}

bool IsValidCell(int x, int y, int a[n][n])
{
  return x >= 0 && x < n && y >= 0 && y < n && a[x][y] == 0;
}

// 初期解生成
void Initialize()
{
  int x = 0, y = 0;
  int d = 0;
  int dir = 1;
  while (true) {
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
        if (x == n - 1) { break; }
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
        if (x == n - 1) { break; }
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
    if (IsRowEmpty(x)) {
      if (x == 0) { break; }
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
        if (!IsRowEmpty(x)) {
          dir *= -1;
          continue;
        }
        if (x == 0) { break; }
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
        if (!IsRowEmpty(x)) {
          dir *= -1;
          continue;
        }
        if (x == 0) { break; }
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

  CopyToBest();
  best_ansScore = CalcScore();
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

void MoveTo(int& currentX, int& currentY, int targetX, int targetY)
{
  while (targetX > currentX) {
    Down();
    currentX++;
  }
  while (targetX < currentX) {
    Up();
    currentX--;
  }
  while (targetY > currentY) {
    Right();
    currentY++;
  }
  while (targetY < currentY) {
    Left();
    currentY--;
  }
}

vector<P> route;
void InitRoute1()
{
  route.clear();
  for (int i = 0; i < n; ++i) {
    if (i % 2 == 0) {
      for (int j = 0; j < n; ++j) { route.emplace_back(i, j); }
    }
    else {
      for (int j = n - 1; j >= 0; --j) { route.emplace_back(i, j); }
    }
  }
}

void InitRoute2()
{
  route.clear();
  for (int j = 0; j < n; ++j) {
    if (j % 2 == 0) {
      for (int i = 0; i < n; ++i) { route.emplace_back(i, j); }
    }
    else {
      for (int i = n - 1; i >= 0; --i) { route.emplace_back(i, j); }
    }
  }
}

void GenerateSpiralRoute(int& x, int& y, int a[n][n], int spiralDepth)
{
  for (int i = 0; i < spiralDepth; ++i) {
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
}

void InitRoute3()
{
  route.clear();
  int spiralDepthRand = Rand() % 9 + 1;
  int a[n][n];
  for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) a[i][j] = 0;
  int x = 0, y = 0;
  a[x][y] = 1;
  route.emplace_back(x, y);

  GenerateSpiralRoute(x, y, a, spiralDepthRand);

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
  int spiralDepthRand = Rand() % 9 + 1;
  int a[n][n];
  for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) a[i][j] = 0;
  int x = 0, y = 0;
  a[x][y] = 1;
  route.emplace_back(x, y);

  GenerateSpiralRoute(x, y, a, spiralDepthRand);

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
void Method1(int threshold1 = 300, int threshold2 = 300)
{
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) { spot[i][j].clear(); }
  }
  InitH();
  ansSize = 0;
  {
    Amount a;
    int now = 0;
    int nowX = route[now].first;
    int nowY = route[now].second;
    int remain = max(0, -h[nowX][nowY]);
    for (auto p : route) {
      int i = p.first;
      int j = p.second;
      if (h[i][j] > 0) {
        int height = h[i][j];
        while (height > 0) {
          if (remain >= height) {
            a.idx = now;
            a.d = height;
            spot[i][j].emplace_back(a);
            remain -= height;
            height = 0;
          }
          else {
            if (remain > 0) {
              a.idx = now;
              a.d = remain;
              spot[i][j].emplace_back(a);
              height -= remain;
              remain = 0;
            }

            now++;
            nowX = route[now].first;
            nowY = route[now].second;
            remain = max(0, -h[nowX][nowY]);
          }
        }
      }
    }
  }

  int d = 0;
  for (int i = 0; i < route.size(); ++i) {
    int x = route[i].first;
    int y = route[i].second;
    if (spot[x][y].size() > 0) {
      for (auto amt : spot[x][y]) {
        if (amt.idx > i) {
          GetLoad(x, y, amt.d, d);
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

    if (d >= threshold1) {
      int currentIdx = i;
      while (d > 0) {
        int nidx = currentIdx + 1;
        while (nidx <= route.size() - 1) {
          int checkX = route[nidx].first;
          int checkY = route[nidx].second;
          if (h[checkX][checkY] < 0 && d > 0) {
            break;
          }
          nidx++;
        }

        int nx = route[nidx].first;
        int ny = route[nidx].second;
        MoveTo(x, y, nx, ny);

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
        int checkX = route[ii].first;
        int checkY = route[ii].second;
        if (spot[checkX][checkY].size() > 0) {
          for (auto amt : spot[checkX][checkY]) {
            if (amt.idx > ii) {
              skip = 0;
            }
          }
        }
        if (h[checkX][checkY] < 0 && d > 0) {
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
      MoveTo(x, y, nx, ny);

      i = nidx - 1;
    }
  }

  for (int i = route.size() - 1; i >= 0; --i) {
    int x = route[i].first;
    int y = route[i].second;
    if (spot[x][y].size() > 0) {
      for (auto amt : spot[x][y]) {
        if (amt.idx < i) {
          GetLoad(x, y, amt.d, d);
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

    if (d >= threshold1) {
      int currentIdx = i;
      while (d > 0) {
        int nidx = currentIdx - 1;
        while (nidx >= 0) {
          int checkX = route[nidx].first;
          int checkY = route[nidx].second;
          if (h[checkX][checkY] < 0 && d > 0) {
            break;
          }
          nidx--;
        }

        int nx = route[nidx].first;
        int ny = route[nidx].second;
        MoveTo(x, y, nx, ny);

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
      for (int j = i; j >= 0; --j) {
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
        int checkX = route[ii].first;
        int checkY = route[ii].second;
        if (spot[checkX][checkY].size() > 0) {
          for (auto amt : spot[checkX][checkY]) {
            if (amt.idx < ii) {
              skip = 0;
            }
          }
        }
        if (h[checkX][checkY] < 0 && d > 0) {
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
      MoveTo(x, y, nx, ny);

      i = nidx + 1;
    }
  }

  int score = CalcScore();
  if (score > best_ansScore) {
    CopyToBest();
    best_ansScore = score;
  }
}

// 解答出力
void Output(ofstream& ofs)
{
  if (mode == 0) {
    for (int i = 0; i < ansSize; ++i) {
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
    for (int i = 0; i < ansSize; ++i) {
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
  start_timer();

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

  for (int loop = 0; loop < 100000; ++loop) {
    int randThreshold1 = Rand() % 1000;
    int randThreshold2 = Rand() % 1000;
    int routeType = Rand() % 100;
    if (routeType < 10) {
      InitRoute1();
    }
    else if (routeType < 20) {
      InitRoute2();
    }
    else if (routeType < 60) {
      InitRoute3();
    }
    else {
      InitRoute4();
    }
    Method1(randThreshold1, randThreshold2);
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
  mode = 1;

  if (mode == 0) {
    Solve(0);
  }
  else if (mode == 1) {
    ll sum = 0;
    for (int i = 0; i < 10; ++i) {
      ll score = Solve(i);
      sum += score;
      cout << "num = " << i << ", ";
      cout << "score = " << score << ", ";
      cout << "sum = " << sum << ", ";
      cout << "time = " << get_elapsed_time() << " sec" << endl;
    }
  }

  return 0;
}
