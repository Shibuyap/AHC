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

const int INF = 1001001001;

namespace /* 乱数ライブラリ */
{
  static uint32_t randxor()
  {
    static uint32_t x = 123456789;
    static uint32_t y = 362436069;
    static uint32_t z = 521288629;
    static uint32_t w = 88675123;
    uint32_t t;

    t        = x ^ (x << 11);
    x        = y;
    y        = z;
    z        = w;
    return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
  }

  // 0以上1未満の小数をとる乱数
  static double rand01()
  {
    return (randxor() + 0.5) * (1.0 / UINT_MAX);
  }

  // 配列シャッフル
  void FisherYates(int* data, int n)
  {
    for (int i = n - 1; i >= 0; i--) {
      int j   = randxor() % (i + 1);
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

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

double TL = 2.8;
int mode;
clock_t startTime, endTime;

const int N = 600;
const int MAX_LB = 24;
int M, T, LA, LB;
int U[2000], V[2000];
int t[610];
int X[610], Y[610];
int a[2000];

vector<int> G[610];
vector<int> routes[610][610];

int ans[110000][4];
int ansCount;
int cCount;
int real_ans[110000][4];
int real_ansCount;
int real_cCount;

void CopyToReal()
{
  real_ansCount = ansCount;
  real_cCount   = cCount;
  rep(i, ansCount)
  {
    rep(j, 4)
    {
      real_ans[i][j] = ans[i][j];
    }
  }
}

void CopyToAns()
{
  ansCount = real_ansCount;
  cCount   = real_cCount;
  rep(i, ansCount)
  {
    rep(j, 4)
    {
      ans[i][j] = real_ans[i][j];
    }
  }
}

double GetNowTime()
{
  endTime        = clock();
  double nowTime = ((double)endTime - startTime) / CLOCKS_PER_SEC;
  return nowTime;
}

// 複数ケース回すときに内部状態を初期値に戻す
void SetUp()
{
  rep(i, 610)
  {
    G[i].clear();
  }
  rep(i, 610)
  {
    rep(j, 610)
    {
      routes[i][j].clear();
    }
  }
}

void InitRoutes()
{
  int visited[610];
  int prev[610];
  queue<int> que;
  rep(i, N)
  {
    rep(j, N)
    {
      visited[j] = 0;
      prev[j]    = -1;
    }
    visited[i] = 1;
    que.push(i);
    while (que.size()) {
      int x = que.front();
      que.pop();
      for (auto y : G[x]) {
        if (visited[y] == 0) {
          visited[y] = 1;
          prev[y]    = x;
          que.push(y);
        }
      }
    }

    rep(j, N)
    {
      int x = j;
      while (x != i) {
        routes[i][j].push_back(x);
        x = prev[x];
      }
      reverse(routes[i][j].begin(), routes[i][j].end());
    }
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
    int nnn;
    cin >> nnn >> M >> T >> LA >> LB;
    rep(i, M)
    {
      cin >> U[i] >> V[i];
    }
    rep(i, T)
    {
      cin >> t[i];
    }
    rep(i, N)
    {
      cin >> X[i] >> Y[i];
    }
  }
  // ファイル入力する
  else {
    int nnn;
    ifs >> nnn >> M >> T >> LA >> LB;
    rep(i, M)
    {
      ifs >> U[i] >> V[i];
    }
    rep(i, T)
    {
      ifs >> t[i];
    }
    rep(i, N)
    {
      ifs >> X[i] >> Y[i];
    }
  }

  rep(i, M)
  {
    G[U[i]].push_back(V[i]);
    G[V[i]].push_back(U[i]);
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
ll CalcScore()
{
  ll res = cCount;
  return res;
}

void MakeA2(const vector<int>& route)
{
  int i = 0;
  while (i < LA) {
    int f[N] = {};
    for (auto x : route) {
      if (i == LA) break;
      if (f[x] == 0) {
        a[i] = x;
        f[x] = 1;
        i++;
      }
    }
    rep(j, N) f[j] = 0;
    drep(j, route.size())
    {
      if (i == LA) break;
      int x = route[j];
      if (f[x] == 0) {
        a[i] = x;
        f[x] = 1;
        i++;
      }
    }
  }
}

// 初期解生成
void Initialize(const vector<int>& route)
{
  ansCount = 0;
  cCount   = 0;

  int b[MAX_LB];
  rep(i, LB)
  {
    b[i] = -1;
  }

  rep(i, route.size())
  {
    int r      = route[i];
    int exists = 0;
    rep(j, LB)
    {
      if (b[j] == r) {
        exists = 1;
        break;
      }
    }

    if (exists == 0) {
      rep(j, LA)
      {
        if (a[j] == r) {
          int len          = min(LB, LA - j);
          ans[ansCount][0] = 0;
          ans[ansCount][1] = len;
          ans[ansCount][2] = j;
          ans[ansCount][3] = 0;
          ansCount++;
          cCount++;
          break;
        }
      }
    }

    ans[ansCount][0] = 1;
    ans[ansCount][1] = r;
    ansCount++;
  }
}

vector<int> argN[610];
void Method1(const vector<int>& route)
{
  ansCount = 0;
  cCount   = 0;

  int b[MAX_LB];
  rep(i, LB)
  {
    b[i] = -1;
  }

  rep(i, route.size())
  {
    int r      = route[i];
    int exists = 0;
    rep(j, LB)
    {
      if (b[j] == r) {
        exists = 1;
        break;
      }
    }

    if (exists == 0) {
      int pos = argN[r][0];
      int len = min(LB, LA - pos);
      rep(j, len)
      {
        b[j] = a[pos + j];
      }
      ans[ansCount][0] = 0;
      ans[ansCount][1] = len;
      ans[ansCount][2] = pos;
      ans[ansCount][3] = 0;
      ansCount++;
      cCount++;
    }

    ans[ansCount][0] = 1;
    ans[ansCount][1] = r;
    ansCount++;
  }
}

void Method2(const vector<int>& route)
{
  ansCount = 0;
  cCount   = 0;

  int b[MAX_LB];
  rep(i, LB)
  {
    b[i] = -1;
  }

  rep(i, route.size())
  {
    int r      = route[i];
    int exists = 0;
    rep(j, LB)
    {
      if (b[j] == r) {
        exists = 1;
        break;
      }
    }

    if (exists == 0) {
      int ma   = -1;
      int pos2 = -1;
      int len2 = -1;
      {
        int pos = argN[r][0];
        while (LA - pos < LB) {
          pos--;
        }
        int len = min(LB, LA - pos);

        int now = i;
        while (now < route.size() - 1) {
          int rr     = route[now + 1];
          int exists = 0;
          rep(j, LB)
          {
            if (a[pos + j] == rr) {
              exists = 1;
              break;
            }
          }
          if (exists) {
            now++;
          }
          else {
            break;
          }
        }
        pos2 = pos;
        len2 = len;
        ma   = now;
      }

      {
        int pos = argN[r][0] - LB + 1;
        pos     = max(pos, 0);
        int len = min(LB, LA - pos);

        int now = i;
        while (now < route.size() - 1) {
          int rr     = route[now + 1];
          int exists = 0;
          rep(j, LB)
          {
            if (a[pos + j] == rr) {
              exists = 1;
              break;
            }
          }
          if (exists) {
            now++;
          }
          else {
            break;
          }
        }
        if (ma < now) {
          pos2 = pos;
          len2 = len;
          ma   = now;
        }
      }

      rep(j, len2)
      {
        b[j] = a[pos2 + j];
      }
      ans[ansCount][0] = 0;
      ans[ansCount][1] = len2;
      ans[ansCount][2] = pos2;
      ans[ansCount][3] = 0;
      ansCount++;
      cCount++;
    }

    ans[ansCount][0] = 1;
    ans[ansCount][1] = r;
    ansCount++;
  }
}

// 解答出力
void Output(ofstream& ofs)
{
  if (mode == 0) {
    rep(i, LA)
    {
      cout << a[i] << ' ';
    }
    cout << endl;
    rep(i, ansCount)
    {
      if (ans[i][0] == 0) {
        cout << "s " << ans[i][1] << ' ' << ans[i][2] << ' ' << ans[i][3] << endl;
      }
      else {
        cout << "m " << ans[i][1] << endl;
      }
    }
  }
  else {
    rep(i, LA)
    {
      ofs << a[i] << ' ';
    }
    ofs << endl;
    rep(i, ansCount)
    {
      if (ans[i][0] == 0) {
        ofs << "s " << ans[i][1] << ' ' << ans[i][2] << ' ' << ans[i][3] << endl;
      }
      else {
        ofs << "m " << ans[i][1] << endl;
      }
    }
  }
}

vector<int> GetRoute()
{
  int now = 0;
  vector<int> res;
  rep(i, T)
  {
    for (auto x : routes[now][t[i]]) {
      res.push_back(x);
    }
    now = t[i];
  }
  return res;
}

void MakeA1DFS(int x, vector<int>& route, vector<int>& visited, int order = 0)
{
  if (visited[x]) return;
  visited[x] = 1;
  route.push_back(x);
  if (order == 0) {
    for (auto y : G[x]) {
      MakeA1DFS(y, route, visited);
    }
  }
  else {
    drep(i, G[x].size())
    {
      int y = G[x][i];
      MakeA1DFS(y, route, visited);
    }
  }
}

void MakeA1()
{
  int center = -1;
  int mi     = INF;
  rep(i, N)
  {
    int tmp = abs(500 - X[i]) + abs(500 - Y[i]);
    if (tmp < mi) {
      center = i;
      mi     = tmp;
    }
  }

  vector<int> route;
  vector<int> visited(N, 0);
  MakeA1DFS(center, route, visited);

  int i = 0;
  for (auto x : route) {
    a[i] = x;
    i++;
  }

  int lastPoint = route.back();
  route.clear();
  rep(i, N) visited[i] = 0;
  MakeA1DFS(center, route, visited, 1);
  srep(j, 1, route.size())
  {
    if (i == LA) break;
    a[i] = route[j];
    i++;
  }

  rep(j, route.size())
  {
    if (i == LA) break;
    a[i] = route[j];
    i++;
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

  InitRoutes();

  vector<int> route = GetRoute();

  MakeA1();
  // MakeA2(route);

  // 初期解生成
  Initialize(route);

  CopyToReal();

  int loop  = 0;
  int loop2 = 0;
  while (GetNowTime() < TL) {
    loop++;
    int raMode = randxor() % 2;
    int ra1    = randxor() % LA;
    int ra2    = randxor() % LA;
    int raLen  = randxor() % 10 + 2;

    if (raMode == 0) {
      swap(a[ra1], a[ra2]);
    }
    else {
      if (ra1 + raLen > LA) continue;
      reverse(a + ra1, a + ra1 + raLen);
    }

    rep(i, N) argN[i].clear();
    rep(i, LA)
    {
      argN[a[i]].push_back(i);
    }

    // Method1(route);
    Method2(route);

    int tmpScore = CalcScore();
    if (tmpScore <= real_cCount) {
      if (tmpScore < real_cCount) {
        // cout << tmpScore << endl;
      }
      loop2++;
      CopyToReal();
    }
    else {
      if (raMode == 0) {
        swap(a[ra1], a[ra2]);
      }
      else {
        reverse(a + ra1, a + ra1 + raLen);
      }
    }
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
    cout << "loop = " << loop << ", loop2 = " << loop2;
    cout << endl;
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
