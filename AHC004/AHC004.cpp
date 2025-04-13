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
static double Rand01() {
  return (Rand() + 0.5) * (1.0 / UINT_MAX);
}

static const ll perfect_score = 100000000;
static const int n = 20;        
static const int max_patterns = 1000;      

int m;
string s[max_patterns];
vector<int> a[max_patterns];
int b[max_patterns] = {};
int ans[20][20];
ll ma;

bool f[max_patterns][20][20][2];
int countOK[max_patterns];

ll calc()
{
  int cnt = 0;
  rep(k, m)
  {
    int ok = 0;
    rep(i, n)
    {
      rep(j, n)
      {
        ok = 1;
        rep(l, b[k])
        {
          if (ans[i][(j + l) % n] != a[k][l]) {
            ok = 0;
            break;
          }
        }
        if (ok) break;
        ok = 1;
        rep(l, b[k])
        {
          if (ans[(i + l) % n][j] != a[k][l]) {
            ok = 0;
            break;
          }
        }
        if (ok) break;
      }
      if (ok) break;
    }
    cnt += ok;
  }
  ll res = 0;
  if (cnt < m) {
    res = perfect_score * cnt / m;
  }
  else {
    int cnt2 = 0;
    rep(i, n)
    {
      rep(j, n)
      {
        if (ans[i][j] == 0) cnt2++;
      }
    }
    res = perfect_score * 2 * n * n / (2 * n * n - cnt2);
  }
  return res;
}

ll calc2(int x, int y)
{
  rep(i, m)
  {
    rep(j, n)
    {
      countOK[i] -= f[i][j][y][0];
      countOK[i] -= f[i][j][y][1];
      f[i][j][y][0] = 0;
      f[i][j][y][1] = 0;
    }
    rep(k, n)
    {
      countOK[i] -= f[i][x][k][0];
      countOK[i] -= f[i][x][k][1];
      f[i][x][k][0] = 0;
      f[i][x][k][1] = 0;
    }
  }

  rep(k, m)
  {
    int ok = 0;
    srep(i, x, x + 1)
    {
      rep(j, n)
      {
        ok = 1;
        rep(l, b[k])
        {
          if (ans[i][(j + l) % n] != a[k][l]) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          countOK[k]++;
          f[k][i][j][0] = 1;
        }
        ok = 1;
        rep(l, b[k])
        {
          if (ans[(i + l) % n][j] != a[k][l]) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          countOK[k]++;
          f[k][i][j][1] = 1;
        }
      }
    }
  }

  rep(k, m)
  {
    int ok = 0;
    rep(i, n)
    {
      srep(j, y, y + 1)
      {
        if (i == x) continue;
        ok = 1;
        rep(l, b[k])
        {
          if (ans[i][(j + l) % n] != a[k][l]) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          countOK[k]++;
          f[k][i][j][0] = 1;
        }
        ok = 1;
        rep(l, b[k])
        {
          if (ans[(i + l) % n][j] != a[k][l]) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          countOK[k]++;
          f[k][i][j][1] = 1;
        }
      }
    }
  }

  int cnt = 0;
  rep(i, m) if (countOK[i]) cnt++;

  ll res = 0;
  if (cnt < m) {
    res = perfect_score * cnt / m;
  }
  else {
    int cnt2 = 0;
    rep(i, n)
    {
      rep(j, n)
      {
        if (ans[i][j] == 0) cnt2++;
      }
    }
    res = perfect_score * 2 * n * n / (2 * n * n - cnt2);
  }
  return res;
}

int main()
{
  string fileNameIfs = "in\\0000.txt";
  ifstream ifs(fileNameIfs.c_str());
  if (!ifs.is_open()) {  // 標準入力する
    int _n;
    cin >> _n >> m;
    rep(i, m)
    {
      cin >> s[i];
      b[i] = s[i].size();
      rep(j, b[i]) {
        a[i].push_back(s[i][j] - 'A' + 1);
      }
    }
  }
  else {  // ファイル入力する
    int _n;
    ifs >> _n >> m;
    rep(i, m)
    {
      ifs >> s[i];
      b[i] = s[i].size();
      rep(j, b[i]) { 
        a[i].push_back(s[i][j] - 'A' + 1);
      }
    }
  }

  rep(i, 20)
  {
    rep(j, 20) { 
      ans[i][j] = Rand() % 8 + 1; 
    }
  }

  clock_t start, end;
  start = clock();

  ma = calc();
  rep(i, n) {
    calc2(i, i);
  }

  int loop = 0;
  while (true) {
    loop++;
    int x = Rand() % n;
    int y = Rand() % n;
    int val = Rand() % 9;

    int keep = ans[x][y];
    ans[x][y] = val;

    ll tmp = calc2(x, y);
    if (tmp >= ma) {
      ma = tmp;
    }
    else {
      ans[x][y] = keep;
      calc2(x, y);
    }

    end = clock();
    if ((double)(end - start) / CLOCKS_PER_SEC > 2.9) break;
  }

  rep(i, n)
  {
    rep(j, n)
    {
      if (ans[i][j] == 0) {
        cout << '.';
      }
      else {
        cout << (char)(ans[i][j] + 'A' - 1);
      }
    }
    cout << endl;
  }
  return 0;
}
