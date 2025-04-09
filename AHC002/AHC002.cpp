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
#define rep(i,n) for(int i = 0; i < (n); ++i)
#define srep(i,s,t) for(int i = s; i < t; ++i)
#define drep(i,n) for(int i = (n)-1; i >= 0; --i)
using namespace std;
// using namespace atcoder;
typedef long long int ll;
typedef pair<int, int> P;
#define yn {puts("Yes");}else{puts("No");}
#define MAX_N 200005

static uint32_t Rand()
{
  static uint32_t x = 123456789;
  static uint32_t y = 362436069;
  static uint32_t z = 521288629;
  static uint32_t w = 88675123;
  uint32_t t;

  t = x ^ (x << 11);
  x = y; y = z; z = w;
  return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
}

int nxt[24][4] = { {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 1, 3}, {0, 2, 3, 1}, {0, 3, 1, 2}, {0, 3, 2, 1},
                   {1, 0, 2, 3}, {1, 0, 3, 2}, {1, 2, 0, 3}, {1, 2, 3, 0}, {1, 3, 0, 2}, {1, 3, 2, 0},
                   {2, 0, 1, 3}, {2, 0, 3, 1}, {2, 1, 0, 3}, {2, 1, 3, 0}, {2, 3, 0, 1}, {2, 3, 1, 0},
                   {3, 0, 1, 2}, {3, 0, 2, 1}, {3, 1, 0, 2}, {3, 1, 2, 0}, {3, 2, 0, 1}, {3, 2, 1, 0} };

int dx[4] = { -1, 0, 1, 0 };
int dy[4] = { 0, -1, 0, 1 };
char nxtc[4] = { 'U','L','D','R' };

const int n = 50;
int si, sj;
int f[60][60];
int a[60][60];
int g[60][60];

P h[60][60];

string ans;
int ma = 0;

int calcScore(string s)
{
  int x = si, y = sj;
  int res = a[x][y];
  int m = s.size();
  rep(i, m)
  {
    if (s[i] == 'U') x--;
    if (s[i] == 'D') x++;
    if (s[i] == 'L') y--;
    if (s[i] == 'R') y++;
    res += a[x][y];
  }
  return res;
}



int main()
{
  srand((unsigned)time(NULL));
  while (rand() % 10 != 0) {
    Rand();
  }

  string fileNameIfs = "1120.txt";
  const char* cstrIfs = fileNameIfs.c_str();
  ifstream ifs(cstrIfs);
  if (!ifs.is_open()) { // 標準入力する
    cin >> si >> sj;
    rep(i, n)
    {
      rep(j, n)
      {
        cin >> f[i][j];
      }
    }
    rep(i, n)
    {
      rep(j, n)
      {
        cin >> a[i][j];
      }
    }
  }
  else { // ファイル入力する
    ifs >> si >> sj;
    rep(i, n)
    {
      rep(j, n)
      {
        ifs >> f[i][j];
      }
    }
    rep(i, n)
    {
      rep(j, n)
      {
        ifs >> a[i][j];
      }
    }
  }



  clock_t start, end;
  start = clock();

  rep(i, n)
  {
    rep(j, n)
    {
      int num = f[i][j];
      h[i][j].first = -1;
      h[i][j].second = -1;
      rep(k, 4)
      {
        int nx = i + dx[k];
        int ny = j + dy[k];
        if (0 <= nx && nx < n && 0 <= ny && ny < n && f[nx][ny] == num) {
          h[i][j].first = nx;
          h[i][j].second = ny;
        }
      }
    }
  }

  int loop = 0;
  // 初期解
  rep(_, 10000)
  {
    loop++;
    rep(i, n) rep(j, n) g[i][j] = 0;
    int x = si;
    int y = sj;
    int tmp = a[x][y];

    g[x][y] = 1;
    if (h[x][y].first != -1) g[h[x][y].first][h[x][y].second] = 1;

    string t;

    while (true) {
      int ra = Rand() % 24;
      int ok = 0;

      rep(i, 4)
      {
        int nx = x + dx[nxt[ra][i]];
        int ny = y + dy[nxt[ra][i]];
        if (0 <= nx && nx < n && 0 <= ny && ny < n && g[nx][ny] == 0) {
          ok = 1;
          tmp += a[nx][ny];
          g[nx][ny] = 1;
          if (h[nx][ny].first != -1) g[h[nx][ny].first][h[nx][ny].second] = 1;
          t += nxtc[nxt[ra][i]];
          x = nx;
          y = ny;
          break;
        }
      }

      if (ok == 0) break;
    }


    if (tmp > ma) {
      ma = tmp;
      ans = t;
    }

    end = clock();
    if ((double)(end - start) / CLOCKS_PER_SEC > 1.9)break;
  }


  while (true) {
    loop++;
    rep(i, n) rep(j, n) g[i][j] = 0;
    int x = si;
    int y = sj;
    int tmp = a[x][y];

    g[x][y] = 1;
    if (h[x][y].first != -1) g[h[x][y].first][h[x][y].second] = 1;

    string t;
    int m = rand() % ans.size() + 1;
    t = ans.substr(0, m);
    rep(i, m)
    {
      if (t[i] == 'U') x--;
      if (t[i] == 'D') x++;
      if (t[i] == 'L') y--;
      if (t[i] == 'R') y++;
      tmp += a[x][y];
      g[x][y] = 1;
      if (h[x][y].first != -1) g[h[x][y].first][h[x][y].second] = 1;
    }

    while (true) {
      int ra = Rand() % 24;
      int ok = 0;

      rep(i, 4)
      {
        int nx = x + dx[nxt[ra][i]];
        int ny = y + dy[nxt[ra][i]];
        if (0 <= nx && nx < n && 0 <= ny && ny < n && g[nx][ny] == 0) {
          ok = 1;
          tmp += a[nx][ny];
          g[nx][ny] = 1;
          if (h[nx][ny].first != -1) g[h[nx][ny].first][h[nx][ny].second] = 1;
          t += nxtc[nxt[ra][i]];
          x = nx;
          y = ny;
          break;
        }
      }

      if (ok == 0) break;
    }


    if (tmp > ma) {
      ma = tmp;
      ans = t;
    }

    end = clock();
    if ((double)(end - start) / CLOCKS_PER_SEC > 1.0)break;
  }


  while (true) {
    loop++;

    rep(i, n) rep(j, n) g[i][j] = 0;
    int x = si;
    int y = sj;
    int tmp = a[x][y];

    g[x][y] = 1;
    if (h[x][y].first != -1) g[h[x][y].first][h[x][y].second] = 1;

    string t;
    int m = ans.size();
    t = ans;

    vector<int> xx, yy;
    xx.push_back(x);
    yy.push_back(y);

    rep(i, m)
    {
      if (t[i] == 'U') x--;
      if (t[i] == 'D') x++;
      if (t[i] == 'L') y--;
      if (t[i] == 'R') y++;
      tmp += a[x][y];
      g[x][y] = 1;
      if (h[x][y].first != -1) g[h[x][y].first][h[x][y].second] = 1;
      xx.push_back(x);
      yy.push_back(y);
    }

    int left = Rand() % (m - 40) + 10;
    int right = left + 1 + Rand() % 20;
    int tmp2 = 0;
    srep(i, left + 1, right)
    {
      x = xx[i], y = yy[i];
      tmp -= a[x][y];
      tmp2 += a[x][y];
      g[x][y] = 0;
      if (h[x][y].first != -1) g[h[x][y].first][h[x][y].second] = 0;
    }

    string t1, t2, t3;
    rep(i, left) t1 += t[i];
    srep(i, left, right) t2 += t[i];
    srep(i, right, m) t3 += t[i];

    int sx = xx[left], sy = yy[left];
    int gx = xx[right], gy = yy[right];
    // cout << ans << endl;
    rep(_, 100)
    {
      x = sx; y = sy;
      vector<int> xxx, yyy;
      int tmp3 = 0;
      string ttt;
      while (x != gx || y != gy) {
        int ra = Rand() % 24;
        int ok = 0;

        rep(i, 4)
        {
          int nx = x + dx[nxt[ra][i]];
          int ny = y + dy[nxt[ra][i]];
          if (nx == gx && ny == gy && f[x][y] != f[nx][ny]) {
            ok = 2;
            x = gx;
            y = gy;
            ttt += nxtc[nxt[ra][i]];
            break;
          }
          if (0 <= nx && nx < n && 0 <= ny && ny < n && g[nx][ny] == 0) {
            ok = 1;
            tmp3 += a[nx][ny];
            g[nx][ny] = 1;
            xxx.push_back(nx);
            yyy.push_back(ny);
            if (h[nx][ny].first != -1) g[h[nx][ny].first][h[nx][ny].second] = 1;
            ttt += nxtc[nxt[ra][i]];
            x = nx;
            y = ny;
            break;
          }
        }

        if (ok == 0) break;
      }

      if (x == gx && y == gy && tmp3 > tmp2) {
        tmp2 = tmp3;
        t2 = ttt;
      }
      rep(i, xxx.size())
      {
        x = xxx[i]; y = yyy[i];
        g[x][y] = 0;
        if (h[x][y].first != -1) g[h[x][y].first][h[x][y].second] = 0;
      }
    }


    if (tmp + tmp2 > ma) {
      ma = tmp + tmp2;
      ans = t1 + t2 + t3;
    }

    // cout << ans << endl;


    end = clock();
    if ((double)(end - start) / CLOCKS_PER_SEC > 1.9)break;
  }

  // cout << loop << endl;
  // cout << ma << ' ' << calcScore(ans) << endl;
  cout << ans << endl;
  return 0;
}


