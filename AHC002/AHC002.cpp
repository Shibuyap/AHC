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

#define rep(i,n) for(int i = 0; i < (n); ++i)
#define srep(i,s,t) for(int i = s; i < t; ++i)
#define drep(i,n) for(int i = (n)-1; i >= 0; --i)
using namespace std;
typedef long long int ll;
typedef pair<int, int> P;

#define MAX_N 200005

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

struct Point {
  int x;
  int y;
};

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
int value_grid[60][60];
P h[60][60];

int visited2[n * n];
int v2_counter;

const int KOUHO_SIZE = 100;
string kouho_paths[KOUHO_SIZE];
int kouho_scores[KOUHO_SIZE];

string best_path;
int best_score = 0;

int calcScore(const string& path)
{
  int x = si, y = sj;
  int res = value_grid[x][y];
  int m = path.size();
  rep(i, m)
  {
    if (path[i] == 'U') x--;
    if (path[i] == 'D') x++;
    if (path[i] == 'L') y--;
    if (path[i] == 'R') y++;
    res += value_grid[x][y];
  }
  return res;
}

bool is_out_of_range(int x, int y)
{
  if (x < 0 || n <= x || y < 0 || n <= y) return true;
  return false;
}

void sort_kouho(int size)
{
  vector<pair<int, string>> kouho;
  rep(i, size)
  {
    if (kouho_scores[i] > 0) {
      kouho.push_back(make_pair(kouho_scores[i], kouho_paths[i]));
    }
  }
  sort(kouho.begin(), kouho.end(), greater<pair<int, string>>());
  rep(i, size)
  {
    if (i < kouho.size()) {
      kouho_scores[i] = kouho[i].first;
      kouho_paths[i] = kouho[i].second;
    }
    else {
      kouho_scores[i] = 0;
      kouho_paths[i] = "";
    }
  }
}

static void input_data() {
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
        cin >> value_grid[i][j];
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
        ifs >> value_grid[i][j];
      }
    }
  }
}

int main()
{
  input_data();

  start_timer();

  // h初期化
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
        if (!is_out_of_range(nx, ny) && f[nx][ny] == num) {
          h[i][j].first = nx;
          h[i][j].second = ny;
        }
      }
    }
  }

  rep(i, KOUHO_SIZE) {
    kouho_scores[i] = 0;
  }

  int loop1 = 0;
  // 初期解
  rep(_, 100000)
  {
    loop1++;
    v2_counter++;
    int x = si;
    int y = sj;
    int tmp = value_grid[x][y];

    visited2[f[x][y]] = v2_counter;

    string t;

    while (true) {
      int ra = Rand() % 24;
      int ok = 0;

      rep(i, 4)
      {
        int nx = x + dx[nxt[ra][i]];
        int ny = y + dy[nxt[ra][i]];
        if (!is_out_of_range(nx, ny) && visited2[f[nx][ny]] != v2_counter) {
          ok = 1;
          tmp += value_grid[nx][ny];
          visited2[f[nx][ny]] = v2_counter;
          t += nxtc[nxt[ra][i]];
          x = nx;
          y = ny;
          break;
        }
      }

      if (ok == 0) break;
    }

    int num = _ % KOUHO_SIZE;
    if (tmp > kouho_scores[num]) {
      kouho_scores[num] = tmp;
      kouho_paths[num] = t;
    }

    if (get_elapsed_time() > 1.9)
    {
      break;
    } 
  }

  int loop2 = 0;
  while (true) {
    loop2++;
    v2_counter++;
    int x = si;
    int y = sj;
    int tmp = value_grid[x][y];

    visited2[f[x][y]] = v2_counter;

    int num = Rand() % KOUHO_SIZE;

    string t;
    int m = rand() % kouho_paths[num].size() + 1;
    t = kouho_paths[num].substr(0, m);
    rep(i, m)
    {
      if (t[i] == 'U') x--;
      if (t[i] == 'D') x++;
      if (t[i] == 'L') y--;
      if (t[i] == 'R') y++;
      tmp += value_grid[x][y];
      visited2[f[x][y]] = v2_counter;
    }

    while (true) {
      int ra = Rand() % 24;
      int ok = 0;

      rep(i, 4)
      {
        int nx = x + dx[nxt[ra][i]];
        int ny = y + dy[nxt[ra][i]];
        if (!is_out_of_range(nx, ny) && visited2[f[nx][ny]] != v2_counter) {
          ok = 1;
          tmp += value_grid[nx][ny];
          visited2[f[nx][ny]] = v2_counter;
          t += nxtc[nxt[ra][i]];
          x = nx;
          y = ny;
          break;
        }
      }

      if (ok == 0) break;
    }

    if (tmp > kouho_scores[num]) {
      kouho_scores[num] = tmp;
      kouho_paths[num] = t;
    }

    if (get_elapsed_time() > 1.0)
    {
      break;
    }
  }

  sort_kouho(KOUHO_SIZE);
  best_path = kouho_paths[0];
  best_score = kouho_scores[0];

  bool is_sorted = false;

  int xxx[10000], yyy[10000];
  int xxx_count, yyy_count;
  int loop3 = 0;
  while (true) {
    loop3++;

    v2_counter++;
    int x = si;
    int y = sj;
    int tmp = value_grid[x][y];

    visited2[f[x][y]] = v2_counter;

    int num = Rand() % 10;
    if (is_sorted) {
      num = 0;
    }

    int m = kouho_paths[num].size();
    string t = kouho_paths[num];

    vector<int> xx, yy;
    xx.push_back(x);
    yy.push_back(y);

    rep(i, m)
    {
      if (t[i] == 'U') x--;
      if (t[i] == 'D') x++;
      if (t[i] == 'L') y--;
      if (t[i] == 'R') y++;
      tmp += value_grid[x][y];
      visited2[f[x][y]] = v2_counter;
      xx.push_back(x);
      yy.push_back(y);
    }

    int left = Rand() % (m - 40) + 10;
    int right = left + 1 + Rand() % 20;
    int tmp2 = 0;
    srep(i, left + 1, right)
    {
      x = xx[i], y = yy[i];
      tmp -= value_grid[x][y];
      tmp2 += value_grid[x][y];
      visited2[f[x][y]] = -1;
    }

    string t1, t2, t3;
    rep(i, left) t1 += t[i];
    srep(i, left, right) t2 += t[i];
    srep(i, right, m) t3 += t[i];

    int sx = xx[left], sy = yy[left];
    int gx = xx[right], gy = yy[right];

    rep(_, 100)
    {
      xxx_count = 0;
      yyy_count = 0;
      x = sx; y = sy;
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
          if (!is_out_of_range(nx, ny) && visited2[f[nx][ny]] != v2_counter) {
            ok = 1;
            tmp3 += value_grid[nx][ny];
            visited2[f[nx][ny]] = v2_counter;
            xxx[xxx_count] = nx;
            yyy[yyy_count] = ny;
            xxx_count++;
            yyy_count++;
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
      rep(i, xxx_count)
      {
        x = xxx[i];
        y = yyy[i];
        visited2[f[x][y]] = -1;
      }
    }

    if (tmp + tmp2 > kouho_scores[num]) {
      kouho_scores[num] = tmp + tmp2;
      kouho_paths[num] = t1 + t2 + t3;
      if (tmp + tmp2 > best_score) {
        best_score = tmp + tmp2;
        best_path = t1 + t2 + t3;
      }
    }

    if (get_elapsed_time() > 1.9)
    {
      break;
    }
    if (!is_sorted && get_elapsed_time() > 1.5)
    {
      sort_kouho(10);
      is_sorted = true;
    }
  }

  cerr << "loop1 = " << loop1 << endl;
  cerr << "loop2 = " << loop2 << endl;
  cerr << "loop3 = " << loop3 << endl;
  cerr << "best_score = " << best_score << endl;

  cout << best_path << endl;
  return 0;
}
