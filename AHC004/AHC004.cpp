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

int pattern_count;
string s[max_patterns];
vector<int> patterns[max_patterns];
int pattern_length[max_patterns] = {};
int grid[n][n];
ll best_score;

bool matched_flag[max_patterns][n][n][2];
int matched_counts[max_patterns];

ll calc_score()
{
  int matched_pattern_count = 0;

  rep(k, pattern_count) {
    bool found = false;
    rep(i, n) {
      rep(j, n) {
        // 横方向チェック
        bool is_match = true;
        rep(l, pattern_length[k]) {
          if (grid[i][(j + l) % n] != patterns[k][l]) {
            is_match = false;
            break;
          }
        }
        if (is_match) {
          found = true;
          break;
        }
        // 縦方向チェック
        is_match = true;
        rep(l, pattern_length[k]) {
          if (grid[(i + l) % n][j] != patterns[k][l]) {
            is_match = false;
            break;
          }
        }
        if (is_match) {
          found = true;
          break;
        }
      }
      if (found) break;
    }
    if (found) matched_pattern_count++;
  }

  ll result = 0;
  if (matched_pattern_count < pattern_count) {
    result = perfect_score * matched_pattern_count / pattern_count;
  }
  else {
    int empty_count = 0;
    rep(i, n)
    {
      rep(j, n)
      {
        if (grid[i][j] == 0) {
          empty_count++;
        }
      }
    }
    result = perfect_score * 2LL * n * n / (2LL * n * n - empty_count);
  }
  return result;
}

ll update_local_score(int x, int y)
{
  rep(i, pattern_count)
  {
    rep(j, n)
    {
      matched_counts[i] -= matched_flag[i][j][y][0];
      matched_counts[i] -= matched_flag[i][j][y][1];
      matched_flag[i][j][y][0] = 0;
      matched_flag[i][j][y][1] = 0;
    }
    rep(k, n)
    {
      matched_counts[i] -= matched_flag[i][x][k][0];
      matched_counts[i] -= matched_flag[i][x][k][1];
      matched_flag[i][x][k][0] = 0;
      matched_flag[i][x][k][1] = 0;
    }
  }

  rep(k, pattern_count)
  {
    int ok = 0;
    srep(i, x, x + 1)
    {
      rep(j, n)
      {
        ok = 1;
        rep(l, pattern_length[k])
        {
          if (grid[i][(j + l) % n] != patterns[k][l]) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          matched_counts[k]++;
          matched_flag[k][i][j][0] = 1;
        }
        ok = 1;
        rep(l, pattern_length[k])
        {
          if (grid[(i + l) % n][j] != patterns[k][l]) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          matched_counts[k]++;
          matched_flag[k][i][j][1] = 1;
        }
      }
    }
  }

  rep(k, pattern_count)
  {
    int ok = 0;
    rep(i, n)
    {
      srep(j, y, y + 1)
      {
        if (i == x) continue;
        ok = 1;
        rep(l, pattern_length[k])
        {
          if (grid[i][(j + l) % n] != patterns[k][l]) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          matched_counts[k]++;
          matched_flag[k][i][j][0] = 1;
        }
        ok = 1;
        rep(l, pattern_length[k])
        {
          if (grid[(i + l) % n][j] != patterns[k][l]) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          matched_counts[k]++;
          matched_flag[k][i][j][1] = 1;
        }
      }
    }
  }

  int cnt = 0;
  rep(i, pattern_count) {
    if (matched_counts[i]) {
      cnt++;
    }
  }

  ll res = 0;
  if (cnt < pattern_count) {
    res = perfect_score * cnt / pattern_count;
  }
  else {
    int cnt2 = 0;
    rep(i, n)
    {
      rep(j, n)
      {
        if (grid[i][j] == 0) {
          cnt2++;
        }
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
    cin >> _n >> pattern_count;
    rep(i, pattern_count)
    {
      cin >> s[i];
      pattern_length[i] = s[i].size();
      rep(j, pattern_length[i]) {
        patterns[i].push_back(s[i][j] - 'A' + 1);
      }
    }
  }
  else {  // ファイル入力する
    int _n;
    ifs >> _n >> pattern_count;
    rep(i, pattern_count)
    {
      ifs >> s[i];
      pattern_length[i] = s[i].size();
      rep(j, pattern_length[i]) {
        patterns[i].push_back(s[i][j] - 'A' + 1);
      }
    }
  }

  rep(i, 20)
  {
    rep(j, 20) {
      grid[i][j] = Rand() % 8 + 1;
    }
  }

  clock_t start_clock = clock();

  best_score = calc_score();

  rep(i, n) {
    update_local_score(i, i);
  }

  int iteration_count = 0;
  while (true) {
    iteration_count++;
    int row = Rand() % n;
    int col = Rand() % n;
    int candidate_value = Rand() % 9;

    int old_value = grid[row][col];
    grid[row][col] = candidate_value;

    ll new_score = update_local_score(row, col);
    if (new_score >= best_score) {
      best_score = new_score;
    }
    else {
      grid[row][col] = old_value;
      update_local_score(row, col);
    }

    clock_t end_clock = clock();
    if ((double)(end_clock - start_clock) / CLOCKS_PER_SEC > 2.9) {
      break;
    }
  }

  rep(i, n)
  {
    rep(j, n)
    {
      if (grid[i][j] == 0) {
        cout << '.';
      }
      else {
        cout << (char)(grid[i][j] + 'A' - 1);
      }
    }
    cout << endl;
  }
  return 0;
}
