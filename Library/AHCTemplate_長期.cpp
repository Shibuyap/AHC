#include <algorithm>
#include <array>
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

typedef pair<int, int> P;
typedef long long int ll;

const int dx[4] = { -1, 0, 1, 0 };
const int dy[4] = { 0, -1, 0, 1 };

// タイマー
class Timer
{
private:
  std::chrono::steady_clock::time_point start_time_clock;

public:
  void start()
  {
    start_time_clock = std::chrono::steady_clock::now();
  }

  double get_elapsed_time()
  {
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time_clock;
    return elapsed.count();
  }
};

Timer timer;

// 乱数
namespace
{
  static uint32_t rand_xorshift()
  {
    static uint32_t x = 123456789;
    static uint32_t y = 362436069;
    static uint32_t z = 521288629;
    static uint32_t w = 88675123;
    uint32_t t = x ^ (x << 11);
    x = y;
    y = z;
    z = w;
    w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
    return w;
  }

  static double rand_01()
  {
    return (rand_xorshift() + 0.5) * (1.0 / UINT_MAX);
  }

  static double rand_range(double l, double r)
  {
    return l + (r - l) * rand_01();
  }

  static uint32_t rand_range(uint32_t l, uint32_t r)
  {
    return l + rand_xorshift() % (r - l + 1); // [l, r]
  }

  void shuffle_array(int* arr, int n)
  {
    for (int i = n - 1; i >= 0; i--) {
      int j = rand_xorshift() % (i + 1);
      int swa = arr[i];
      arr[i] = arr[j];
      arr[j] = swa;
    }
  }
}

const double TIME_LIMIT = 1.9;
int exec_mode;

const int n = 40;
int m;
vector<vector<int>> init_grid(n, vector<int>(n, 0));

vector<vector<int>> input_data(int case_num)
{
  std::ostringstream oss;
  oss << "./in/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
  ifstream ifs(oss.str());

  int _n;
  if (!ifs.is_open()) {
    // 標準入力
    cin >> _n >> m;
    for (int i = 0; i < _n; ++i) {
      string s;
      cin >> s;
      for (int j = 0; j < _n; ++j) {
        init_grid[i][j] = (s[j] == '#') ? 1 : 0;
      }
    }
  }
  else {
    // ファイル入力
    ifs >> _n >> m;
    for (int i = 0; i < _n; ++i) {
      string s;
      ifs >> s;
      for (int j = 0; j < _n; ++j) {
        init_grid[i][j] = (s[j] == '#') ? 1 : 0;
      }
    }
    ifs.close();
  }

  auto grid = init_grid;
  return grid;
}

void output_data(int case_num, vector<P>& ans)
{
  if (exec_mode == 0) {
    // 標準出力
    for (auto p : ans) {
      cout << p.first << ' ' << p.second << endl;
    }
  }
  else {
    // ファイル出力
    std::ostringstream oss;
    oss << "./out/" << std::setw(4) << std::setfill('0') << case_num << ".txt";
    ofstream ofs(oss.str());

    for (auto p : ans) {
      ofs << p.first << ' ' << p.second << endl;
    }

    if (ofs.is_open()) {
      ofs.close();
    }
  }
}

ll calculate_score()
{
  ll res = 0;
  return res;
}

void Move(const vector<vector<int>>& g, const vector<vector<double>>& cnt, vector<vector<double>>& cnt3, int dir)
{
  // 左→右
  if (dir == 0) {
    for (int i = 0; i < n; i++) {
      double now = 0;
      for (int j = 0; j < n; j++) {
        now += cnt[i][j];
        if (j == n - 1) {
          cnt3[i][j] = now;
          now = 0;
        }
        else if (g[i][j + 1] == 1) {
          cnt3[i][j] = now;
          now = 0;
        }
      }
    }
  }
  // 右→左
  if (dir == 1) {
    for (int i = 0; i < n; i++) {
      double now = 0;
      for (int j = n - 1; j >= 0; j--) {
        now += cnt[i][j];
        if (j == 0) {
          cnt3[i][j] = now;
          now = 0;
        }
        else if (g[i][j - 1] == 1) {
          cnt3[i][j] = now;
          now = 0;
        }
      }
    }
  }
  // 左→右
  if (dir == 2) {
    for (int j = 0; j < n; j++) {
      double now = 0;
      for (int i = 0; i < n; i++) {
        now += cnt[i][j];
        if (i == n - 1) {
          cnt3[i][j] = now;
          now = 0;
        }
        else if (g[i + 1][j] == 1) {
          cnt3[i][j] = now;
          now = 0;
        }
      }
    }
  }
  // 左→右
  if (dir == 3) {
    for (int j = 0; j < n; j++) {
      double now = 0;
      for (int i = n - 1; i >= 0; i--) {
        now += cnt[i][j];
        if (i == 0) {
          cnt3[i][j] = now;
          now = 0;
        }
        else if (g[i - 1][j] == 1) {
          cnt3[i][j] = now;
          now = 0;
        }
      }
    }
  }
}

int Method1(vector<vector<int>>& g, vector<P>& ans)
{
  double score = 0;
  double survive = n * n - m;

  vector<vector<double>> cnt(n, vector<double>(n, 0));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (g[i][j] == 0) {
        cnt[i][j] = 1;
      }
    }
  }

  vector<vector<double>> cnt2(n, vector<double>(n, 0));
  vector<vector<double>> cnt3(n, vector<double>(n, 0));

  int q = n * n - m;
  for (int _ = 0; _ < q; _++) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        cnt2[i][j] = 0;
      }
    }

    for (int k = 0; k < 4; k++) {
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          cnt3[i][j] = 0;
        }
      }

      Move(g, cnt, cnt3, k);

      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          cnt2[i][j] += cnt3[i][j] / 4.0;
        }
      }
    }

    int x = -1;
    int y = -1;
    double mi1 = 1e9;
    double ma2 = -1;

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (cnt2[i][j] < mi1) {
          x = i;
          y = j;
          mi1 = cnt2[i][j];
          ma2 = 0;
          for (int k = 0; k < 4; k++) {
            int nx = i + dx[k];
            int ny = j + dy[k];
            if (0 <= nx && nx < n && 0 <= ny && ny < n) {
              ma2 += cnt2[nx][ny];
            }
          }
        }
        else if (cnt2[i][j] == mi1) {
          double tmp = 0;
          for (int k = 0; k < 4; k++) {
            int nx = i + dx[k];
            int ny = j + dy[k];
            if (0 <= nx && nx < n && 0 <= ny && ny < n) {
              tmp += cnt2[nx][ny];
            }
          }
          if (tmp > ma2) {
            x = i;
            y = j;
            mi1 = cnt2[i][j];
            ma2 = tmp;
          }
        }
      }
    }
    cout << x << ' ' << y << endl;
    ans.push_back(P(x, y));
    survive -= cnt2[x][y];
    score += survive;
    cnt2[x][y] = 0;
    g[x][y] = 1;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        cnt[i][j] = cnt2[i][j];
      }
    }
  }

  score /= n * n - m;
  int scoreInt = round(1e6 * score / (n * n - m - 1));
  return scoreInt;
}

ll solve_case(int case_num)
{
  timer.start();

  auto g = input_data(case_num);
  vector<P> ans;

  ll score = Method1(g, ans);

  output_data(case_num, ans);

  return score;
}

int main_new()
{
  exec_mode = 2;

  if (exec_mode == 0) {
    solve_case(0);
  }
  else if (exec_mode <= 2) {
    ll sum_score = 0;
    for (int i = 0; i < 15; i++) {
      ll score = solve_case(i);
      sum_score += score;
      if (exec_mode == 1) {
        cerr << score << endl;
      }
      else {
        cerr << "case = " << setw(2) << i << ", "
          << "score = " << setw(4) << score << ", "
          << "sum = " << setw(5) << sum_score << ", "
          << "time = " << setw(5) << timer.get_elapsed_time() << ", "
          << endl;
      }
    }
  }

  return 0;
}
